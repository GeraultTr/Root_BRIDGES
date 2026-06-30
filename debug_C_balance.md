# C-balance assertion debug (`root_CN.py` line 616-619)

## What the assertion checks

`check_balance` (a `@totalstate` called each step) verifies:

```
actual_C  ≡  previous_C + dt × total_boundary_rate + incremental_deficit
```

where

| symbol | meaning |
|---|---|
| `actual_C` | `compute_root_system_C_content().sum()` – sum of (6×hexose + 12×sucrose + 6×reserve + 5×65×storage_protein + 5×AA + 5×phloem_AA + 5×xylem_AA) × lsm |
| `total_boundary_rate` | sum of all C-carrying fluxes that cross the root-system boundary (mol C s⁻¹, positive = inflow) |
| `incremental_deficit` | correction for pools that were clipped to zero |

---

## Bug 1 — CRITICAL: wrong deficit correction (root_CN.py ~line 550)

### What the code does

```python
current_deficit_amount  = dt * clipped_deficit_rate      # mol C clipped this step
incremental_deficit     = current_deficit_amount - self.previous_deficit_amount
...
expected_C = self.previous_C_amount_in_the_root_system + dt * total_boundary_rate + incremental_deficit
```

`incremental_deficit` is the **change** in clipped amount relative to the previous step.

### Why it is wrong

When a pool is clipped to zero, the pool state that `previous_C` was built from **already reflects the clipping** (the pool is at 0, not at –D). So the starting point for the current step already incorporates last step's deficit.

The correct balance identity for a single step is:

```
actual_C = previous_C + dt × boundary_rate + dt × D_current
         = previous_C + dt × boundary_rate + current_deficit_amount
```

where `D_current` is the C clipped *this* step (pool would have gone negative by this amount without clamping). The previous step's deficit `D_previous` must **not** be subtracted, because it was already baked into `previous_C`.

Substituting the buggy formula:

```
expected_C  = previous_C + dt × boundary_rate + (D_current − D_previous)
residual    = actual_C − expected_C = D_previous × dt
```

The residual equals `previous_deficit_amount`. Whenever deficits persist across consecutive steps (steady-state C-starvation in any pool), the residual is non-zero and grows, eventually breaching `rtol=1e-3` or `atol=1e-5`.

### Fix

```python
# WRONG
incremental_deficit = current_deficit_amount - self.previous_deficit_amount

# CORRECT — previous_C already includes last step's clipping
incremental_deficit = current_deficit_amount
```

Equivalently, drop the subtraction and use `current_deficit_amount` directly in the `expected_C` line:

```python
expected_C = (self.previous_C_amount_in_the_root_system
              + dt * total_boundary_rate
              + current_deficit_amount)
```

`self.previous_deficit_amount` can be kept for diagnostic tracking but is no longer needed for the assertion.

---

## Bug 2 — SECONDARY: overwrite of `sucrose_root_to_shoot_phloem` / `AA_root_to_shoot_phloem` (root_nitrogen.py line 1529)

In `axial_transport_N_arrays`, the collar exchange rates are first set correctly:

```python
# line 1517
props["sucrose_root_to_shoot_phloem"][1] = k_collar × (Cv_sol[root] − cv_shoot_sucrose)
# line 1521
props["AA_root_to_shoot_phloem"][1]       = k_collar × (Cv_sol[root] − cv_shoot_amino_acids)
```

Then, when `boundary_from_reached_segments = True` (which happens if `axial_export_water_up_phloem[root] ≥ 0`, i.e., phloem water reverses), line 1529 overwrites both:

```python
props[cfg["solute_flux_to_shoot"]][1] = (n0.sum() + dt*R_total_actual.sum() − n_sol.sum()) / dt
```

After the algebra the overwrite evaluates to **zero** for phloem solutes (because `R_total_actual` already had the collar exchange subtracted at line 1516). Setting both shoot-exchange registers to 0 causes `check_balance` to see no sucrose/AA exchange at the collar even though the solver applied a non-zero collar flux.

Under normal downward phloem flow `axial_export_water_up_phloem[root] < 0`, so this branch is not entered and the bug is dormant. It activates during transient phloem water reversal (e.g., early morning shoot pressure rise).

**Fix**: guard the line-1529 overwrite so it does not execute for `C_sucrose_root` or `phloem_AA`:

```python
if boundary_from_reached_segments:
    if name not in ("C_sucrose_root", "phloem_AA"):   # ← add this guard
        props[cfg["solute_flux_to_shoot"]][1] = (n0.sum() + dt*R_total_actual.sum() − n_sol.sum()) / dt
```

---

## Boundary term audit (check_balance total_boundary_rate)

All terms verified against `_C_hexose_root`, `_AA`, `_C_hexose_reserve`, `_storage_protein`, and the phloem/xylem radial/axial fluxes.

| flux | C weight | direction | included? |
|---|---|---|---|
| `maintenance_respiration` | −1 mol C / mol C | outflow (CO₂) | ✓ `s1_maint_resp` |
| `N_metabolic_respiration` | −1 mol C / mol C | outflow (CO₂) | ✓ `s1_Nresp` |
| `hexose_exudation` | −6 | outflow to soil | ✓ `s1_hex_exud` |
| `phloem_hexose_exudation` | −6 | outflow to soil | ✓ `s1_ph_hex_exud` (§1 categorisation but correct for total) |
| `hexose_uptake_from_soil` | +6 | inflow from soil | ✓ `s1_hex_uptake` |
| `phloem_hexose_uptake_from_soil` | +6 | inflow from soil | ✓ `s1_ph_hex_uptk` |
| `mucilage_secretion` | −6 | outflow | ✓ `s1_mucilage` |
| `cells_release` | −6 | outflow | ✓ `s1_cells` |
| `hexose_consumption_by_growth` | −6 | → struct mass (untracked) | ✓ `s1_hex_growth` |
| `hexose_consumption_by_fungus` | −6 | outflow to fungus | ✓ `s1_hex_fungus` |
| `amino_acids_consumption_by_growth` | −5 | → struct mass (untracked) | ✓ `s1_AA_growth` |
| `import_AA` | +5 | inflow from soil | ✓ `s1_import_AA` |
| `diffusion_AA_soil` | −5 | outflow to soil | ✓ `s1_diff_AA_sl` |
| `apoplastic_AA_soil_xylem` | −5 | outflow to soil | ✓ `s1_aplastic_AA` |
| `sucrose_root_to_shoot_phloem` | ×−12 net | collar net exchange | ✓ `s2_suc_shoot` |
| `AA_root_to_shoot_phloem` | ×−5 net | collar net exchange | ✓ `s2_phAA_shoot` |
| `AA_root_to_shoot_xylem` | ×−5 net | collar net export | ✓ `s2_xyAA_shoot` |
| `AA_synthesis` / `AA_catabolism` | 0 net | internal (hex↔AA, r=5/6) | ✓ absent (cancel) |
| `sucrose_loading_in_phloem` | 0 net | internal (hex→phloem) | ✓ absent (cancel) |
| `hexose_diffusion_from_phloem` | 0 net | internal (phloem→hex) | ✓ absent (cancel) |
| `hexose_mobilization/immobilization_as_reserve` | 0 net | internal | ✓ absent |
| `storage_synthesis / catabolism` | 0 net | internal (AA↔storage) | ✓ absent |
| `import_Nm` / `export_Nm` | 0 | N only, no C | ✓ absent |

No missing or double-counted boundary terms detected.

---

## How to enable the diagnostic print block

Line 566 has an unconditional `and False` that silences the per-step diagnostic even when a gap is detected:

```python
if gap_detected and False:   # ← False suppresses all output
```

Remove `and False` during debugging so that the section-by-section residuals are printed whenever the gap fires.

---

---

## Bug 3 — CRITICAL: new root segments initialised with wrong CN concentrations (root_growth.py `adding_a_child`)

### What happens

When a lateral root emerges, `adding_a_child` is called with `identical_properties=True`. The function copies `C_hexose_root` and `AA` from the mother element but omits the CN-specific concentrations introduced by `RootCNUnified`:

| property | declared default | typical mature value |
|---|---|---|
| `C_sucrose_root` | 1e-5 mol/g | ~1e-3 mol/g |
| `C_hexose_reserve` | 2e-3 mol/g | ~0.05 mol/g |
| `phloem_AA` | 2e-3 mol/g | ~0.05 mol/g |
| `xylem_AA` | 5e-7 mol/g | ~1e-4 mol/g |
| `storage_protein` | 0 mol/g | ~1e-3 mol/g |

Each omitted property falls back to its `declare(default=...)` value instead of the mother's actual concentration.

### Why the assertion fails

For the step when new segments appear, `check_balance` computes `actual_C` over **all** segments, including the newly created ones. `expected_C` was built from `previous_C` (which covered only the segments that existed last step) plus `dt × boundary_rate`. The sudden concentration mismatch for new segments is:

```
residual ≈ Σ_new_segments  (default_concentration − actual_mother_concentration) × lsm
```

For ~80 new segments with `lsm ≈ 1e-6 g`, the dominant term is `C_hexose_reserve` and `phloem_AA`:

```
residual ≈ 80 × (2e-3 − 0.05) × 1e-6 × 6 C/mol ≈ −7.5e-6 mol C   (per property)
total residual ≈ −7.5e-5 mol C
```

This matches the observed assertion error:

```
AssertionError: Quantities mismatch, expected 0.00125 and actual 0.00118, residual -7.39e-05
```

The assertion tolerance is `atol=1e-5 + rtol×actual_C ≈ 1.1e-5`, which is 7× smaller than the residual.

A TODO comment at `root_growth.py` lines 2524-2525 explicitly flagged this gap: *"FOR TRISTAN: When working with a dynamic root structure, you will need to specify in this function 'ADDING_A_CHILD' your new variables that will either be set to 0 (nil properties) or be equal to that of the mother element."*

### Fix

In the `nil_properties=True` branch (after `AA=0.`):

```python
C_sucrose_root=0.,
C_hexose_reserve=0.,
phloem_AA=0.,
xylem_AA=0.,
storage_protein=0.,
```

In the `identical_properties=True` branch (after `AA=mother_element.AA`):

```python
C_sucrose_root=mother_element.C_sucrose_root,
C_hexose_reserve=mother_element.C_hexose_reserve,
phloem_AA=mother_element.phloem_AA,
xylem_AA=mother_element.xylem_AA,
storage_protein=mother_element.storage_protein,
```

---

## Summary

| priority | location | description | fix |
|---|---|---|---|
| **Critical** | `root_CN.py` ~550 | `incremental_deficit` subtracts previous step's deficit, making `expected_C` too low when deficits persist | use `current_deficit_amount` directly |
| **Critical** | `root_growth.py` `adding_a_child` | CN-specific concentrations not copied to new lateral-root segments; fall back to defaults, causing ~7.5e-5 mol C jump at emergence steps | copy `C_sucrose_root`, `C_hexose_reserve`, `phloem_AA`, `xylem_AA`, `storage_protein` from mother element |
| Secondary | `root_nitrogen.py` 1529 | overwrites collar exchange with 0 during phloem water reversal | guard with `if name not in ("C_sucrose_root", "phloem_AA")` |
| Minor | `root_CN.py` 566 | `and False` silences gap diagnostic | remove `and False` |
