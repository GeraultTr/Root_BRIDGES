# Public packages
import numpy as np
from dataclasses import dataclass

# Utility packages
from metafspm.component_factory import *
from metafspm.component import declare

# Model packages
from rhizodep.soil_model import RhizoInputsSoilModel


family = "soil"


@dataclass
class SoilModel(RhizoInputsSoilModel):
    """
    Empty doc
    """

    # INPUTS
    # FROM NITROGEN MODEL
    mineralN_uptake: float = declare(default=0., unit="mol.s-1", unit_comment="of nitrates", description="", 
                                                    min_value="", max_value="", value_comment="", references="", DOI="", 
                                                    variable_type="input", by="model_nitrogen", state_variable_type="extensive", edit_by="user")
    amino_acids_uptake: float = declare(default=0., unit="mol.s-1", unit_comment="of amino acids", description="", 
                                                    min_value="", max_value="", value_comment="", references="", DOI="",
                                                    variable_type="input", by="model_nitrogen", state_variable_type="extensive", edit_by="user")
    mineralN_diffusion_from_roots: float =  declare(default=0., unit="mol.s-1", unit_comment="of nitrates", description="", 
                                                    min_value="", max_value="", value_comment="", references="", DOI="",
                                                    variable_type="input", by="model_nitrogen", state_variable_type="extensive", edit_by="user")
    amino_acids_diffusion_from_roots: float =  declare(default=0., unit="mol.s-1", unit_comment="of amino acids", 
                                                    min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                                    variable_type="input", by="model_nitrogen", state_variable_type="extensive", edit_by="user")
    mineralN_diffusion_from_xylem: float =  declare(default=0., unit="mol.s-1", unit_comment="of nitrates", 
                                                    min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                                    variable_type="input", by="model_nitrogen", state_variable_type="extensive", edit_by="user")
    amino_acids_diffusion_from_xylem: float =  declare(default=0., unit="mol.s-1", unit_comment="of amino_acids", 
                                                    min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                                    variable_type="input", by="model_nitrogen", state_variable_type="extensive", edit_by="user")
    
    # FROM WATER MODEL
    water_uptake: float =  declare(default=0., unit="mol.s-1", unit_comment="of water", 
                                                    min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                                    variable_type="input", by="model_water", state_variable_type="extensive", edit_by="user")
    
    # FROM METEO
    water_irrigation: float =  declare(default=10/(24*3600), unit="g.s-1", unit_comment="of water", 
                                                    min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                                    variable_type="plant_scale_state", by="meteo", state_variable_type="extensive", edit_by="user")
    water_evaporation: float =  declare(default=5/(24*3600), unit="g.s-1", unit_comment="of water", 
                                                    min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                                    variable_type="plant_scale_state", by="meteo", state_variable_type="extensive", edit_by="user")
    water_drainage: float =  declare(default=5/(24*3600), unit="g.s-1", unit_comment="of water", 
                                                    min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                                    variable_type="plant_scale_state", by="meteo", state_variable_type="extensive", edit_by="user")
    mineral_N_fertilization: float =  declare(default=0., unit="g.s-1", unit_comment="of water", 
                                                    min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                                    variable_type="state_variable", by="meteo", state_variable_type="extensive", edit_by="user")

    # STATE VARIABLES

    # STATES INITIALIZATION
    # C related
    POC: float = declare(default=2.e-3, unit="adim", unit_comment="gC per g of dry soil", description="Particulate Organic Carbon massic concentration in soil", 
                                        value_comment="", references="Fischer et al. 1966", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="intensive", edit_by="user")
    MAOC: float = declare(default=8.e-3, unit="adim", unit_comment="gC per g of dry soil", description="Mineral Associated Organic Carbon in soil", 
                                        value_comment="", references="Fischer et al. 1966", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="intensive", edit_by="user")
    DOC: float = declare(default=2e-7, unit="adim", unit_comment="gC per g of dry soil", description="Dissolved Organic Carbon massic concentration in soil", 
                                        value_comment="", references="Fischer et al. 1966", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="intensive", edit_by="user")
    microbial_C: float = declare(default=0.2e-3, unit="adim", unit_comment="gC per g of dry soil", description="microbial Carbon massic concentration in soil", 
                                        value_comment="", references="Fischer et al. 1966", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="intensive", edit_by="user")
    CO2: float = declare(default=0, unit="adim", unit_comment="gC per g of dry soil", description="Carbon dioxyde massic concentration in soil", 
                                        value_comment="", references="Fischer et al. 1966", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="intensive", edit_by="user")

    # N related
    PON: float = declare(default=0.1e-3, unit="adim", unit_comment="gN per g of dry soil", description="Particulate Organic Nitrogen massic concentration in soil", 
                                        value_comment="", references="Fischer et al. 1966", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="intensive", edit_by="user")
    MAON: float = declare(default=0.8e-3, unit="adim", unit_comment="gN per g of dry soil", description="Mineral-Associated Organic Nitrogen massic concentration in soil", 
                                        value_comment="", references="Fischer et al. 1966", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="intensive", edit_by="user")
    DON: float = declare(default=2e-8, unit="adim", unit_comment="gN per g of dry soil", description="Dissolved Organic Nitrogen massic concentration in soil", 
                                        value_comment="", references="Fischer et al. 1966", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="intensive", edit_by="user")
    microbial_N: float = declare(default=0.03e-3, unit="adim", unit_comment="gN per g of dry soil", description="microbial N massic concentration in soil", 
                                        value_comment="", references="Fischer et al. 1966", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="intensive", edit_by="user")
    dissolved_mineral_N: float = declare(default=20e-6, unit="adim", unit_comment="gN per g of dry soil", description="dissolved mineral N massic concentration in soil",
                                        value_comment="", references="Fischer et al. 1966", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="intensive", edit_by="user")

    C_mineralN_soil: float = declare(default=2.2, unit="mol.m-3", unit_comment="of equivalent mineral nitrogen", description="Mineral nitrogen concentration in soil", 
                                        value_comment="", references="Fischer et al. 1966", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="intensive", edit_by="user")
    C_amino_acids_soil: float = declare(default=8.2e-3, unit="mol.m-3", unit_comment="of equivalent mineral nitrogen", description="Mineral nitrogen concentration in soil", 
                                        value_comment="", references="Fischer et al 2007, water leaching estimation", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="intensive", edit_by="user")

    # Water related
    water_potential_soil: float = declare(default=-0.1e6, unit="Pa", unit_comment="", description="Mean soil water potential", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="intensive", edit_by="user")
    soil_moisture: float = declare(default=0.3, unit="adim", unit_comment="g.g-1", description="Volumetric proportion of water per volume of soil", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="intensive", edit_by="user")
    water_volume: float = declare(default=0.25e-6, unit="m3", unit_comment="", description="Volume of the water in the soil element in contact with a the root segment", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="extensive", edit_by="user")
    
    
    # Structure related
    voxel_volume: float = declare(default=1e-6, unit="m3", unit_comment="", description="Volume of the soil element in contact with a the root segment",
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="extensive", edit_by="user")
    bulk_density: float = declare(default=1.42, unit="g.mL", unit_comment="", description="Volumic density of the dry soil", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="extensive", edit_by="user")
    dry_soil_mass: float = declare(default=1.42, unit="g", unit_comment="", description="dry weight of the considered voxel element", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="extensive", edit_by="user")


    # RATES INITIALIZATION
    # In-voxel rates
    microbial_activity: float = declare(default=0., unit="adim", unit_comment="", description="microbial degradation activity indicator depending on microbial activity locally", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="extensive", edit_by="user")
    degradation_POC: float = declare(default=0., unit=".s-1", unit_comment="gC per g of soil per second", description="degradation rate of POC", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="extensive", edit_by="user")
    degradation_MAOC: float = declare(default=0., unit=".s-1", unit_comment="gC per g of soil per second", description="degradation rate of MAOC", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="extensive", edit_by="user")
    degradation_DOC: float = declare(default=0., unit=".s-1", unit_comment="gC per g of soil per second", description="degradation rate of DOC", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="extensive", edit_by="user")
    degradation_microbial_OC: float = declare(default=0., unit=".s-1", unit_comment="gC per g of soil per second", description="degradation rate of microbial OC", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="extensive", edit_by="user")
    mineral_N_net_mineralization: float = declare(default=0., unit=".s-1", unit_comment="gN per g of soil per second", description="mineral N uptake by micro organisms", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="extensive", edit_by="user")
    
    # Transport rates
    soil_water_flux: float = declare(default=0., unit="m.s-1", unit_comment="m3.m-2.s-1", description="volumetric water flux per surface area from Richards equations", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="extensive", edit_by="user")
    mineral_N_transport: float = declare(default=0., unit="mol.s-1", unit_comment="", description="mineral N advection-dispersion flux derived from Richards water transport", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="extensive", edit_by="user")
    amino_acid_transport: float = declare(default=0., unit="mol.s-1", unit_comment="", description="mineral N advection-dispersion flux derived from Richards water transport", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="extensive", edit_by="user")

    # PARAMETERS

    # C related
    k_POC: float = declare(default=0.5 / (3600 * 24 * 365), unit=".s-1", unit_comment="", description="", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    k_MAOC: float = declare(default=0.01 / (3600 * 24 * 365), unit=".s-1", unit_comment="", description="", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    k_DOC: float = declare(default=20. / (3600 * 24 * 365), unit=".s-1", unit_comment="", description="", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    k_MbOC: float = declare(default=10. / (3600 * 24 * 365), unit=".s-1", unit_comment="", description="", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    
    microbial_C_min: float = declare(default=0.1, unit="adim", unit_comment="gC per g of dry soil", description="", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    microbial_C_max: float = declare(default=0.5, unit="adim", unit_comment="gC per g of dry soil", description="", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    microbial_proportion_of_MAOM: float = declare(default=0.1, unit="adim", unit_comment="gC POM per gC of microbial biomass", description="", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    
    # N related
    CN_ratio_POM: float = declare(default=20, unit="adim", unit_comment="gC per gN", description="", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    CN_ratio_MAOM: float = declare(default=10, unit="adim", unit_comment="gC per gN", description="", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    CN_ratio_microbial_biomass: float = declare(default=11, unit="adim", unit_comment="gC per gN", description="", 
                                        value_comment="", references="Perveen et al. 2014", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    CN_ratio_root_cells: float = declare(default=8, unit="adim", unit_comment="gC per gN", description="", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    
    CUE_POC: float = declare(default=0.2, unit="adim", unit_comment="gC per gC", description="Carbon Use efficiency of microorganism degradation for POC", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    CUE_MAOC: float = declare(default=0.4, unit="adim", unit_comment="gC per gC", description="Carbon Use efficiency of microorganism degradation for MAOC", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    CUE_DOC: float = declare(default=0.4, unit="adim", unit_comment="gC per gC", description="Carbon Use efficiency of microorganism degradation for DOC", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    CUE_MbOC: float = declare(default=0.3, unit="adim", unit_comment="gC per gC", description="Carbon Use efficiency of microorganism degradation for MBC", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    
    # max_N_uptake_per_microbial_C: float = declare(default=1e-7, unit=".s-1", unit_comment="gN per second per gC of microbial C", description="", 
    #                                     value_comment="", references="", DOI="",
    #                                    min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    # Km_microbial_N_uptake: float = declare(default=0.01 / 1e3, unit="adim", unit_comment="gN per g of dry soil", description="", 
    #                                     value_comment="", references="", DOI="",
    #                                    min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    
    ratio_C_per_amino_acid: float = declare(default=6, unit="adim", unit_comment="number of carbon per molecule of amino acid", description="", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    CN_ratio_amino_acids: float = declare(default=6/1.4, unit="adim", unit_comment="", description="CN ratio of amino acids (6 C and 1.4 N on average)", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")

    # Temperature related parameters
    microbial_degradation_rate_max_T_ref: float = declare(default=20, unit="°C", unit_comment="", description="the reference temperature", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    microbial_degradation_rate_max_A: float = declare(default=0., unit="adim", unit_comment="", description="parameter A (may be equivalent to the coefficient of linear increase)", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    microbial_degradation_rate_max_B: float = declare(default=3.98, unit="adim", unit_comment="", description="parameter B (may be equivalent to the Q10 value)", 
                                        value_comment="", references="The value for B (Q10) has been fitted from the evolution of Vmax measured by Coody et al. (1986, SBB), who provided the evolution of the maximal uptake of glucose by soil microorganisms at 4, 12 and 25 degree C.", DOI="",
                                       min_value="", max_value="", variable_type="parametyer", by="model_soil", state_variable_type="", edit_by="user")
    microbial_degradation_rate_max_C: float = declare(default=1, unit="adim", unit_comment="", description="parameter C (either 0 or 1)", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")

    # Water-related parameters
    water_volumic_mass: float = declare(default=1e6, unit="g.m-3", unit_comment="", description="Constant water volumic mass", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="extensive", edit_by="user")
    g_acceleration: float = declare(default=9.806, unit="m.s-2", unit_comment="", description="gravitationnal acceleration constant", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="extensive", edit_by="user")
    saturated_hydraulic_conductivity: float = declare(default=1e-4 / (24*3600), unit="adim", unit_comment="m.s-1", description="staturated hydraulic conductivity parameter", 
                                        value_comment="", references="clay loam estimated with Hydrus, bulk density = 1.42", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    theta_R: float = declare(default=0.0835, unit="adim", unit_comment="m3.m-3", description="Soil retention moisture", 
                                        value_comment="", references="clay loam estimated with Hydrus, bulk density = 1.42", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    theta_S: float = declare(default=0.4383, unit="adim", unit_comment="m3.m-3", description="Soil saturation moisture", 
                                        value_comment="", references="clay loam estimated with Hydrus, bulk density = 1.42", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    water_alpha: float = declare(default=0.0138, unit="cm-3", unit_comment="", description="alpha is the inverse of the air-entry value (or bubbling pressure)", 
                                        value_comment="", references="clay loam estimated with Hydrus, bulk density = 1.42", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    water_n: float = declare(default=1.3945, unit="cm-3", unit_comment="", description="alpha is the inverse of the air-entry value (or bubbling pressure)", 
                                        value_comment="", references="clay loam estimated with Hydrus, bulk density = 1.42", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    field_capacity: float = declare(default=0.36, unit="adim", unit_comment="", description="Soil moisture at which soil doesn't retain water anymore.", 
                                        value_comment="", references="Cornell university, case of sandy loam soil", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    permanent_wilting_point: float = declare(default=0.065, unit="adim", unit_comment="", description="Soil moisture at which soil doesn't retain water anymore.", 
                                        value_comment="", references="Cornell university, case of sandy loam soil", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    water_dt: float = declare(default=3600, unit="s", unit_comment="", description="Initialized time_step to try converging the soil water potential profile", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    min_water_dt: float = declare(default=10, unit="s", unit_comment="", description="min value for adaptative time-step in the convergence cycle for water potential profile", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    max_water_dt: float = declare(default=3600, unit="s", unit_comment="", description="max value for adaptative time-step in the convergence cycle for water potential profile", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    water_potential_tolerance: float = declare(default=1, unit="Pa", unit_comment="", description="tolerance for soil water potential gradient profile convergence", 
                                        value_comment="estimated from general usual for pressure head expression (1e-4 m) * rho * g_acceleration = 0.91 Pa", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    max_iterations: int = declare(default=20, unit="adim", unit_comment="", description="Maximal convergence cycle for water potential profile", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    
    

    # W Initialization parameters
    water_moisture_patch: float = declare(default=0.2, unit="mol.m-3", unit_comment="of equivalent mineral nitrogen", description="Mineral nitrogen concentration in a located patch in soil", 
                                        value_comment="", references="Drew et al. 1975", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    patch_depth_water_moisture: float = declare(default=0., unit="m", unit_comment="", description="Depth of a nitrate patch in soil", 
                                        value_comment="", references="Drew et al. 1975", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    patch_uniform_width_water_moisture: float = declare(default=2*0.1, unit="m", unit_comment="", description="Width of the zone of the patch with uniform concentration of nitrate", 
                                        value_comment="", references="Drew et al. 1975", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    patch_transition_water_moisture: float = declare(default=1e-3, unit="m", unit_comment="", description="Variance of the normal law smooting the boundary transition of a nitrate patch with the background concentration", 
                                        value_comment="", references="Drew et al. 1975", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")


    # N Initialization parameters
    
    dissolved_mineral_N_patch: float = declare(default=20e-6, unit="mol.g-1", unit_comment="of equivalent mineral nitrogen", description="Mineral nitrogen concentration in a located patch in soil", 
                                        value_comment="", references="Drew et al. 1975", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    patch_depth_mineralN: float = declare(default=10e-2, unit="m", unit_comment="", description="Depth of a nitrate patch in soil", 
                                        value_comment="", references="Drew et al. 1975", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    patch_uniform_width_mineralN: float = declare(default=4e-2, unit="m", unit_comment="", description="Width of the zone of the patch with uniform concentration of nitrate", 
                                        value_comment="", references="Drew et al. 1975", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    patch_transition_mineralN: float = declare(default=1e-3, unit="m", unit_comment="", description="Variance of the normal law smooting the boundary transition of a nitrate patch with the background concentration", 
                                        value_comment="", references="Drew et al. 1975", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")


    def __init__(self, g, time_step, **scenario):
        """Pass to inherited init, necessary with data classes"""
        super().__init__(g, time_step, **scenario)


    # SERVICE FUNCTIONS

    # Just ressource for now
    def initiate_voxel_soil(self):
        """
        Note : not tested for now, just computed to support discussions.
        """
        super().initiate_voxel_soil()

        # Set an heterogeneity uppon the mean background
        # Nitrogen
        self.add_patch_repartition_to_soil(property_name="dissolved_mineral_N", patch_value=self.dissolved_mineral_N_patch, 
                                           z_loc=self.patch_depth_mineralN, 
                                           z_width=self.patch_uniform_width_mineralN, 
                                           z_dev=self.patch_transition_mineralN)
        # Water
        self.add_patch_repartition_to_soil(property_name="soil_moisture", patch_value=self.water_moisture_patch, 
                                           z_loc=self.patch_depth_water_moisture, 
                                           z_width=self.patch_uniform_width_water_moisture, 
                                           z_dev=self.patch_transition_water_moisture)
        
        # Initialize volumic concentrations
        self.voxels["dry_soil_mass"] = 1e6 * self.voxels["voxel_volume"] * self.voxels["bulk_density"]
        self.voxels["water_volume"] = self.voxels["voxel_volume"] * self.voxels["soil_moisture"]
        self.voxels["C_mineralN_soil"] = self.voxels["dissolved_mineral_N"] * self.voxels["dry_soil_mass"] / self.voxels["water_volume"] / 14
        self.voxels["C_amino_acids_soil"] = self.voxels["DON"] * self.voxels["dry_soil_mass"] / self.voxels["water_volume"] / 14
        self.voxels["C_hexose_soil"] = self.voxels["DOC"] * self.voxels["dry_soil_mass"] / self.voxels["water_volume"] / 6 / 12
    
    def add_patch_repartition_to_soil(self, property_name: str, patch_value: float, x_loc=None, y_loc=None, z_loc=None, 
                                                                        x_width=0, y_width=0, z_width=0, 
                                                                        x_dev=1e-3, y_dev=1e-3, z_dev=1e-3,
                                                                        spherical_normal_patch = False, normal_boundaries = False):
        
        if spherical_normal_patch:
            y_dev = x_dev
            z_dev = x_dev
            x_width = 0
            y_width = 0
            z_width = 0
        
        # Start with Z
        if z_loc is not None:
            
            z_mean = (self.voxels["z1"] + self.voxels["z2"]) / 2

            test = np.logical_and(z_loc - z_width/2 < z_mean, z_mean < z_loc + z_width/2)
            self.voxels[property_name][test] = patch_value
            if normal_boundaries:
                test = z_mean > z_loc + z_width/2
                new_values = self.voxels[property_name] + (patch_value - self.voxels[property_name]) / (z_dev * np.sqrt(2 * np.pi)) * np.exp(-((z_mean - (z_loc + z_width/2)) ** 2) / (2 * z_dev ** 2))
                self.voxels[property_name][test] = new_values[test]
                test = z_mean < z_loc - z_width/2
                new_values = self.voxels[property_name] + (patch_value - self.voxels[property_name]) / (z_dev * np.sqrt(2 * np.pi)) * np.exp(-((z_mean - (z_loc - z_width/2)) ** 2) / (2 * z_dev ** 2))
                self.voxels[property_name][test] = new_values[test]

        # Then x and y
        if x_loc is not None:
            x_mean = (self.voxels["x1"] + self.voxels["x2"]) / 2

            self.voxels[property_name][x_loc - x_width/2 < x_mean < x_loc + x_width/2] = patch_value
            self.voxels[property_name][x_mean > x_loc + x_width/2] = self.voxels[property_name] + (patch_value - self.voxels[property_name]) / (x_dev * np.sqrt(2 * np.pi)) * np.exp(-((x_mean - (x_loc + x_width/2)) ** 2) / (2 * x_dev ** 2))
            self.voxels[property_name][x_mean < x_loc - x_width/2] = self.voxels[property_name] + (patch_value - self.voxels[property_name]) / (x_dev * np.sqrt(2 * np.pi)) * np.exp(-((x_mean - (x_loc - x_width/2)) ** 2) / (2 * x_dev ** 2))

        if y_loc is not None:
            y_mean = (self.voxels["y1"] + self.voxels["y2"]) / 2

            self.voxels[property_name][y_loc - y_width/2 < y_mean < y_loc + y_width/2] = patch_value
            self.voxels[property_name][y_mean > y_loc + y_width/2] = self.voxels[property_name] + (patch_value - self.voxels[property_name]) / (y_dev * np.sqrt(2 * np.pi)) * np.exp(-((y_mean - (y_loc + y_width/2)) ** 2) / (2 * y_dev ** 2))
            self.voxels[property_name][y_mean < y_loc - y_width/2] = self.voxels[property_name] + (patch_value - self.voxels[property_name]) / (y_dev * np.sqrt(2 * np.pi)) * np.exp(-((y_mean - (y_loc - y_width/2)) ** 2) / (2 * y_dev ** 2))

    # MODEL EQUATIONS

    # RATES

    @potential
    @rate
    def _microbial_activity(self, microbial_C, soil_temperature):
        temperature_regulation = self.temperature_modification(soil_temperature=soil_temperature,
                                                                   T_ref=self.microbial_degradation_rate_max_T_ref,
                                                                    A=self.microbial_degradation_rate_max_A,
                                                                    B=self.microbial_degradation_rate_max_B,
                                                                    C=self.microbial_degradation_rate_max_C)
        test = (microbial_C > self.microbial_C_min) & (microbial_C < self.microbial_C_max)
        results = np.ones_like(microbial_C)
        results[test] = (1. + (microbial_C[test] - self.microbial_C_min)  / (self.microbial_C_max - self.microbial_C_min))
        return results * temperature_regulation

    @actual
    @rate
    def _degradation_POC(self, microbial_activity, POC):
        return self.k_POC * microbial_activity * POC
    
    @actual
    @rate
    def _degradation_MAOC(self, microbial_activity, MAOC):
        return self.k_MAOC * microbial_activity * MAOC
    
    @actual
    @rate
    def _degradation_DOC(self, microbial_activity, DOC):
        return self.k_DOC * microbial_activity * DOC
    
    @actual
    @rate
    def _degradation_microbial_OC(self, microbial_activity, microbial_C):
        return self.k_MbOC * microbial_activity * microbial_C
    

    def soil_moisture_capacity(self, psi):
        """
        Specific moisture capacity function, C(psi)
        Derivarion of the soil_moisture function
        """
        m = 1 - 1/self.water_n
        return (self.theta_S - self.theta_R) * self.water_alpha * self.water_n * m * (((self.water_alpha * np.abs(psi))**(self.water_n - 1)) /
                                                                                     ((1 + (self.water_alpha * np.abs(psi))**self.water_n)**(m + 1)))

    def soil_water_conductivity(self, theta):
        """
        Compute water conductivity at each point as function of soil moisture according to the van Genuchten-Mualem Model
        """
        m = 1-1/self.water_n
        Se = (theta - self.theta_R) / (self.theta_S - self.theta_R)
        return self.saturated_hydraulic_conductivity * Se**0.5 * (1 - (1 - Se**(1/m))**m)**2
    
    def _soil_moisture(self, water_potential_soil):
        m = 1 - (1/self.water_n)
        return self.theta_R + (self.theta_S - self.theta_R) / (1 + np.abs(self.water_alpha * water_potential_soil)**self.water_n) ** m

    #TP@potential
    #TP@rate
    def richards_1D_water_flux(self):
        """
        Richards_1D_water_flux
        """
        water_potential_soil = self.voxels["water_potential_soil"]
        theta = self.voxels["soil_moisture"]

        volumetric_source_water_flow = - self.voxels["water_uptake"] * 18
        print(self.water_evaporation)
        volumetric_source_water_flow[:, :, 1] += self.water_irrigation[1] - self.water_evaporation[1]
        volumetric_source_water_flow[:, :, -1] -= self.water_drainage[1]
        volumetric_source_water_flow /= self.voxels["dry_soil_mass"] # to convert g of W per g of soil per s-1

        converged = False
        iteration = 0

        while not converged and iteration < self.max_iterations:
            water_potential_soil_prev = water_potential_soil.copy()

            # Update Richards equation and water flux

            # Compute water conductivity at each point as function of soil moisture according to the van Genuchten-Mualem Model
            K = self.soil_water_conductivity(theta)

             # Calculate fluxes between points, including gravitational term
            flux_forward = K[:, :, 1:-1] * ((water_potential_soil[:, :, 2:] - water_potential_soil[:, :, 1:-1]) / self.delta_z + self.water_volumic_mass * self.g_acceleration)
            flux_backward = K[:, :, 1:-1] * ((water_potential_soil[:, :, 1:-1] - water_potential_soil[:, :, :-2]) / self.delta_z + self.water_volumic_mass * self.g_acceleration)

            water_potential_soil[:, :, 1:-1] += self.water_dt * (1 / self.soil_moisture_capacity(water_potential_soil[:, :, 1:-1])) * ((flux_forward - flux_backward) / self.delta_z + volumetric_source_water_flow[:, :, 1:-1])

            # Convergence check
            error = np.max(np.abs(water_potential_soil - water_potential_soil_prev))
            if error < self.water_potential_tolerance:
                converged = True
                
                self.voxels["water_potential_soil"] = water_potential_soil

                # Calculate Darcy flux for use in other equations
                self.voxels["soil_water_flux"] = -K * (np.gradient(water_potential_soil, self.delta_z, axis=2) + self.water_volumic_mass * self.g_acceleration)

                # Update water content
                self.voxels["soil_moisture"] = self._soil_moisture(water_potential_soil)

            else:
                iteration += 1
        
        if iteration < 3:
            self.water_dt = min(self.water_dt * 1.2, self.max_water_dt)  # Increase dt if convergence was fast
        elif iteration == self.max_iterations:
            self.water_dt = max(self.water_dt * 0.5, self.min_water_dt) # Lower it if convergence took too much time, we lack precision there


    def solute_water_advection(self, soil_water_flux, solute_concentration):
        # Advection term: calculate forward and backward flux components
        advection_forward = soil_water_flux[:, :, 1:-1] * (solute_concentration[:, :, 2:] + solute_concentration[:, :, 1:-1]) / 2
        advection_backward = soil_water_flux[:, :, 1:-1] * (solute_concentration[:, :, 1:-1] + solute_concentration[:, :, :-2]) / 2

        # Combine forward and backward components and multiply by voxel cross sectionnal area to get the molar flux
        return (advection_forward - advection_backward) * self.voxels_Z_section_area
    
    def solute_diffusion(self, soil_water_flux, soil_moisture, solute_concentration):
        D_m = 1e-9  # Molecular diffusion coefficient (m^2/s) should differ for each solute
        alpha_L = 0.1  # Longitudinal dispersivity (m)
        dispersion_coefficient = D_m + alpha_L * np.abs(soil_water_flux) / soil_moisture # m^2/s

        # Dispersion term: calculate forward and backward components
        dispersion_forward = dispersion_coefficient[:, :, 1:-1] * (solute_concentration[:, :, 2:] - solute_concentration[:, :, 1:-1]) / self.delta_z
        dispersion_backward = dispersion_coefficient[:, :, 1:-1] * (solute_concentration[:, :, 1:-1] - solute_concentration[:, :, :-2]) / self.delta_z

        # Combine forward and backward components and multiply by voxel cross sectionnal area to get the molar flux
        return (dispersion_forward - dispersion_backward) * self.voxels_Z_section_area

    #TP@actual
    #TP@rate
    def _mineral_N_transport(self, mineral_N_transport, soil_water_flux, soil_moisture, C_mineralN_soil):
        mineral_N_transport[:, :, 1:-1] = self.solute_diffusion(soil_water_flux, soil_moisture, C_mineralN_soil) - self.solute_water_advection(soil_water_flux, C_mineralN_soil)
        return mineral_N_transport
    
    #TP@actual
    #TP@rate
    def _amino_acid_transport(self, amino_acid_transport, soil_water_flux, soil_moisture, C_amino_acids_soil):
        amino_acid_transport[:, :, 1:-1] = self.solute_diffusion(soil_water_flux, soil_moisture, C_amino_acids_soil) - self.solute_water_advection(soil_water_flux, C_amino_acids_soil)
        return amino_acid_transport
    
    #TP@actual
    #TP@rate
    def _mineral_N_fertilization(self, dry_soil_mass):
        result = np.zeros_like(dry_soil_mass)
        result[:, :, 1] = self.mineral_N_fertilization[1]
        return result

    # STATES

    @state
    def _POC(self, POC, dry_soil_mass, degradation_POC, cells_release):
        return POC + (self.time_step_in_seconds / dry_soil_mass) * (
            cells_release
            - degradation_POC * dry_soil_mass
        )
    
    @state
    def _PON(self, POC, dry_soil_mass, degradation_POC, cells_release):
        return POC + (self.time_step_in_seconds / dry_soil_mass) * (
            cells_release / self.CN_ratio_root_cells
            - degradation_POC *dry_soil_mass / self.CN_ratio_POM
        )
    
    @state
    def _MAOC(self, MAOC, dry_soil_mass, degradation_microbial_OC, degradation_MAOC):
        return MAOC + (self.time_step_in_seconds / dry_soil_mass) * (
            degradation_microbial_OC * dry_soil_mass * self.microbial_proportion_of_MAOM
            - degradation_MAOC * dry_soil_mass
        )
    
    @state
    def _MAON(self, MAON, dry_soil_mass, degradation_microbial_OC, degradation_MAOC):
        return MAON + (self.time_step_in_seconds / dry_soil_mass) * (
            degradation_microbial_OC * dry_soil_mass * self.microbial_proportion_of_MAOM / self.CN_ratio_microbial_biomass
            - degradation_MAOC * dry_soil_mass / self.CN_ratio_MAOM
        )

    @state
    def _DOC(self, DOC, dry_soil_mass, degradation_microbial_OC, degradation_DOC, hexose_exudation, phloem_hexose_exudation, mucilage_secretion, amino_acids_diffusion_from_roots,  amino_acids_diffusion_from_xylem, amino_acids_uptake, amino_acid_transport):
        return DOC + (self.time_step_in_seconds / dry_soil_mass) * (
            degradation_microbial_OC * dry_soil_mass * (1 - self.microbial_proportion_of_MAOM)
            - degradation_DOC * dry_soil_mass
            + hexose_exudation
            + phloem_hexose_exudation
            + mucilage_secretion
            + amino_acids_diffusion_from_roots
            + amino_acids_diffusion_from_xylem
            - amino_acids_uptake
            + amino_acid_transport
        )
    
    @state
    def _DON(self, DON, DOC, dry_soil_mass, degradation_microbial_OC, degradation_DOC, amino_acids_diffusion_from_roots,  amino_acids_diffusion_from_xylem, amino_acids_uptake, amino_acid_transport):
        return DON + (self.time_step_in_seconds / dry_soil_mass) * (
            degradation_microbial_OC * dry_soil_mass * (1 - self.microbial_proportion_of_MAOM) / self.CN_ratio_microbial_biomass
            - degradation_DOC * dry_soil_mass * DON / DOC
            + (amino_acids_diffusion_from_roots
            + amino_acids_diffusion_from_xylem 
            - amino_acids_uptake
            + amino_acid_transport) / self.CN_ratio_amino_acids
        )
    
    @state
    def _microbial_C(self, microbial_C, dry_soil_mass, degradation_POC, degradation_MAOC, degradation_DOC, degradation_microbial_OC):
        """
        For microbial biomass C, the new concentrations results i) from the turnover of this pool, ii) from a fraction of
        the degradation products of each SOC pool, including microbial biomass itself, and iii) from new net inputs, if any:
        """
        return microbial_C + (self.time_step_in_seconds / dry_soil_mass) * (
            - degradation_microbial_OC
            + degradation_POC * self.CUE_POC
            + degradation_MAOC * self.CUE_MAOC
            + degradation_DOC * self.CUE_DOC
            + degradation_microbial_OC * self.CUE_MbOC
        )
    
    @state
    def _microbial_N(self, microbial_N, microbial_C, dry_soil_mass, DOC, DON, degradation_POC, degradation_MAOC, degradation_DOC, degradation_microbial_OC):
        """
        For microbial biomass C, the new concentrations results i) from the turnover of this pool, ii) from a fraction of
        the degradation products of each SOC pool, including microbial biomass itself, and iii) from new net inputs, if any:
        """

        balance = microbial_N + (self.time_step_in_seconds / dry_soil_mass) * (
            - degradation_microbial_OC / self.CN_ratio_microbial_biomass
            + degradation_POC * self.CUE_POC / self.CN_ratio_POM
            + degradation_MAOC * self.CUE_MAOC / self.CN_ratio_MAOM
            + degradation_DOC * self.CUE_DOC * DON / DOC
            + degradation_microbial_OC * self.CUE_MbOC/ self.CN_ratio_microbial_biomass
        )

        self.voxels["mineral_N_net_mineralization"] = dry_soil_mass * (balance - microbial_C / self.CN_ratio_microbial_biomass)
        return microbial_C / self.CN_ratio_microbial_biomass
    
    @state
    def _CO2(self, CO2, dry_soil_mass, degradation_POC, degradation_MAOC, degradation_DOC, degradation_microbial_OC):
        return CO2 + (self.time_step_in_seconds / dry_soil_mass) * (
                degradation_POC * (1 - self.CUE_POC)
                + degradation_MAOC * (1 - self.CUE_MAOC)
                + degradation_microbial_OC * (1 - self.CUE_MbOC)
                + degradation_DOC * (1 - self.CUE_DOC)
            )
    
    #TP@state
    def _dissolved_mineral_N(self, dissolved_mineral_N, dry_soil_mass, mineral_N_net_mineralization, 
                             mineralN_diffusion_from_roots, mineralN_diffusion_from_xylem, mineralN_uptake, mineral_N_fertilization, mineral_N_transport):
        return dissolved_mineral_N + (self.time_step_in_seconds / dry_soil_mass) * (
            mineral_N_net_mineralization
            + mineralN_diffusion_from_roots
            + mineralN_diffusion_from_xylem
            - mineralN_uptake
            + mineral_N_fertilization
            + mineral_N_transport
        )

    #TP@segmentation
    #TP@state
    def _C_mineralN_soil(self, dissolved_mineral_N, dry_soil_mass, soil_moisture, voxel_volume):
        return dissolved_mineral_N * (dry_soil_mass / (soil_moisture * voxel_volume)) / 14

    @segmentation
    @state
    def _C_amino_acids_soil(self, DOC, dry_soil_mass, soil_moisture, voxel_volume):
        return DOC * (dry_soil_mass * (soil_moisture / voxel_volume)) / 12 / self.ratio_C_per_amino_acid
    
    @segmentation
    @state
    def _C_hexose_soil(self, DOC, dry_soil_mass, soil_moisture, voxel_volume):
        return DOC * (dry_soil_mass * (soil_moisture * voxel_volume)) / 14 / 6
    
    @state
    def _water_volume(self, soil_moisture, voxel_volume):
        return soil_moisture * voxel_volume
    
    @state
    def _water_potential_soil(self, voxel_volume, water_volume):
        """
        Water retention curve from van Genuchten 1980
        """
        m = 1 - (1/self.water_n)
        return - (1 / self.water_alpha) * (
                                            ((self.theta_S - self.theta_R) / ((water_volume / voxel_volume) - self.theta_R)) ** (1 / m) - 1 
                                        )** (1 / self.water_n)
