from dataclasses import dataclass
import numpy as np
import pandas as pd
from openalea.metafspm.component_factory import *
from openalea.metafspm.component import declare

from openalea.rhizodep import RootCarbonModel
from openalea.rootcynaps import RootNitrogenModel

# Deported class inheritance to include this information in the __globals__, so that it can be picked by decorators to merge the steps of all classes
inheriting = (RootCarbonModel, RootNitrogenModel)


@dataclass
class RootCNUnified(*inheriting):

    # @note INPUTS

    # FROM SHOOT MODEL
    Cv_sucrose_phloem_collar: float = declare(default=950, unit="mol.m-3", unit_comment="", description="Sucrose volumic concentration in phloem at collar point", 
                                       min_value=0, max_value=1200, value_comment="", references="Winter et al. 1992", DOI="",
                                        variable_type="input", by="model_shoot", state_variable_type="", edit_by="user")
    mstruct_axis_shoot: float = declare(default=0.0541, unit="g", unit_comment="", description="Total axis initial structural mass, shoot because it cannot have the same variable name as the shoot model to avoid confusion", 
                                       min_value=0, max_value=1200, value_comment="", references="", DOI="",
                                        variable_type="input", by="model_shoot", state_variable_type="", edit_by="user")
    sucrose_phloem_shoot: float = declare(default=17 / 12 / 1e6, unit="mol", unit_comment="of sucrose", description="", 
                                       min_value=0, max_value=1200, value_comment="", references="", DOI="",
                                        variable_type="input", by="model_shoot", state_variable_type="", edit_by="user")
    AA_phloem_shoot: float = declare(default=1 / 12 / 1e6, unit="mol", unit_comment="of amino acids", description="", 
                                       min_value=0, max_value=1200, value_comment="", references="", DOI="",
                                        variable_type="input", by="model_shoot", state_variable_type="", edit_by="user")
    sucrose_phloem_contributors_flow: float = declare(default=30e-6 / 3600 / 12, unit="mol.s-1", unit_comment="of sucrose", description="", 
                                       min_value=0, max_value=1200, value_comment="", references="", DOI="",
                                        variable_type="input", by="model_shoot", state_variable_type="", edit_by="user")
    AA_phloem_contributors_flow: float = declare(default=1e-6 / 3600 / 1.4, unit="mol.s-1", unit_comment="of amino acids", description="", 
                                       min_value=0, max_value=1200, value_comment="", references="", DOI="",
                                        variable_type="input", by="model_shoot", state_variable_type="", edit_by="user")

    # FROM GROWTH MODEL
    amino_acids_consumption_by_growth: float = declare(default=0., unit="mol.s-1", unit_comment="", description="amino_acids consumption rate by growth processes", 
                                                 min_value="", max_value="", value_comment="", references="", DOI="",
                                                  variable_type="input", by="model_growth", state_variable_type="", edit_by="user")

    # @note STATE VARIABLES

    N_metabolic_respiration: float = declare(default=0., unit="mol.s-1", unit_comment="of carbon", description="Respiration related to nitrogen metabolism", 
                                            min_value="", max_value="", value_comment="", references="", DOI="",
                                             variable_type="state_variable", by="model_cn", state_variable_type="NonInertialExtensive", edit_by="user")
    nitrate_transporters_affinity_factor: float = declare(default=1., unit="mol.s-1", unit_comment="of nitrates", description="nitrate_transporters_affinity_factor, introduced to account for NRT1 signalling function when going through LATS regime", 
                                                    min_value="", max_value="", value_comment="", references="Remans et al 2006", DOI="", 
                                                    variable_type="state_variable", by="model_cn", state_variable_type="NonInertialIntensive", edit_by="user")
    total_hexose_diffusion_from_phloem: float = declare(default=0., unit="umol of C.g-1 mstruc.h-1", unit_comment="", description="Property computed to compare with shoot model unloading",
                                    min_value="", max_value="", value_comment="", references="", DOI="",
                                    variable_type="plant_scale_state", by="model_cn", state_variable_type="", edit_by="user")
    
    # @note SUMMED STATE VARIABLES

    sucrose_root_to_shoot_phloem: float =       declare(default=-1e-10, unit="mol.s-1", unit_comment="of sucrose", description="",
                                                min_value="", max_value="", value_comment="WARNING: value has to be non 0 for proper water flow computation at init", references="", DOI="",
                                                variable_type="plant_scale_state", by="model_cn", state_variable_type="", edit_by="user")
    Cv_sucrose_average: float =                   declare(default=1., unit="mol.m-3", unit_comment="of amino acids", description="", 
                                                min_value="", max_value="", value_comment="", references="", DOI="",
                                                variable_type="plant_scale_state", by="model_cn", state_variable_type="", edit_by="user")
    Cv_hexose_average: float =                   declare(default=1., unit="mol.m-3", unit_comment="of amino acids", description="", 
                                                min_value="", max_value="", value_comment="", references="", DOI="",
                                                variable_type="plant_scale_state", by="model_cn", state_variable_type="", edit_by="user")
    Cv_sucrose_root: float =                   declare(default=1., unit="mol.m-3", unit_comment="of amino acids", description="", 
                                                min_value="", max_value="", value_comment="", references="", DOI="",
                                                variable_type="state_variable", by="model_cn", state_variable_type="", edit_by="user")
    Cv_hexose_root: float =                   declare(default=1., unit="mol.m-3", unit_comment="of amino acids", description="", 
                                                min_value="", max_value="", value_comment="", references="", DOI="",
                                                variable_type="state_variable", by="model_cn", state_variable_type="", edit_by="user")
    
    # @note PARAMETERS
    r_hexose_AA: float = declare(default=5/6, unit="adim", unit_comment="mol of hexose per mol of amino acids in roots", description="stoechiometric ratio during amino acids synthesis for hexose consumption", 
                                min_value="", max_value="", value_comment="", references="we hypothesize from Yemm and Willis 1956 that synthetized soluble amino acids are mainly composed of glutamine and asparagine", DOI="",
                                variable_type="parameter", by="model_carbon", state_variable_type="", edit_by="user")
    r_Nm_AA: float =     declare(default=1.4, unit="adim", unit_comment="mol of N per mol of amino acids", description="concentration stoechiometric ratio between mineral nitrogen and amino acids in roots", 
                                min_value="", max_value="", value_comment="", references="", DOI="",
                                variable_type="parameter", by="model_nitrogen", state_variable_type="", edit_by="user")
    respi_costs_mineralN_import: float = declare(default=0.397, unit="adim", unit_comment="mol of C per mol of N", description="Respiratory of active imports in root.", 
                                min_value="", max_value="", value_comment="", references="Barillot et al., 2016", DOI="",
                                variable_type="parameter", by="model_carbon", state_variable_type="", edit_by="user")
    respi_costs_mineralN_reduction: float = declare(default=1.98, unit="adim", unit_comment="mol of C per mol of N", description="Respiratory of active imports in root.", 
                                min_value="", max_value="", value_comment="", references="Robinson 2001; Barillot et al., 2016", DOI="",
                                variable_type="parameter", by="model_carbon", state_variable_type="", edit_by="user")
    massic_reference_rate_of_AA_consumption_by_growth: float = declare(default=2.78e-10, unit="mol.s-1.g-1", unit_comment="of hexose", description="Coefficient of permeability of unloading phloem", 
                                                min_value="", max_value="", value_comment="From RhizoDep parameter, applied 5e-13 * 6 * 12 / 0.44 * 0.015 / 14 / 1.4", references="Reference consumption rate of hexose for growth for a given root element (used to multiply the reference unloading rate when growth has consumed hexose)", DOI="",
                                                variable_type="parameter", by="model_carbon", state_variable_type="", edit_by="user")
    massic_reference_rate_of_hexose_consumption_by_growth: float = declare(default=1.39e-9, unit="mol.s-1", unit_comment="of hexose", description="Coefficient of permeability of unloading phloem", 
                                                min_value="", max_value="", value_comment="", references="Reference consumption rate of hexose for growth for a given root element (used to multiply the reference unloading rate when growth has consumed hexose)", DOI="",
                                                variable_type="parameter", by="model_carbon", state_variable_type="", edit_by="user")
    
    
    def __init__(self, g, time_step: int,  **scenario: dict):
        """
        DESCRIPTION
        -----------
        __init__ method

        :param g: the root MTG
        :param time_step: time step of the simulation (s)
        :param scenario: mapping of existing variable initialization and parameters to superimpose.
        :return:
        """
        # Before any other operation, we apply the provided scenario by changing default parameters and initialization
        self.apply_scenario(**scenario)

        self.g = g
        self.props = self.g.properties()
        self.time_step = time_step
        self.choregrapher.add_time_and_data(instance=self, sub_time_step=self.time_step, data=self.props)
        self.vertices = self.g.vertices(scale=self.g.max_scale())

        # Before any other operation, we apply the provided scenario by changing default parameters and initialization
        self.link_self_to_mtg()
        self.initiate_heterogeneous_variables()

        # self.total_root_sucrose_and_living_struct_mass() # Needed otherwise first shoot unloading will be unrealistic

        self.solute_configs["C_sucrose_root"] = {
        "solute_massic_concentration_prop": "C_sucrose_root",
        "solute_massic_concentration_symplasm": "C_hexose_root",
        "diffusive_flux_name": "hexose_diffusion_from_phloem",
        "diffusive_flux_conversion": - 1 / 2,
        "diffusion_parameter": "diffusion_phloem",
        "conductive_element_volume_prop": "phloem_volume",
        "water_flux_prop": "axial_export_water_up_phloem",
        "radial_solute_flux": lambda hexose_diffusion_from_phloem, hexose_active_production_from_phloem, phloem_hexose_exudation, sucrose_loading_in_phloem, phloem_hexose_uptake_from_soil : (
                                        - hexose_diffusion_from_phloem / 2.
                                        - hexose_active_production_from_phloem / 2.
                                        - phloem_hexose_exudation / 2.
                                        + sucrose_loading_in_phloem
                                        + phloem_hexose_uptake_from_soil / 2.),
        "flux_shoot_boundary": lambda props: props["sucrose_input_rate"][1],
        "boundary_shoot_solute_concentration": lambda props: props["Cv_sucrose_phloem_collar"][1],
        "solute_flux_to_shoot": "sucrose_root_to_shoot_phloem",
        "solute_volumic_concentration_bounds": (1e-4, 3e3),
        }

        self.cumulated_time = 0.

        # struct_mass_N_content = 0.005 / 14
        # struct_mass_C_content = 0.44 / 12

        # for vid in self.vertices:
        #     ini = self.props["AA"][vid]
        #     self.props["AA"][vid] = - self.props["C_hexose_root"][vid] * 6 / (5 - (1.4 * struct_mass_C_content / struct_mass_N_content))
        #     print(self.props["AA"][vid], ini / self.props["AA"][vid])

    # @note PROCESSES

    # Note, here the decorator naming doesn't make much sense, but it was placed so that resolution of this flux is made after every other one.
    # Indeed, the expected behovior is to have rates computed from previous time step states. However, if we didn't waited for all import / export to compute,
    # This respiration would have reflected states of two time-steps ago.
    @actual
    @rate
    def _N_metabolic_respiration(self, import_Nm, export_Nm, import_AA, export_AA, AA_synthesis):
        """
        Here we explicit a respiration cost associated to the exchange of N in the root (e.g. cost for uptake and xylem loading + anabolism costs of soluble nitrogen).        
        In line with Thornley and Cannell 2000; Robinson 2001 and Barillot et al. 2016, respiratory costs are considered as proportionnal to transport and synthesis costs.
        Moreover, such synthesis can have side effects of carbohydrates catabolism affecting the respiration through the release of labile hexoses (Yemm and Willis 1956).
        """
        transport_respiration = self.respi_costs_mineralN_import * (import_Nm + export_Nm + import_AA + export_AA)
        anabolism_respiration = self.respi_costs_mineralN_reduction * AA_synthesis * self.r_Nm_AA
        return transport_respiration + anabolism_respiration
    

    @rate
    def _hexose_diffusion_from_phloem(self, type, length, label, root_order, phloem_exchange_surface, C_sucrose_root, C_hexose_root,
                                             hexose_consumption_by_growth, deficit_hexose_root, living_struct_mass, symplasmic_volume, phloem_volume, soil_temperature):
        """
        Superimposing original, staying with a massic concentration gradient as fist approximation to avoid changing parameters
        """

        Cv_sucrose_root = C_sucrose_root * living_struct_mass / phloem_volume
        Cv_hexose_root = C_hexose_root * living_struct_mass / symplasmic_volume

        reference_rate_of_hexose_consumption_by_growth = self.reference_rate_of_hexose_consumption_by_growth
        # reference_rate_of_hexose_consumption_by_growth = np.where(label==self.label_Apex, reference_rate_of_hexose_consumption_by_growth/1, reference_rate_of_hexose_consumption_by_growth)

        phloem_permeability = self.diffusion_phloem * (1 + np.where(type == self.type_Base_of_the_root_system, (deficit_hexose_root) / (reference_rate_of_hexose_consumption_by_growth),
                                                                    (hexose_consumption_by_growth + deficit_hexose_root) / (reference_rate_of_hexose_consumption_by_growth)))
        # phloem_permeability = self.diffusion_phloem * (1 + hexose_consumption_by_growth /
        #                                                     (living_struct_mass * self.massic_reference_rate_of_hexose_consumption_by_growth))
        # phloem_permeability = self.diffusion_phloem

        phloem_permeability *= self.temperature_modification(soil_temperature=soil_temperature,
                                                                T_ref=self.phloem_unloading_T_ref,
                                                                A=self.phloem_unloading_A,
                                                                B=self.phloem_unloading_B,
                                                                C=self.phloem_unloading_C)
        
        flux = phloem_permeability * (np.maximum(0, Cv_sucrose_root) - np.maximum(0, Cv_hexose_root / 2.)) * phloem_exchange_surface

        # return np.where(flux > 0., flux, 0.)
        return flux

        # return np.where((length <= 0.) | (type == self.type_Just_dead) | (type == self.type_Dead), 0.,
        #         2. * phloem_permeability * (Cv_sucrose_root - Cv_hexose_root / 2.) * phloem_exchange_surface)



   
    @rate
    def _hexose_active_production_from_phloem(self, type, C_sucrose_root, length, phloem_exchange_surface,
                                              hexose_consumption_by_growth, soil_temperature):
        """
        Superimposing original, staying with a massic concentration gradient as fist approximation to avoid changing parameters
        """
        # Removed condition to limit based on deficit compared to RhizoDep
        max_unloading_rate = self.max_unloading_rate * (1 + hexose_consumption_by_growth /
                                                        self.reference_rate_of_hexose_consumption_by_growth)
        max_unloading_rate *= self.temperature_modification(soil_temperature=soil_temperature,
                                                            T_ref=self.phloem_unloading_T_ref,
                                                            A=self.phloem_unloading_A,
                                                            B=self.phloem_unloading_B,
                                                            C=self.phloem_unloading_C)
        
        return np.where((length <= 0.) | (phloem_exchange_surface <= 0.) | (type == self.type_Just_dead) | (type == self.type_Dead), 0.,
                        np.maximum(2. * max_unloading_rate * C_sucrose_root * phloem_exchange_surface / (
                        self.Km_unloading + C_sucrose_root), 0))
    

    @rate
    def _diffusion_AA_phloem(self, label, amino_acids_consumption_by_growth, deficit_AA, AA, phloem_AA, phloem_exchange_surface, soil_temperature, living_struct_mass, symplasmic_volume, phloem_volume):
        """ Passive radial diffusion between phloem and cortex through plasmodesmata """

        reference_rate_of_AA_consumption_by_growth = self.reference_rate_of_AA_consumption_by_growth
        # reference_rate_of_AA_consumption_by_growth = np.where(label==self.label_Apex, reference_rate_of_AA_consumption_by_growth/1, reference_rate_of_AA_consumption_by_growth)

        # permeability_phloem_AA = self.permeability_phloem_AA * (1 + amino_acids_consumption_by_growth / (living_struct_mass * self.massic_reference_rate_of_AA_consumption_by_growth))
        permeability_phloem_AA = self.permeability_phloem_AA * (1 + np.where(type == self.type_Base_of_the_root_system, (deficit_AA) / (reference_rate_of_AA_consumption_by_growth),
                                                                             (amino_acids_consumption_by_growth + deficit_AA) / (reference_rate_of_AA_consumption_by_growth)))
        # permeability_phloem_AA = self.permeability_phloem_AA 

        permeability_phloem_AA *= self.temperature_modification(soil_temperature=soil_temperature,
                                                                    T_ref=self.passive_processes_T_ref,
                                                                    A=self.passive_processes_A,
                                                                    B=self.passive_processes_B,
                                                                    C=self.passive_processes_C)

        flux = permeability_phloem_AA * (np.maximum(0, (phloem_AA * living_struct_mass) / phloem_volume) - np.maximum(0, (AA * living_struct_mass) / symplasmic_volume)) * phloem_exchange_surface

        # return np.where(flux > 0., flux, 0.)
        return flux


    @rate
    def _unloading_AA_phloem(self, phloem_AA, amino_acids_consumption_by_growth, phloem_exchange_surface, soil_temperature, living_struct_mass, phloem_volume):
        Cv_AA_phloem = (phloem_AA * living_struct_mass) / phloem_volume
        
        vmax_unloading_AA_phloem = self.vmax_unloading_AA_phloem * (1 + amino_acids_consumption_by_growth / self.reference_rate_of_AA_consumption_by_growth)
        vmax_unloading_AA_phloem *= self.temperature_modification(soil_temperature=soil_temperature,
                                                            T_ref=self.active_processes_T_ref,
                                                            A=self.active_processes_A,
                                                            B=self.active_processes_B,
                                                            C=self.active_processes_C)
        
        return np.where(vmax_unloading_AA_phloem > 0., np.minimum(vmax_unloading_AA_phloem * Cv_AA_phloem * phloem_exchange_surface / (
                    self.km_unloading_AA_phloem + Cv_AA_phloem), phloem_AA * living_struct_mass / 2), 
                    0.)


    # @note CONCENTRATIONS BALANCE

    @state
    def _C_hexose_root(self, C_hexose_root, living_struct_mass, hexose_exudation, hexose_uptake_from_soil,
                           mucilage_secretion, cells_release, maintenance_respiration,
                           hexose_consumption_by_growth, hexose_consumption_by_fungus, hexose_diffusion_from_phloem,
                           hexose_active_production_from_phloem, sucrose_loading_in_phloem,
                           hexose_mobilization_from_reserve, hexose_immobilization_as_reserve, deficit_hexose_root, 
                           AA_synthesis, AA_catabolism, N_metabolic_respiration) -> tuple[float, str, float]:
        """
        Added the following flows to the balance :
        - Amino acid synthesis hexose consumption
        - Amino acid catabolism releasing hexose
        - Nitrogen metabolism related respiration costs
        """

        f = 1e13 # arbitrary
        _hexose_exudation = hexose_exudation * f
        _hexose_uptake_from_soil = hexose_uptake_from_soil * f
        _mucilage_secretion = mucilage_secretion * f
        _cells_release = cells_release * f
        _maintenance_respiration = maintenance_respiration * f
        _hexose_consumption_by_growth = hexose_consumption_by_growth * f
        _hexose_consumption_by_fungus = hexose_consumption_by_fungus * f
        _hexose_diffusion_from_phloem = hexose_diffusion_from_phloem * f
        _hexose_active_production_from_phloem = hexose_active_production_from_phloem * f
        _sucrose_loading_in_phloem = sucrose_loading_in_phloem * f
        _hexose_mobilization_from_reserve = hexose_mobilization_from_reserve * f
        _hexose_immobilization_as_reserve = hexose_immobilization_as_reserve * f
        _deficit_hexose_root = deficit_hexose_root * f
        _AA_synthesis = AA_synthesis * f
        _AA_catabolism = AA_catabolism * f
        _N_metabolic_respiration = N_metabolic_respiration * f

        inflow =  (_hexose_uptake_from_soil
                + _hexose_diffusion_from_phloem
                + _hexose_active_production_from_phloem
                + _hexose_mobilization_from_reserve
                + _AA_catabolism * self.r_hexose_AA)
        
        outflow = (_hexose_exudation
                + _mucilage_secretion
                + _cells_release
                + _maintenance_respiration / 6.
                + _hexose_consumption_by_growth
                + _hexose_consumption_by_fungus
                + 2. * _sucrose_loading_in_phloem
                + _hexose_immobilization_as_reserve
                + _deficit_hexose_root
                + _AA_synthesis * self.r_hexose_AA
                + _N_metabolic_respiration / 6.)
        
        netflow = inflow - outflow

        _living_struct_mass = 1e6 * living_struct_mass # µg

        derivative = (self.time_step / _living_struct_mass) * netflow
        derivative = derivative * 1e-7
        raw_balance = C_hexose_root + derivative
        
        is_neg = raw_balance < 0.0
        deficit = np.where(is_neg, -raw_balance * (living_struct_mass / self.time_step), 0.0)

        balance = np.where(is_neg, 0.0, raw_balance)

        assert not np.any(np.isnan(balance))
        assert not np.any(np.isinf(balance))
        assert not np.any(np.isnan(deficit))
        assert not np.any(np.isinf(deficit))
        assert np.all(np.abs(raw_balance - (balance - deficit * (self.time_step / living_struct_mass))) < 1e-15)

        return balance, 'deficit_hexose_root', deficit


    @state
    def _AA(self, vertex_index, AA, living_struct_mass, diffusion_AA_phloem, unloading_AA_phloem, loading_AA_phloem, import_AA, diffusion_AA_soil, diffusion_AA_xylem, export_AA, AA_synthesis,
                  amino_acids_consumption_by_growth, storage_synthesis, storage_catabolism, AA_catabolism, deficit_AA) -> tuple[float, str, float]:
        
        f = 1e13 # arbitrary
        _diffusion_AA_phloem = diffusion_AA_phloem * f
        _unloading_AA_phloem = unloading_AA_phloem * f
        _loading_AA_phloem = loading_AA_phloem * f
        _import_AA = import_AA * f
        _AA_synthesis = AA_synthesis * f
        _storage_catabolism = storage_catabolism * f
        _diffusion_AA_soil = diffusion_AA_soil * f
        _diffusion_AA_xylem = diffusion_AA_xylem * f
        _export_AA = export_AA * f
        _amino_acids_consumption_by_growth = amino_acids_consumption_by_growth * f
        _storage_synthesis = storage_synthesis * f
        _AA_catabolism = AA_catabolism * f
        _deficit_AA = deficit_AA * f

        inflow = (_diffusion_AA_phloem
                + _unloading_AA_phloem
                + _diffusion_AA_xylem
                + _import_AA
                + _AA_synthesis
                + _storage_catabolism * self.r_AA_stor)
        
        outflow = (_diffusion_AA_soil
                + _loading_AA_phloem
                + _export_AA
                + _amino_acids_consumption_by_growth
                + _storage_synthesis * self.r_AA_stor
                + _AA_catabolism
                + _deficit_AA)
        
        netflow = inflow - outflow

        _living_struct_mass = 1e6 * living_struct_mass # µg

        derivative = (self.time_step / _living_struct_mass) * netflow
        derivative = derivative * 1e-7
        raw_balance = AA + derivative

        is_neg = raw_balance < 0.0
        deficit = np.where(is_neg, -raw_balance * (living_struct_mass / self.time_step), 0.0)
        # deficit = np.where(deficit > 1e-20, deficit, 0.0)

        balance = np.where(is_neg, 0.0, raw_balance)

        # if np.any(is_neg):
        #     print(vertex_index, AA, diffusion_AA_phloem, unloading_AA_phloem, import_AA, diffusion_AA_soil, export_AA, AA_synthesis,
        #           amino_acids_consumption_by_growth, storage_synthesis, storage_catabolism, AA_catabolism, deficit_AA)

        return balance, 'deficit_AA', deficit

    
    @state
    def _C_solutes_phloem(self, C_sucrose_root, phloem_AA):
        """
        Sucrose could not be included before because phloem massic concentrations were not computed with volumic considerations
        """
        ions_proportion = 0.4 # To account for high 300 mM concentrations of potassium in phloem sap, related to sucrose symport co-transport Diant et al. 2010
        return (C_sucrose_root) / (1 - ions_proportion) + phloem_AA
    

    @totalstate
    def _Cv_sucrose_average(self, C_sucrose_root, living_struct_mass, phloem_volume):
        total_amount = 0
        total_volume = 0
        for vid in living_struct_mass.keys():
            if living_struct_mass[vid] > 0:
                total_amount += C_sucrose_root[vid] * living_struct_mass[vid]
                total_volume += phloem_volume[vid]
        return total_amount / total_volume

    @totalstate
    def _Cv_hexose_average(self, C_hexose_root, living_struct_mass, symplasmic_volume):
        total_amount = 0
        total_volume = 0
        for vid in living_struct_mass.keys():
            if living_struct_mass[vid] > 0:
                total_amount += C_hexose_root[vid] * living_struct_mass[vid]
                total_volume += symplasmic_volume[vid]
        return total_amount / total_volume

    # @note Disable processes from base components

    @state
    def _C_sucrose_root(self):
        """
        Handled by the heterogeneous axial transport model now
        """
        return
    
    @stepinit
    def shoot_sucrose_supply_and_spreading(self):
        """
        Handled by the heterogeneous axial transport model now
        """
        return

    @state
    def _Cv_sucrose_root(self, C_sucrose_root, living_struct_mass, phloem_volume):
        return np.where(phloem_volume > 0., C_sucrose_root * living_struct_mass / np.where(phloem_volume > 0., phloem_volume, 1.), 0.)
    
    
    @state
    def _Cv_hexose_root(self, C_hexose_root, living_struct_mass, symplasmic_volume):
        return np.where(symplasmic_volume > 0., C_hexose_root * living_struct_mass / np.where(symplasmic_volume > 0., symplasmic_volume, 1.), 0.)
    

    # @note balance check functions

    @totalstate
    def check_balance(self):
        """
        Three-section carbon balance diagnostic.

        Section 1 – Explicit symplasm pools (@state-updated):
            C_hexose_root, C_hexose_reserve, AA, storage_protein
        Section 2 – Implicit axial vessel pools (implicit-solver-updated):
            C_sucrose_root, phloem_AA, xylem_AA
        Section 3 – Total system (Sections 1+2 combined)

        For each section, residual = actual_Δpool - net_boundary_C_flows × dt.
        A non-zero residual in a section points to that resolution path.
        """
        props = self.props
        dt    = self.time_step
        self.cumulated_time += dt

        if not hasattr(self, "track_residuals"):
            self.track_residuals = []

        vertex_index = props["vertex_index"] 
        focus_vids  = np.asarray(props["focus_elements"], dtype=np.int64)        # (n,)
        focus_glob_idx  = vertex_index.indices_of(focus_vids)

        # ── Structural mass ──────────────────────────────────────────────────
        lsm  = props["living_struct_mass"].values_array()[focus_glob_idx]

        # ── Pool snapshots (mol-C per pool) ──────────────────────────────────
        # mol-C = C_weight × massic_conc × lsm
        # C-weights: hexose=7, sucrose=12, reserve=6, storage_protein≈5×65, AA=5
        _hex_arr  = props["C_hexose_root"].values_array()[focus_glob_idx]
        _suc_arr  = props["C_sucrose_root"].values_array()[focus_glob_idx]
        _res_arr  = props["C_hexose_reserve"].values_array()[focus_glob_idx]
        _stor_arr = props["storage_protein"].values_array()[focus_glob_idx]
        _AA_arr   = props["AA"].values_array()[focus_glob_idx]
        _phAA_arr = props["phloem_AA"].values_array()[focus_glob_idx]
        _xyAA_arr = props["xylem_AA"].values_array()[focus_glob_idx]

        c_hex  = (6    * _hex_arr  * lsm).sum()
        c_suc  = (12   * _suc_arr  * lsm).sum()
        c_res  = (6    * _res_arr  * lsm).sum()
        c_stor = (5*65 * _stor_arr * lsm).sum()
        c_AA   = (5    * _AA_arr   * lsm).sum()
        c_phAA = (5    * _phAA_arr * lsm).sum()
        c_xyAA = (5    * _xyAA_arr * lsm).sum()

        _cur = dict(hex=c_hex, suc=c_suc, res=c_res, stor=c_stor,
                    AA=c_AA, phAA=c_phAA, xyAA=c_xyAA)

        if not hasattr(self, '_prev_pool_C'):
            self._prev_pool_C = _cur
        prev = self._prev_pool_C

        # ── Section 2: explicit symplasm (mol-C/s, positive = inflow to root) ─
        #   hex         : via _C_hexose_root @state
        #   reserve     : via _C_hexose_reserve @state
        #   AA symplasm : via _AA @state
        #   storage_prot: via _storage_protein @state
        #
        #   Boundary terms that cross the root system boundary:
        stp_symp_maint_resp  = -1  * props["maintenance_respiration"].values_array()[focus_glob_idx].sum()
        stp_symp_Nresp       = -1  * props["N_metabolic_respiration"].values_array()[focus_glob_idx].sum()
        stp_symp_hex_growth  = -6  * props["hexose_consumption_by_growth"].values_array()[focus_glob_idx].sum()
        stp_symp_hex_fungus  = -6  * props["hexose_consumption_by_fungus"].values_array()[focus_glob_idx].sum()
        stp_symp_hex_exud    = -6  * props["hexose_exudation"].values_array()[focus_glob_idx].sum()
        stp_symp_ph_hex_exud = -6  * props["phloem_hexose_exudation"].values_array()[focus_glob_idx].sum()
        stp_symp_hex_uptake  = +6  * props["hexose_uptake_from_soil"].values_array()[focus_glob_idx].sum()
        stp_symp_ph_hex_uptk = +6  * props["phloem_hexose_uptake_from_soil"].values_array()[focus_glob_idx].sum()
        stp_symp_mucilage    = -6  * props["mucilage_secretion"].values_array()[focus_glob_idx].sum()
        stp_symp_cells       = -6  * props["cells_release"].values_array()[focus_glob_idx].sum()
        stp_symp_import_AA   = +5  * props["import_AA"].values_array()[focus_glob_idx].sum()
        stp_symp_diff_AA_sl  = -5  * props["diffusion_AA_soil"].values_array()[focus_glob_idx].sum()
        stp_symp_AA_growth   = -5  * props["amino_acids_consumption_by_growth"].values_array()[focus_glob_idx].sum()
        stp_symp_aplastic_AA = -5  * props["apoplastic_AA_soil_xylem"].values_array()[focus_glob_idx].sum()

        # Internal transport flows between symplasm and xylem / phloem pools
        # Sucrose
        stp_symp_diff_suc_tosymp  = 6 * props["hexose_diffusion_from_phloem"].values_array()[focus_glob_idx].sum()
        stp_symp_active_suc_tosymp  = 6 * props["hexose_active_production_from_phloem"].values_array()[focus_glob_idx].sum()
        stp_symp_active_suc_toph  = - 12 * props["sucrose_loading_in_phloem"].values_array()[focus_glob_idx].sum()

        # Phloem AA
        stp_symp_active_phaa_toph = - 5 * props["loading_AA_phloem"].values_array()[focus_glob_idx].sum()
        stp_symp_diff_phaa_tosymp  = 5 * props["diffusion_AA_phloem"].values_array()[focus_glob_idx].sum()
        stp_symp_active_phaa_tosymp  = 5 * props["unloading_AA_phloem"].values_array()[focus_glob_idx].sum()

        # Xylem AA
        stp_symp_active_xyaa_toxy = - 5 * props["export_AA"].values_array()[focus_glob_idx].sum()
        stp_symp_diff_xyaa_tosymp = 5 * props["diffusion_AA_xylem"].values_array()[focus_glob_idx].sum()

        # internal transfers from axial vessels into symplasm are NOT boundaries
        # (they cancel between sections), so we leave them out here.
        stp_symp_root_soil_boundary_rate = (stp_symp_maint_resp + stp_symp_Nresp +
                            stp_symp_hex_growth + stp_symp_hex_fungus +
                            stp_symp_hex_exud + stp_symp_ph_hex_exud +
                            stp_symp_hex_uptake + stp_symp_ph_hex_uptk +
                            stp_symp_mucilage + stp_symp_cells +
                            stp_symp_import_AA + stp_symp_diff_AA_sl +
                            stp_symp_AA_growth + stp_symp_aplastic_AA)
        
        stp_symp_symplasm_boundary_rate = (stp_symp_maint_resp + stp_symp_Nresp +
                            stp_symp_hex_growth + stp_symp_hex_fungus +
                            stp_symp_hex_exud + stp_symp_hex_uptake +
                            stp_symp_mucilage + stp_symp_cells +
                            stp_symp_import_AA + stp_symp_diff_AA_sl +
                            stp_symp_AA_growth + 
                            stp_symp_diff_suc_tosymp + stp_symp_active_suc_tosymp + stp_symp_active_suc_toph +
                            stp_symp_active_phaa_toph + stp_symp_diff_phaa_tosymp + stp_symp_active_phaa_tosymp + 
                            stp_symp_active_xyaa_toxy + stp_symp_diff_xyaa_tosymp
                            )

        clipped_deficit_rate = (6  * props["deficit_hexose_root"].values_array()[focus_glob_idx]
                              + 6  * props["deficit_hexose_reserve"].values_array()[focus_glob_idx]
                              + 12 * props["deficit_sucrose_root"].values_array()[focus_glob_idx]
                              + 5  * props["deficit_AA"].values_array()[focus_glob_idx]
                              + 5  * props["deficit_AA_phloem"].values_array()[focus_glob_idx]
                              + 5  * props["deficit_AA_xylem"].values_array()[focus_glob_idx]
                              ).sum()
        clipped_deficit_rate_symp = (6  * props["deficit_hexose_root"].values_array()[focus_glob_idx]
                              + 6  * props["deficit_hexose_reserve"].values_array()[focus_glob_idx]
                              + 5  * props["deficit_AA"].values_array()[focus_glob_idx]
                              ).sum()
        current_deficit_amount = dt * clipped_deficit_rate
        current_deficit_amount_symp = dt * clipped_deficit_rate_symp
        if not hasattr(self, 'previous_deficit_amount'):
            self.previous_deficit_amount = 0.0
        if not hasattr(self, 'previous_deficit_amount_symp'):
            self.previous_deficit_amount_symp = 0.0

        incremental_deficit = current_deficit_amount - self.previous_deficit_amount
        incremental_deficit_symp = current_deficit_amount_symp - self.previous_deficit_amount_symp

        stp_symp_pool_delta = (c_hex - prev['hex']) + (c_res - prev['res']) + \
                        (c_AA  - prev['AA'])  + (c_stor - prev['stor'])
        stp_symp_residual   = stp_symp_pool_delta - dt * stp_symp_symplasm_boundary_rate - incremental_deficit_symp
        stp_symp_C_amount = c_hex + c_AA + c_res + c_stor
        stp_symp_residual_percentage = 100 * stp_symp_residual / stp_symp_C_amount

        # ── Section 3: implicit axial vessels ──────────────────────────────
        #   sucrose phloem (C_sucrose_root), phloem AA, xylem AA
        stp_vessels_suc_shoot_inflow  = -12 * props["sucrose_root_to_shoot_phloem"][1]
        stp_vessels_phAA_shoot_inflow = -5  * props["AA_root_to_shoot_phloem"][1]
        stp_vessels_xyAA_shoot_outflow = 5  * props["AA_root_to_shoot_xylem"][1]

        radial_suc_inflows = 12 * (- props["hexose_diffusion_from_phloem"].values_array()[focus_glob_idx] / 2.
                            - props["hexose_active_production_from_phloem"].values_array()[focus_glob_idx] / 2.
                            - props["phloem_hexose_exudation"].values_array()[focus_glob_idx] / 2.
                            + props["sucrose_loading_in_phloem"].values_array()[focus_glob_idx]
                            + props["phloem_hexose_uptake_from_soil"].values_array()[focus_glob_idx] / 2.).sum()
        radial_xyAA_inflows = 5 * (props["export_AA"].values_array()[focus_glob_idx] - props["apoplastic_AA_soil_xylem"].values_array()[focus_glob_idx] - props["diffusion_AA_xylem"].values_array()[focus_glob_idx]).sum()
        radial_phAA_inflows = 5 * (props["loading_AA_phloem"].values_array()[focus_glob_idx] - props["diffusion_AA_phloem"].values_array()[focus_glob_idx] - props["unloading_AA_phloem"].values_array()[focus_glob_idx]).sum()

        stp_vessels_shoot_root_boundary_rate = stp_vessels_suc_shoot_inflow + stp_vessels_phAA_shoot_inflow - stp_vessels_xyAA_shoot_outflow
        stp_vessels_radial_vessels_boundary_rate = radial_suc_inflows + radial_phAA_inflows + radial_xyAA_inflows

        stp_vessels_pool_delta = ((c_suc  - prev['suc']) +
                         (c_phAA - prev['phAA']) +
                         (c_xyAA - prev['xyAA']))

        stp_vessels_residual   = stp_vessels_pool_delta - dt * (stp_vessels_shoot_root_boundary_rate + stp_vessels_radial_vessels_boundary_rate)
        stp_vessels_C_amount = c_suc + c_phAA + c_xyAA
        stp_vessels_residual_percentage = 100 * stp_vessels_residual / stp_vessels_C_amount 
        stp_vessels_suc_residual = (c_suc - prev['suc']) - dt * (stp_vessels_suc_shoot_inflow + radial_suc_inflows)
        stp_vessels_phAA_residual = (c_phAA - prev['phAA']) - dt * (stp_vessels_phAA_shoot_inflow + radial_phAA_inflows)
        stp_vessels_xyAA_residual = (c_xyAA - prev['xyAA']) - dt * (-stp_vessels_xyAA_shoot_outflow + radial_xyAA_inflows)

        # ── Section 4: total system ──────────────────────────────────────────
        # Deficit accounting: cumulated deficit represents C that should have
        # left the system but was clipped to zero. The incremental deficit is
        # the change in this obligation between steps.

        total_boundary_rate = stp_symp_root_soil_boundary_rate + stp_vessels_shoot_root_boundary_rate
        total_pool_delta = sum(_cur.values()) - sum(prev.values())
        total_residual = total_pool_delta - dt * total_boundary_rate - incremental_deficit
        total_C_amount = sum(_cur.values())
        total_residual_percentage = 100 * total_residual / total_C_amount

        self.track_residuals.append({"t": self.cumulated_time, "r_total": total_residual, "t_total": total_C_amount, "r_stp_symp":stp_symp_residual, "r_stp_vessels": stp_vessels_residual,
            "r_stp_vessels_suc": stp_vessels_suc_residual, "r_stp_vessels_phAA": stp_vessels_phAA_residual, "r_stp_vessels_xyAA": stp_vessels_xyAA_residual})

        if total_residual_percentage > 5.:
            print("     SIGNIFICATIVE ERROR IN THE C BALANCE:")
            print(f"    Total pool delta           : {total_pool_delta:+.7e}")
            print(f"    Total boundary x dt        : {dt*total_boundary_rate:+.7e}")
            print(f"    Deficit delta (cur-prev)   : {incremental_deficit:+.7e}  (cur={current_deficit_amount:.3e}, prev={self.previous_deficit_amount:.3e})")
            print(f"    TOTAL RESIDUAL             : {total_residual:+.7e}, p = {total_residual_percentage:+.6e}% (=S1+S2 - internal transfers cancel)")
            print(f"    TOTAL SYMPLASM RESIDUAL    : {stp_symp_residual:+.6e}, p = {stp_symp_residual_percentage:+.6e}%")
            print(f"    TOTAL AXIAL TRANSPORT RESIDUAL : {stp_vessels_residual:+.6e}, p = {stp_vessels_residual_percentage:+.6e}%")
            print(f"    SUC AXIAL TRANSPORT RESIDUAL : {stp_vessels_suc_residual:+.6e}, {100 * stp_vessels_suc_residual / c_suc:+.6e}%")
            print(f"    phAA AXIAL TRANSPORT RESIDUAL : {stp_vessels_phAA_residual:+.6e}, {100 * stp_vessels_phAA_residual / c_phAA:+.6e}%")
            print(f"    xyAA AXIAL TRANSPORT RESIDUAL : {stp_vessels_xyAA_residual:+.6e}, {100 * stp_vessels_xyAA_residual / c_xyAA:+.6e}%")
            print("===========================")
       
        df = pd.DataFrame(self.track_residuals)
        rmse = np.sqrt(np.sum(df["r_total"].to_numpy() ** 2))
        current_percentage = 100 * rmse / total_C_amount
        if current_percentage > 5.:
            print(f"WARNING, RMSE {rmse:+.3e} mol C, in percentage total C balance: {current_percentage:+.3e}%")

        assert current_percentage < 15., f"ERROR, RMSE exceeds 15% of current C poool: RMSE {rmse:+.3e} accounts for {current_percentage}%"
        self._prev_pool_C             = _cur
        self.previous_deficit_amount  = current_deficit_amount
        self.previous_deficit_amount_symp = current_deficit_amount_symp


    def compute_root_system_C_content(self):
        segment_C_content = (6*self.props["C_hexose_root"].values_array()
                             + 12*self.props["C_sucrose_root"].values_array()
                             + 6*self.props["C_hexose_reserve"].values_array()
                             + 5*65*self.props["storage_protein"].values_array()
                             + 5*self.props["AA"].values_array()
                             + 5*self.props["phloem_AA"].values_array()
                             + 5*self.props["xylem_AA"].values_array()) * self.props["living_struct_mass"].values_array()

        assert (not np.any(np.isnan(segment_C_content))) and (not np.any(np.isinf(segment_C_content))) and (not np.any(segment_C_content < 0.)), "Some segments have nan, infinite or negative C balance"

        return segment_C_content
