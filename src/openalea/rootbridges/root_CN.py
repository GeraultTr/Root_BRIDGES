from dataclasses import dataclass
import numpy as np
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

    sucrose_root_to_shoot_phloem: float =       declare(default=0, unit="mol.time_step-1", unit_comment="of sucrose", description="",
                                                min_value="", max_value="", value_comment="", references="", DOI="",
                                                variable_type="plant_scale_state", by="model_nitrogen", state_variable_type="", edit_by="user")
    Cv_sucrose_average: float =                   declare(default=1., unit="mol.m-3", unit_comment="of amino acids", description="", 
                                                min_value="", max_value="", value_comment="", references="", DOI="",
                                                variable_type="plant_scale_state", by="model_nitrogen", state_variable_type="", edit_by="user")
    Cv_hexose_average: float =                   declare(default=1., unit="mol.m-3", unit_comment="of amino acids", description="", 
                                                min_value="", max_value="", value_comment="", references="", DOI="",
                                                variable_type="plant_scale_state", by="model_nitrogen", state_variable_type="", edit_by="user")
    
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

        self.previous_C_amount_in_the_root_system = self.compute_root_system_C_content()
        # self.total_root_sucrose_and_living_struct_mass() # Needed otherwise first shoot unloading will be unrealistic

        self.solute_configs["C_sucrose_root"] = {
        "solute_massic_concentration_prop": "C_sucrose_root",
        "solute_massic_concentration_symplasm": "hexose_diffusion_from_phloem",
        "diffusive_flux_name": "hexose_diffusion_from_phloem",
        "diffusive_flux_conversion": - 1 /2,
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
    def _hexose_diffusion_from_phloem(self, type, length, phloem_exchange_surface, C_sucrose_root, C_hexose_root,
                                             hexose_consumption_by_growth, living_struct_mass, symplasmic_volume, phloem_volume, soil_temperature):
        """
        Superimposing original, staying with a massic concentration gradient as fist approximation to avoid changing parameters
        """

        Cv_sucrose_root = C_sucrose_root * living_struct_mass / phloem_volume
        Cv_hexose_root = C_hexose_root * living_struct_mass / symplasmic_volume

        # if Cv_sucrose_root <= Cv_hexose_root / 2:
        #     print("sucrose limits", Cv_sucrose_root, Cv_hexose_root / 2., living_struct_mass, hexose_consumption_by_growth )
        #     return 0
        # else:
        # print(phloem_volume, symplasmic_volume)
        # print("gradient", Cv_sucrose_root,  Cv_hexose_root / 2., living_struct_mass, hexose_consumption_by_growth)
        phloem_permeability = self.diffusion_phloem * (1 + hexose_consumption_by_growth /
                                                            self.reference_rate_of_hexose_consumption_by_growth)

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
        max_unloading_rate = 0 * self.max_unloading_rate * (1 + hexose_consumption_by_growth /
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
    def _diffusion_AA_phloem(self, amino_acids_consumption_by_growth, AA, phloem_AA, phloem_exchange_surface, soil_temperature, living_struct_mass, symplasmic_volume, phloem_volume):
        """ Passive radial diffusion between phloem and cortex through plasmodesmata """

        permeability_phloem_AA = self.permeability_phloem_AA * (1 + amino_acids_consumption_by_growth / (self.reference_rate_of_AA_consumption_by_growth))

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
        # print({k: v for k, v in locals().items() if k != 'self'})

        balance = C_hexose_root + (self.time_step / living_struct_mass) * (
                - hexose_exudation
                + hexose_uptake_from_soil
                - mucilage_secretion
                - cells_release
                - maintenance_respiration / 6.
                - hexose_consumption_by_growth
                - hexose_consumption_by_fungus
                + hexose_diffusion_from_phloem
                + hexose_active_production_from_phloem
                - 2. * sucrose_loading_in_phloem
                + hexose_mobilization_from_reserve
                - hexose_immobilization_as_reserve
                - deficit_hexose_root
                - AA_synthesis * self.r_hexose_AA
                + AA_catabolism * self.r_hexose_AA
                - N_metabolic_respiration / 6.)
        
        deficit = - balance * living_struct_mass / self.time_step
        deficit = np.where(deficit > 1e-20, deficit, 0.)
        balance = np.maximum(balance, 0.)

        return balance, 'deficit_hexose_root', deficit

        
    @state
    def _AA(self, AA, living_struct_mass, diffusion_AA_phloem, unloading_AA_phloem, import_AA, diffusion_AA_soil, export_AA, AA_synthesis,
                  amino_acids_consumption_by_growth, storage_synthesis, storage_catabolism, AA_catabolism, deficit_AA) -> tuple[float, str, float]:
        
        f = 1e13 # arbitrary
        diffusion_AA_phloem = diffusion_AA_phloem * f
        unloading_AA_phloem = unloading_AA_phloem * f
        import_AA = import_AA * f
        AA_synthesis = AA_synthesis * f
        storage_catabolism = storage_catabolism * f
        diffusion_AA_soil = diffusion_AA_soil * f
        export_AA = export_AA * f
        amino_acids_consumption_by_growth = amino_acids_consumption_by_growth * f
        storage_synthesis = storage_synthesis * f
        AA_catabolism = AA_catabolism * f
        deficit_AA = deficit_AA * f

        inflow = (diffusion_AA_phloem
                + unloading_AA_phloem
                + import_AA
                + AA_synthesis
                + storage_catabolism * self.r_AA_stor)
        
        outflow = (diffusion_AA_soil
                + export_AA
                + amino_acids_consumption_by_growth
                + storage_synthesis * self.r_AA_stor
                + AA_catabolism
                + deficit_AA)
        
        netflow = inflow - outflow

        living_struct_mass = 1e6 * living_struct_mass # µg

        derivative = (self.time_step / living_struct_mass) * netflow
        derivative = derivative * 1e-7
        raw_balance = AA + derivative

        is_neg = raw_balance < 0.0
        deficit = np.where(is_neg, -raw_balance * (living_struct_mass / self.time_step), 0.0)
        # deficit = np.where(deficit > 1e-20, deficit, 0.0)

        balance = np.where(is_neg, 0.0, raw_balance)

        return balance, 'deficit_AA', deficit

    # @state
    # def _AA(self, AA, living_struct_mass, diffusion_AA_phloem, unloading_AA_phloem, import_AA, diffusion_AA_soil, export_AA, AA_synthesis,
    #               amino_acids_consumption_by_growth, storage_synthesis, storage_catabolism, AA_catabolism, deficit_AA) -> tuple[float, str, float]:
        
    #     balance =  AA + (self.time_step / living_struct_mass) * (
    #             diffusion_AA_phloem
    #             + unloading_AA_phloem
    #             + import_AA
    #             - diffusion_AA_soil
    #             - export_AA
    #             + AA_synthesis
    #             - amino_acids_consumption_by_growth
    #             - storage_synthesis * self.r_AA_stor
    #             + storage_catabolism * self.r_AA_stor
    #             - AA_catabolism
    #             - deficit_AA)
        
    #     deficit = - balance * living_struct_mass / self.time_step
    #     deficit = np.where(deficit > 1e-20, deficit, 0.)
    #     balance = np.maximum(balance, 0.)

    #     return balance, 'deficit_AA', deficit

    
    @state
    def _C_solutes_phloem(self, C_sucrose_root, phloem_AA):
        """
        Sucrose could not be included before because phloem massic concentrations were not computed with volumic considerations
        """
        ions_proportion = 0.4 # To account for high 300 mM concentrations of potassium in phloem sap, related to sucrose symport co-transport Diant et al. 2010
        return (C_sucrose_root + phloem_AA) / (1 - ions_proportion)
    

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
    

        

    