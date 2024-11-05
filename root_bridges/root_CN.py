from dataclasses import dataclass
from metafspm.component_factory import *
from metafspm.component import declare

from rhizodep.root_carbon import RootCarbonModel
from root_cynaps.root_nitrogen import RootNitrogenModel


family = "metabolic"


@dataclass
class RootCNUnified(RootCarbonModel, RootNitrogenModel):

    # INPUTS
    # FROM GROWTH MODEL
    amino_acids_consumption_by_growth: float = declare(default=0., unit="mol.s-1", unit_comment="", description="amino_acids consumption rate by growth processes", 
                                                 min_value="", max_value="", value_comment="", references="", DOI="",
                                                  variable_type="input", by="model_growth", state_variable_type="", edit_by="user")

    # STATE VARIABLES
    N_metabolic_respiration: float = declare(default=0., unit="mol.s-1", unit_comment="of carbon", description="Respiration related to nitrogen metabolism", 
                                            min_value="", max_value="", value_comment="", references="", DOI="",
                                             variable_type="state_variable", by="model_carbon", state_variable_type="extensive", edit_by="user")
    total_hexose_diffusion_from_phloem: float = declare(default=0., unit="umol of C.g-1 mstruc.h-1", unit_comment="", description="Property computed to compare with shoot model unloading",
                                    min_value="", max_value="", value_comment="", references="", DOI="",
                                    variable_type="plant_scale_state", by="model_carbon", state_variable_type="", edit_by="user")
    deficit_AA_root: float = declare(default=0., unit="mol.s-1", unit_comment="of amino acids", description="Amino acids deficit rate in root", 
                                        min_value="", max_value="", value_comment="", references="Hypothesis of no initial deficit", DOI="",
                                         variable_type="state_variable", by="model_nitrogen", state_variable_type="extensive", edit_by="user")
    
    # PARAMETERS
    r_hexose_AA: float = declare(default=4.5/6, unit="adim", unit_comment="mol of hexose per mol of amino acids in roots", description="stoechiometric ratio during amino acids synthesis for hexose consumption", 
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
        self.g = g
        self.props = self.g.properties()
        self.time_step = time_step
        self.choregrapher.add_time_and_data(instance=self, sub_time_step=self.time_step, data=self.props)
        self.vertices = self.g.vertices(scale=self.g.max_scale())

        # Before any other operation, we apply the provided scenario by changing default parameters and initialization
        self.apply_scenario(**scenario)
        self.link_self_to_mtg()
        self.initiate_heterogeneous_variables()

        self.previous_C_amount_in_the_root_system = self.compute_root_system_C_content()

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


    @potential
    @state
    def _C_hexose_root(self, C_hexose_root, struct_mass, living_root_hairs_struct_mass, hexose_exudation, hexose_uptake_from_soil,
                           mucilage_secretion, cells_release, maintenance_respiration,
                           hexose_consumption_by_growth, hexose_diffusion_from_phloem,
                           hexose_active_production_from_phloem, sucrose_loading_in_phloem,
                           hexose_mobilization_from_reserve, hexose_immobilization_as_reserve, deficit_hexose_root, 
                           AA_synthesis, AA_catabolism, N_metabolic_respiration):
        """
        Added the following flows to the balance :
        - Amino acid synthesis hexose consumption
        - Amino acid catabolism releasing hexose
        - Nitrogen metabolism related respiration costs
        """
        return C_hexose_root + (self.time_step / (struct_mass + living_root_hairs_struct_mass)) * (
                - hexose_exudation
                + hexose_uptake_from_soil
                - mucilage_secretion
                - cells_release
                - maintenance_respiration / 6.
                - hexose_consumption_by_growth
                + hexose_diffusion_from_phloem
                + hexose_active_production_from_phloem
                - 2. * sucrose_loading_in_phloem
                + hexose_mobilization_from_reserve
                - hexose_immobilization_as_reserve
                - deficit_hexose_root
                - AA_synthesis * self.r_hexose_AA
                + AA_catabolism / self.r_hexose_AA
                - N_metabolic_respiration / 6.)

    @potential
    @state
    def _AA(self, AA, struct_mass, diffusion_AA_phloem, import_AA, diffusion_AA_soil, export_AA, AA_synthesis, 
            storage_synthesis, storage_catabolism, AA_catabolism, struct_mass_produced, amino_acids_consumption_by_growth, deficit_AA_root):
        """
        EDIT : replaced structural nitrogen synthesis by rhizodep input in balance
        """

        if struct_mass > 0:
            return AA + (self.time_step / struct_mass) * (
                    diffusion_AA_phloem
                    + import_AA
                    - diffusion_AA_soil
                    - export_AA
                    + AA_synthesis
                    - storage_synthesis * self.r_AA_stor
                    + storage_catabolism / self.r_AA_stor
                    - AA_catabolism
                    - amino_acids_consumption_by_growth
                    - deficit_AA_root
            ) 
        else:
            return 0

    @deficit
    @state
    def _deficit_AA_root(self, AA, struct_mass, living_root_hairs_struct_mass):
        if AA < 0:
            return - AA * (struct_mass + living_root_hairs_struct_mass) / self.time_step
        else:
            return 0.
        
    @actual
    @state
    def _threshold_C_sucrose_root(self):
        for vid in self.props["focus_elements"]:
            if self.AA[vid] < 0.:
                self.AA[vid] = 0.

    @rate
    def _hexose_active_production_from_phloem(self, length, phloem_exchange_surface, C_sucrose_root, C_hexose_root,
                                              hexose_consumption_by_growth, soil_temperature):
        # We consider all the cases where no net exchange should be allowed:
        if length <= 0. or phloem_exchange_surface <= 0. or type == "Just_dead" or type == "Dead":
            return 0

        else:
            max_unloading_rate = self.max_unloading_rate * (1 + hexose_consumption_by_growth /
                                                            self.reference_rate_of_hexose_consumption_by_growth)
            max_unloading_rate *= self.temperature_modification(soil_temperature=soil_temperature,
                                                                T_ref=self.phloem_unloading_T_ref,
                                                                A=self.phloem_unloading_A,
                                                                B=self.phloem_unloading_B,
                                                                C=self.phloem_unloading_C)
            
            return max(2. * max_unloading_rate * C_sucrose_root * phloem_exchange_surface / (
                    self.Km_unloading + C_sucrose_root), 0)
