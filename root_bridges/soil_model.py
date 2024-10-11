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


    # STATE VARIABLES

    # STATES INITIALIZATION
    # C related
    POC: float = declare(default=2.e-3, unit="adim", unit_comment="gC per g of dry soil", description="Particulate Organic Carbon massic concentration in soil", 
                                        value_comment="", references="Fischer et al. 1966", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="intensive", edit_by="user")
    MAOC: float = declare(default=8.e-3, unit="adim", unit_comment="gC per g of dry soil", description="Mineral Associated Organic Carbon in soil", 
                                        value_comment="", references="Fischer et al. 1966", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="intensive", edit_by="user")
    DOC: float = declare(default=0.5e-3, unit="adim", unit_comment="gC per g of dry soil", description="Dissolved Organic Carbon massic concentration in soil", 
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
    DON: float = declare(default=0.05e-3, unit="adim", unit_comment="gN per g of dry soil", description="Dissolved Organic Nitrogen massic concentration in soil", 
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
    soil_moisture: float = declare(default=0.25, unit="adim", unit_comment="", description="Volumetric proportion of water per volume of soil", 
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
    dry_doil_mass: float = declare(default=1.3e6, unit="g", unit_comment="", description="dry weight of the considered voxel element", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="extensive", edit_by="user")


    # RATES INITIALIZATION
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
    mineral_N_microbial_uptake: float = declare(default=0., unit=".s-1", unit_comment="gN per g of soil per second", description="mineral N uptake by micro organisms", 
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
    
    CN_ratio_POM: float = declare(default=20, unit="adim", unit_comment="gC per gN", description="", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    CN_ratio_MAOM: float = declare(default=10, unit="adim", unit_comment="gC per gN", description="", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    CN_ratio_microbial_biomass: float = declare(default=8, unit="adim", unit_comment="gC per gN", description="", 
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
    
    max_N_uptake_per_microbial_C: float = declare(default=1e-7, unit=".s-1", unit_comment="gN per second per gC of microbial C", description="", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    Km_microbial_N_uptake: float = declare(default=0.01 / 1e3, unit="adim", unit_comment="gN per g of dry soil", description="", 
                                        value_comment="", references="", DOI="",
                                       min_value="", max_value="", variable_type="parameter", by="model_soil", state_variable_type="", edit_by="user")
    
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
        self.voxels["C_mineralN_soil"] = self.voxels["dissolved_mineral_N"] * self.voxels["dry_soil_mass"] / self.voxels["water_volume"]
        self.voxels["C_amino_acids_soil"] = self.voxels["DON"] * self.voxels["dry_soil_mass"] / self.voxels["water_volume"]
    
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
        if microbial_C > self.microbial_C_min and microbial_C < self.microbial_C_max:
            temperature_regulation = self.temperature_modification(soil_temperature=soil_temperature,
                                                                   T_ref=self.microbial_degradation_rate_max_T_ref,
                                                                    A=self.microbial_degradation_rate_max_A,
                                                                    B=self.microbial_degradation_rate_max_B,
                                                                    C=self.microbial_degradation_rate_max_C)
            return (1. + (microbial_C - self.microbial_C_min)  / (self.microbial_C_max - self.microbial_C_min)) * temperature_regulation
        else:
            return 1. * temperature_regulation

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

    @actual
    @rate
    def _mineral_N_microbial_uptake(self, microbial_C, dissolved_mineral_N):
        return self.max_N_uptake_per_microbial_C * microbial_C * dissolved_mineral_N / (self.Km_microbial_N_uptake + dissolved_mineral_N)
    
    # STATES

    @state
    def _POC(self, POC, dry_doil_mass, degradation_POC, cells_release):
        return POC + (self.time_step_in_seconds / dry_doil_mass) * (
            cells_release
            - degradation_POC * dry_doil_mass
        )
    
    @state
    def _PON(self, POC, dry_doil_mass, degradation_POC, cells_release):
        return POC + (self.time_step_in_seconds / dry_doil_mass) * (
            cells_release / self.CN_ratio_root_cells
            - degradation_POC *dry_doil_mass / self.CN_ratio_POM
        )
    
    @state
    def _MAOC(self, MAOC, dry_doil_mass, degradation_MAOC):
        return MAOC + (self.time_step_in_seconds / dry_doil_mass) * (
            - degradation_MAOC * dry_doil_mass
        )
    
    @state
    def _MAON(self, MAON, dry_doil_mass, degradation_MAOC):
        return MAON + (self.time_step_in_seconds / dry_doil_mass) * (
            - degradation_MAOC * dry_doil_mass / self.CN_ratio_MAOM
        )

    @state
    def _DOC(self, DOC, dry_doil_mass, degradation_DOC, hexose_exudation, phloem_hexose_exudation, mucilage_secretion, amino_acids_diffusion_from_roots,  amino_acids_diffusion_from_xylem, amino_acid_uptake):
        return DOC + (self.time_step_in_seconds / dry_doil_mass) * (
            - degradation_DOC * dry_doil_mass
            + hexose_exudation
            + phloem_hexose_exudation
            + mucilage_secretion
            + amino_acids_diffusion_from_roots
            + amino_acids_diffusion_from_xylem
            - amino_acid_uptake
        )
    
    @state
    def _DON(self, DON, DOC, dry_doil_mass, degradation_DOC, amino_acids_diffusion_from_roots,  amino_acids_diffusion_from_xylem, amino_acid_uptake):
        return DON + (self.time_step_in_seconds / dry_doil_mass) * (
            - degradation_DOC * dry_doil_mass * DON / DOC
            + (amino_acids_diffusion_from_roots
            + amino_acids_diffusion_from_xylem - amino_acid_uptake) / self.CN_ratio_amino_acids
        )
    
    @state
    def _microbial_C(self, microbial_C, dry_doil_mass, degradation_POC, degradation_MAOC, degradation_DOC, degradation_microbial_OC):
        """
        For microbial biomass C, the new concentrations results i) from the turnover of this pool, ii) from a fraction of
        the degradation products of each SOC pool, including microbial biomass itself, and iii) from new net inputs, if any:
        """
        return microbial_C + (self.time_step_in_seconds) * (
            degradation_microbial_OC
            - degradation_POC * self.CUE_POC
            - degradation_MAOC * self.CUE_MAOC
            - degradation_DOC * self.CUE_DOC
            - degradation_microbial_OC * self.CUE_MbOC
        )
    
    @state
    def _microbial_N(self, microbial_N, dry_doil_mass, DOC, DON, degradation_POC, degradation_MAOC, degradation_DOC, degradation_microbial_OC):
        """
        For microbial biomass C, the new concentrations results i) from the turnover of this pool, ii) from a fraction of
        the degradation products of each SOC pool, including microbial biomass itself, and iii) from new net inputs, if any:
        """
        return microbial_N + (self.time_step_in_seconds) * (

            degradation_microbial_OC / self.CN_ratio_microbial_biomass
            - degradation_POC * self.CUE_POC / self.CN_ratio_POM
            - degradation_MAOC * self.CUE_MAOC / self.CN_ratio_MAOM
            - degradation_DOC * self.CUE_DOC * DON / DOC
            - degradation_microbial_OC * self.CUE_MbOC/ self.CN_ratio_microbial_biomass
        )
    
    @state
    def _CO2(self, CO2, dry_doil_mass, degradation_POC, degradation_MAOC, degradation_DOC, degradation_microbial_OC):
        return CO2 + (self.time_step_in_seconds) * (
                - degradation_POC * (1 - self.CUE_POC)
                - degradation_MAOC * (1 - self.CUE_MAOC)
                - degradation_microbial_OC * (1 - self.CUE_MbOC)
                - degradation_DOC * (1 - self.CUE_DOC)
            )
    
    @state
    def _dissolved_mineral_N(self, dissolved_mineral_N, dry_doil_mass, degradation_microbial_OC, mineral_N_microbial_uptake, 
                             mineralN_diffusion_from_roots, mineralN_diffusion_from_xylem, mineralN_uptake, fertization_inputs_to_dissolved_mineral_N):
        return dissolved_mineral_N + (self.time_step_in_seconds / dry_doil_mass) * (
            - degradation_microbial_OC / self.CN_ratio_microbial_biomass
            - mineral_N_microbial_uptake 
            + mineralN_diffusion_from_roots
            + mineralN_diffusion_from_xylem
            - mineralN_uptake
            + fertization_inputs_to_dissolved_mineral_N
        )

    @state
    def _C_mineralN_soil(self, dissolved_mineral_N, dry_doil_mass, soil_moisture, voxel_volume):
        return dissolved_mineral_N * (dry_doil_mass * (soil_moisture * voxel_volume)) / 14

    @state
    def _C_amino_acids_soil(self, DOC, dry_doil_mass, soil_moisture, voxel_volume):
        return DOC * (dry_doil_mass * (soil_moisture * voxel_volume)) / 12 / self.ratio_C_per_amino_acid
    
    @state
    def _water_volume(self, soil_moisture, voxel_volume):
        return soil_moisture * voxel_volume
    
    #TP@state
    def _soil_moisture(self, soil_moisture):
        return soil_moisture
    
    @state
    def _water_potential_soil(self, voxel_volume, water_volume):
        """
        Water retention curve from van Genuchten 1980
        """
        m = 1 - (1/self.water_n)
        return - (1 / self.water_alpha) * (
                                            ((self.theta_S - self.theta_R) / ((water_volume / voxel_volume) - self.theta_R)) ** (1 / m) - 1 
                                        )** (1 / self.water_n)
    