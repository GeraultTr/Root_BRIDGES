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
    # N related
    C_mineralN_soil: float = declare(default=2.2, unit="mol.m-3", unit_comment="of equivalent mineral nitrogen", description="Mineral nitrogen concentration in soil", 
                                        value_comment="", references="Fischer et al. 1966", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="intensive", edit_by="user")
    content_mineralN_soil: float = declare(default=20e-6 / 14, unit="mol.g-1", unit_comment="of equivalent mineral nitrogen", description="Mineral nitrogen concentration in soil", 
                                        value_comment="", references="Fischer et al. 1966", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="intensive", edit_by="user")
    C_amino_acids_soil: float = declare(default=8.2e-3, unit="mol.m-3", unit_comment="of equivalent mineral nitrogen", description="Mineral nitrogen concentration in soil", 
                                        value_comment="", references="Fischer et al 2007, water leaching estimation", DOI="",
                                       min_value="", max_value="", variable_type="state_variable", by="model_soil", state_variable_type="intensive", edit_by="user")
    content_amino_acids_soil: float = declare(default= 5 * 0.2e-6 /14, unit="mol.m-3", unit_comment="of equivalent mineral nitrogen", description="Mineral nitrogen concentration in soil", 
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

    # PARAMETERS

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
    
    content_mineralN_patch: float = declare(default=20e-6 / 14, unit="mol.g-1", unit_comment="of equivalent mineral nitrogen", description="Mineral nitrogen concentration in a located patch in soil", 
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
        self.add_patch_repartition_to_soil(property_name="content_mineralN_soil", patch_value=self.content_mineralN_patch, 
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
        self.voxels["C_mineralN_soil"] = self.voxels["content_mineralN_soil"] * self.voxels["dry_soil_mass"] / self.voxels["water_volume"]
        self.voxels["C_amino_acids_soil"] = self.voxels["content_amino_acids_soil"] * self.voxels["dry_soil_mass"] / self.voxels["water_volume"]
    
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

    # MODEL EAQUATIONS

    #TP@state
    def _C_mineralN_soil(self, C_mineralN_soil, soil_moisture, voxel_volume, mineralN_diffusion_from_roots, mineralN_diffusion_from_xylem, mineralN_uptake):
        balance = C_mineralN_soil + (self.time_step_in_seconds / (soil_moisture * voxel_volume)) * (
            mineralN_diffusion_from_roots
            + mineralN_diffusion_from_xylem
            - mineralN_uptake
        )
        balance[balance < 0.] = 0.
        return balance

    #TP@state
    def _C_amino_acids_soil(self, C_amino_acids_soil, soil_moisture, voxel_volume, amino_acids_diffusion_from_roots, amino_acids_diffusion_from_xylem, amino_acids_uptake):
        balance = C_amino_acids_soil + (self.time_step_in_seconds / (soil_moisture * voxel_volume)) * (
            amino_acids_diffusion_from_roots
            + amino_acids_diffusion_from_xylem
            - amino_acids_uptake
        )
        balance[balance < 0.] = 0.
        return balance
    
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
    