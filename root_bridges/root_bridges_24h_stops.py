import os

import root_bridges

# Edited models
from root_bridges.soil_model import SoilModel

# Untouched models
from rhizodep.root_carbon import RootCarbonModel
from root_cynaps.root_nitrogen import RootNitrogenModel
from rhizodep.root_growth import RootGrowthModel
from rhizodep.root_anatomy import RootAnatomy
from root_cynaps.root_water import RootWaterModel

# Utilities
from openalea.metafspm.composite_wrapper import CompositeModel
from openalea.metafspm.component_factory import Choregrapher


class Model(CompositeModel):
    """
    Root-BRIDGES model

    Use guideline :
    1. store in a variable Model(g, time_step) to initialize the model, g being an openalea.MTG() object and time_step an time interval in seconds.

    2. print Model.documentation for more information about editable model parameters (optional).

    3. Use Model.scenario(**dict) to pass a set of scenario-specific parameters to the model (optional).

    4. Use Model.run() in a for loop to perform the computations of a time step on the passed MTG File
    """

    def __init__(self, time_step: int, target_day: int, **scenario):
        """
        DESCRIPTION
        ----------
        __init__ method of the model. Initializes the thematic modules and link them.

        :param g: the openalea.MTG() instance that will be worked on. It must be representative of a root architecture.
        :param time_step: the resolution time_step of the model in seconds.
        """

        # DECLARE GLOBAL SIMULATION TIME STEP
        Choregrapher().add_simulation_time_step(time_step)
        self.time = 0
        parameters = scenario["parameters"]["root_bridges"]["roots"]
        self.input_tables = scenario["input_tables"]

        # INIT INDIVIDUAL MODULES
        if len(scenario["input_mtg"]) > 0:
            self.root_growth = RootGrowthModel(scenario["input_mtg"]["root_mtg_file"], time_step, **parameters)
        else:
            self.root_growth = RootGrowthModel(g=None, time_step=time_step, **parameters)
        self.g = self.root_growth.g
        
        self.root_anatomy = RootAnatomy(self.g, time_step, **parameters)
        self.root_water = RootWaterModel(self.g, time_step, **parameters)
        self.root_carbon = RootCarbonModel(self.g, time_step, **parameters)
        self.root_nitrogen = RootNitrogenModel(self.g, time_step, **parameters)
        self.soil = SoilModel(self.g, time_step, **parameters)
        self.soil_voxels = self.soil.voxels

        self.init_inertials = False

        if self.init_inertials:
            self.root_water_initial_values = {state_var:getattr(self.root_water, state_var)[1] for state_var in self.root_water.state_variables if state_var not in ("K", "xylem_water")}
            self.root_nitrogen_initial_values = {state_var:getattr(self.root_nitrogen, state_var)[1] for state_var in self.root_nitrogen.state_variables}
            self.root_water_total_initial_values = {state_var:getattr(self.root_water, state_var)[1] for state_var in self.root_water.plant_scale_state}
            self.root_nitrogen_total_initial_values = {state_var:getattr(self.root_nitrogen, state_var)[1] for state_var in self.root_nitrogen.plant_scale_state}

        # LINKING MODULES
        self.declare_data_and_couple_components(root=self.g, soil=self.soil_voxels,
                                           translator_path=os.path.join(root_bridges.__path__[0], "coupling_translator_uncoupled"),
                                           components=(self.root_growth, self.root_anatomy, self.root_water, self.root_carbon, self.root_nitrogen, self.soil))

        # Some initialization must be performed after linking modules
        #self.root_water.post_coupling_init()
        self.root_water.collar_children = self.root_growth.collar_children
        self.root_water.collar_skip = self.root_growth.collar_skip
        self.root_nitrogen.collar_children = self.root_growth.collar_children
        self.root_nitrogen.collar_skip = self.root_growth.collar_skip

        self.reinitialize_step = target_day * 24


    def run(self):
        if self.time <= self.reinitialize_step:
            self.apply_input_tables(tables=self.input_tables, to=self.components, when=self.time)
        
        if self.time == self.reinitialize_step and self.init_inertials:
            print("Reinitializing Water and Nitrogen for the next 24h of interest")
            for prop, initial_value in self.root_water_initial_values.items():
                getattr(self.root_water, prop).update({v: initial_value for v in self.root_water.vertices})
            for prop, initial_value in self.root_nitrogen_initial_values.items():
                getattr(self.root_nitrogen, prop).update({v: initial_value for v in self.root_nitrogen.vertices})
            for prop, initial_value in self.root_water_total_initial_values.items():
                getattr(self.root_water, prop).update({1: initial_value})
            for prop, initial_value in self.root_nitrogen_total_initial_values.items():
                getattr(self.root_nitrogen, prop).update({1: initial_value})
                
        if self.time < self.reinitialize_step:
            
            # Compute root growth from resulting states
            self.root_growth(modules_to_update=[c for c in self.components if c.__class__.__name__ != "RootGrowthModel"])

            self.soil.compute_mtg_voxel_neighbors()
            self.soil.get_from_voxel()

            # Update topological surfaces and volumes based on other evolved structural properties
            self.root_anatomy()
        
        # Compute state variations for water and then carbon and nitrogen
        self.root_water()

        if self.time < self.reinitialize_step:
            self.root_carbon()

        self.root_nitrogen()

        if self.time < self.reinitialize_step:
            # Update environment boundary conditions
            self.soil()

        self.time += 1

