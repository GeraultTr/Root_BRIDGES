import os

import root_bridges

# Edited models
from root_bridges.soil_model import SoilModel

# Untouched models
from rhizodep.root_carbon import RootCarbonModel
from root_cynaps.root_cynaps import RootNitrogenModel
from rhizodep.root_growth import RootGrowthModel
from rhizodep.root_anatomy import RootAnatomy
from root_cynaps.root_water import RootWaterModel

# Utilities
from metafspm.composite_wrapper import CompositeModel
from metafspm.component_factory import Choregrapher


class Model(CompositeModel):
    """
    Root-BRIDGES model

    Use guideline :
    1. store in a variable Model(g, time_step) to initialize the model, g being an openalea.MTG() object and time_step an time interval in seconds.

    2. print Model.documentation for more information about editable model parameters (optional).

    3. Use Model.scenario(**dict) to pass a set of scenario-specific parameters to the model (optional).

    4. Use Model.run() in a for loop to perform the computations of a time step on the passed MTG File
    """

    def __init__(self, time_step: int, **scenario):
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
        
        # LINKING MODULES
        self.declare_data_and_couple_components(root=self.g, soil=self.soil_voxels,
                                           translator_path=os.path.join(root_bridges.__path__[0], "coupling_translator_uncoupled"),
                                           components=(self.root_growth, self.root_anatomy, self.root_water, self.root_carbon, self.root_nitrogen, self.soil))

        # Some initialization must be performed after linking modules
        self.root_water.post_coupling_init()


    def run(self):
        self.apply_input_tables(tables=self.input_tables, to=self.components, when=self.time)

        # Update environment boundary conditions
        self.soil()

        # Compute root growth from resulting states
        self.root_growth(modules_to_update=[c for c in self.components if c.__class__.__name__ != "RootGrowthModel"])

        # Extend property dictionaries after growth
        self.soil.post_growth_updating()
        
        # Update topological surfaces and volumes based on other evolved structural properties
        self.root_anatomy()

        # Compute state variations for water and then carbon and nitrogen
        self.root_water()
        self.root_carbon()
        self.root_nitrogen()

        self.time += 1

