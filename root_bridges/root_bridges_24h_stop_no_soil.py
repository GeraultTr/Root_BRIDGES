# Untouched models
from rhizodep.root_carbon import RootCarbonModel
from root_cynaps.root_nitrogen import RootNitrogenModel
from rhizodep.root_growth import RootGrowthModel
from root_cynaps.root_anatomy import RootAnatomy
from root_cynaps.root_water import RootWaterModel

# Utilities
from openalea.metafspm.composite_wrapper import CompositeModel
from openalea.metafspm.component_factory import Choregrapher
from log.visualize import plot_mtg


class RootBRIDGES(CompositeModel):
    """
    Root-BRIDGES model

    Use guideline :
    1. store in a variable Model(g, time_step) to initialize the model, g being an openalea.MTG() object and time_step a time interval in seconds.

    2. print Model.documentation for more information about editable model parameters (optional).

    3. Use Model.scenario(**dict) to pass a set of scenario-specific parameters to the model (optional).

    4. Use Model.run() in a for loop to perform the computations of a time step on the passed MTG File
    """

    def __init__(self, queues_soil_to_plants, queue_plants_to_soil,
                name: str="Plant", time_step: int=3600, coordinates: list=[0, 0, 0], rotation: float=0, translator_path: str = "", **scenario):
        """
        DESCRIPTION
        ----------
        __init__ method of the model. Initializes the thematic modules and link them.

        :param g: the openalea.MTG() instance that will be worked on. It must be representative of a root architecture.
        :param time_step: the resolution time_step of the model in seconds.
        """
        # DECLARE GLOBAL SIMULATION TIME STEP, FOR THE CHOREGRAPHER TO KNOW IF IT HAS TO SUBDIVIDE TIME-STEPS
        self.name = name
        self.coordinates = coordinates
        self.rotation = rotation

        Choregrapher().add_simulation_time_step(time_step)
        self.time = 0

        parameters = scenario["parameters"]
        root_parameters = parameters["root_bridges"]["roots"]
        self.input_tables = scenario["input_tables"]

        # INIT INDIVIDUAL MODULES
        if len(scenario["input_mtg"]) > 0:
            self.root_growth = RootGrowthModel(g=scenario["input_mtg"]["root_mtg_file"], time_step=time_step, **root_parameters)
        else:
            self.root_growth = RootGrowthModel(g=None, time_step=time_step, **root_parameters)
        self.g_root = self.root_growth.g
        self.root_anatomy = RootAnatomy(self.g_root, time_step, **root_parameters)
        self.root_water = RootWaterModel(self.g_root, time_step, **root_parameters)
        self.root_carbon = RootCarbonModel(self.g_root, time_step, **root_parameters)
        self.root_nitrogen = RootNitrogenModel(self.g_root, time_step, **root_parameters)
        
        # LINKING MODULES
        self.declare_data_and_couple_components(root=self.g_root,
                                                translator_path=translator_path,
                                                components=(self.root_growth, self.root_anatomy, self.root_water, self.root_carbon, self.root_nitrogen))
        
        
        # Specific here TODO remove later
        self.root_water.collar_children = self.root_growth.collar_children
        self.root_water.collar_skip = self.root_growth.collar_skip
        self.root_nitrogen.collar_children = self.root_growth.collar_children
        self.root_nitrogen.collar_skip = self.root_growth.collar_skip

        self.reinitialize_step = scenario["target_day"] * 24

        # Provide signature for the MTG
        # Retreive the queues to communicate with environment models
        self.queues_soil_to_plants=queues_soil_to_plants
        self.queue_plants_to_soil=queue_plants_to_soil

        # Get properties from each MTG
        self.root_props = self.g_root.properties()

        self.init_inertials = False

        if self.init_inertials:
            self.root_water_initial_values = {state_var:self.root_props[state_var][1] for state_var in self.root_water.state_variables if state_var not in ("K", "xylem_water")}
            self.root_nitrogen_initial_values = {state_var:self.root_props[state_var][1] for state_var in self.root_nitrogen.state_variables}
            self.root_water_total_initial_values = {state_var:self.root_props[state_var][1] for state_var in self.root_water.plant_scale_state}
            self.root_nitrogen_total_initial_values = {state_var:self.root_props[state_var][1] for state_var in self.root_nitrogen.plant_scale_state}
        
        
        # Performed in initialization and run to update coordinates
        plot_mtg(self.g_root, position=self.coordinates, rotation=self.rotation)

        self.name = name
        # ROOT PROPERTIES INITIAL PASSING IN MTG
        self.root_props["plant_id"] = name
        self.root_props["model_name"] = self.__class__.__name__
        self.root_props["carried_components"] = [component.__class__.__name__ for component in self.components]
        self.queue_plants_to_soil.put({"plant_id": self.name, "data": self.root_props})
        
        # Retreive post environments init states
        self.get_environment_boundaries()

        # Send command to environments models to run first
        self.send_plant_status_to_environment()


    def run(self):

        if self.time <= self.reinitialize_step:
            self.apply_input_tables(tables=self.input_tables, to=self.components, when=self.time)

        if self.time == self.reinitialize_step and self.init_inertials:
            print("Reinitializing Water and Nitrogen for the next 24h of interest")
            for prop, initial_value in self.root_water_initial_values.items():
                self.root_props[prop].update({v: initial_value for v in self.root_water.vertices})
            for prop, initial_value in self.root_nitrogen_initial_values.items():
                self.root_props[prop].update({v: initial_value for v in self.root_nitrogen.vertices})
            for prop, initial_value in self.root_water_total_initial_values.items():
                self.root_props[prop].update({1: initial_value})
            for prop, initial_value in self.root_nitrogen_total_initial_values.items():
                self.root_props[prop].update({1: initial_value})

        if self.time < self.reinitialize_step:

            # Retrieve soil and light status for plant
            self.get_environment_boundaries()
            
            # Compute root growth from resulting states
            self.root_growth(modules_to_update=[c for c in self.components if c.__class__.__name__ != "RootGrowthModel"],
                            soil_boundaries_to_infer=self.soil_outputs)
            
            # Update MTG coordinates accounting for position in the scene
            plot_mtg(self.g_root, position=self.coordinates, rotation=self.rotation)

            # Update topological surfaces and volumes based on other evolved structural properties
            self.root_anatomy()

        # Compute state variations for water and then carbon and nitrogen
        self.root_water()

        if self.time < self.reinitialize_step:
            self.root_carbon()

        self.root_nitrogen()

        # Send plant status to soil and light models
        self.send_plant_status_to_environment()

        self.time += 1


    def get_environment_boundaries(self):
        # Wait for results from both soil and light model before begining
        soil_boundary_props = self.queues_soil_to_plants[self.name].get()

        # NOTE : here you have to perform a per-variable update otherwise dynamic links are broken
        for variable_name in self.soil_outputs + ["voxel_neighbor"]: # TODO : soil_outputs come from declare_data_and_couple_components, not a good structure to keep
            if variable_name not in self.root_props.keys():
                self.root_props[variable_name] = {}
            
            self.root_props[variable_name].update(soil_boundary_props[variable_name])


    def send_plant_status_to_environment(self):
        self.queue_plants_to_soil.put({"plant_id": self.name, "data": self.root_props})