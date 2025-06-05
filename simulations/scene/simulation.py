# Model packages
import root_bridges
from wheat_bridges.rhizospheric_soil import RhizosphericSoil
from root_bridges.root_bridges_no_soil import RootBRIDGES

# Utility packages
from initialize.initialize import MakeScenarios as ms
from log.logging import Logger
from metafspm.scene_wrapper import play_Orchestra


if __name__ == "__main__":
    scenarios = ms.from_table(file_path="inputs/Scenarios_24_11_10.xlsx", which=["RB_ref"])
    for scenario_name, scenario in scenarios.items():
        play_Orchestra(scene_name=scenario_name, output_folder="outputs", plant_models=[RootBRIDGES], plant_scenarios=[scenario], 
                            soil_model=RhizosphericSoil, soil_scenario=scenario,
                            logger_class=Logger, log_settings=Logger.light_log,
                            translator_path=root_bridges.__path__[0],
                            scene_xrange=0.15, scene_yrange=0.15, sowing_density=1,
                            n_iterations=50*24)
