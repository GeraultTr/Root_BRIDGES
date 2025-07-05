# Public packages
import os
import numpy as np
import multiprocessing as mp
import time

# Model packages
import root_bridges
from wheat_bridges.rhizospheric_soil import RhizosphericSoil
from root_bridges.root_bridges_24h_stop_no_soil import RootBRIDGES

# Utility packages
from initialize.initialize import MakeScenarios as ms
from log.logging import Logger
from openalea.metafspm.scene_wrapper import play_Orchestra


if __name__ == '__main__':
    # scenarios = ms.from_table(file_path="inputs/Scenarios_24_11_10.xlsx", which=["RC_ref_big_lats"])
    # scenarios = ms.from_table(file_path="inputs/Scenarios_24_11_10.xlsx", which=["RC_ref_0.1",	"RC_ref_0.01",	"RC_ref_0.05",	"RC_ref_0.5",	"RC_ref_5",	"RC_ref_50"])
    scenarios = ms.from_table(file_path="inputs/Scenarios_24_11_10.xlsx", which=["RC_ref_w"])
    custom_output_folder = "outputs/fig_7.3"

    scene_xrange = 0.15
    scene_yrange = 0.15
    sowing_density = 1
    environment_models_number = 1
    subprocesses_number = int(max(scene_xrange * scene_yrange * sowing_density, 1)) + environment_models_number
    parallel_development = 1 # To keep room in CPUs if launching dev simulations in parallel on the machine
    max_processes = mp.cpu_count() - (subprocesses_number + 1) * (parallel_development + 1) - 1 # -1 for the main process

    # target_days = np.arange(10, 61, 1)
    target_days = np.arange(10, 61, 10)
    target_days = [30]
    # target_concentrations = np.logspace(0, 4, len(target_days)) * 5e-3
    # target_concentrations = np.logspace(0, 4, 11) * 5e-3
    target_concentrations = np.logspace(0, 4, 5) * 5e-3
    target_concentrations = [5e-1]    

    parallel = False

    for scenario_name, scenario in scenarios.items():

        # scenarios = ms.from_table(file_path="inputs/Scenarios_24_11_10.xlsx", which=["RC_debug"])
        # target_days = [10, 20, 30, 40, 50, 60]
        
        # target_days = [125]
        # target_days = [10, 20, 30]
        static_days = 1

        if parallel:
            processes = []
            
            for day in target_days:
                for concentration in target_concentrations:
                    scenario["parameters"]["root_bridges"]["roots"]["dissolved_mineral_N"] = 5e-7 * concentration / 1e-1

                    # Main process creation part
                    while len(processes) * (subprocesses_number + 1) >= max_processes:
                        for proc in processes:
                            if not proc.is_alive():
                                processes.remove(proc)
                        time.sleep(1)
                        
                    current_scenario_name = f"{str(scenario_name)}_{concentration:.2e}_{day}D"

                    scenario["target_day"] = day

                    p = mp.Process(target=play_Orchestra, kwargs=dict(scene_name=current_scenario_name, output_folder=custom_output_folder, plant_models=[RootBRIDGES], plant_scenarios=[scenario], 
                                                                    soil_model=RhizosphericSoil, soil_scenario=scenario,
                                                                    translator_path=os.path.join(root_bridges.__path__[0], "coupling_translator_uncoupled"),
                                                                    logger_class=Logger, log_settings=Logger.light_log,
                                                                    scene_xrange=scene_xrange, scene_yrange=scene_yrange, sowing_density=sowing_density,
                                                                    n_iterations=(day + static_days) * 24))
                    
                    p.start()
                    processes.append(p)

        else:
            for day in target_days:
                for concentration in target_concentrations:
                    scenario["parameters"]["root_bridges"]["roots"]["dissolved_mineral_N"] = 5e-7 * concentration / 1e-1
                    
                    current_scenario_name = f"{str(scenario_name)}_{concentration:.2e}_{day}D"

                    scenario["target_day"] = day

                    play_Orchestra(scene_name=current_scenario_name, output_folder=custom_output_folder, plant_models=[RootBRIDGES], plant_scenarios=[scenario], 
                                        soil_model=RhizosphericSoil, soil_scenario=scenario,
                                        translator_path=os.path.join(root_bridges.__path__[0], "coupling_translator_uncoupled"),
                                        logger_class=Logger, log_settings=Logger.light_log,
                                        scene_xrange=scene_xrange, scene_yrange=scene_yrange, sowing_density=sowing_density,
                                        n_iterations=(day + static_days) * 24) 
