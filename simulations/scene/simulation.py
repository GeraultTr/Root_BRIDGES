# Public packages
import multiprocessing as mp
import time

# Model packages
import openalea.rootbridges
from openalea.rootbridges import RootBRIDGES
from openalea.rhizosoil.model import RhizoSoil

# Utility packages
from openalea.fspm.utility.scenario.initialize import MakeScenarios as ms
from openalea.fspm.utility.writer.logging import Logger
from openalea.metafspm.scene_wrapper import play_Orchestra
from openalea.fspm.utility.plot import analyze_data


if __name__ == "__main__":
    # scenarios = ms.from_table(file_path="inputs/Scenarios_24_11_10.xlsx", which=["RB_ref"])
    scenarios = ms.from_table(file_path="inputs/Scenarios_25_08_08.xlsx", which=["CN_Wheat_init"])
    custom_output_folder = "outputs/recoupling"

    scene_xrange = 0.15
    scene_yrange = 0.15
    sowing_density = 50
    row_spacing = 0.15
    environment_models_number = 1
    subprocesses_number = int(max(scene_xrange * scene_yrange * sowing_density, 1)) + environment_models_number
    parallel_development = 1 # To keep room in CPUs if launching dev simulations in parallel on the machine
    max_processes = mp.cpu_count() - subprocesses_number - parallel_development - 1 # -1 for the main process

    parallel = False
    active_processes = 0 
    processes = []

    for scenario_name, scenario in scenarios.items():
        if parallel:
            # Main process creation part
            # Wait until there is a free slot
            while True:
                # Remove any finished processes and join them to release resources
                alive_processes = []
                for proc in processes:
                    if proc.is_alive():
                        alive_processes.append(proc)
                    else:
                        proc.join()
                processes = alive_processes

                if len(processes) < max_processes:
                    break
                time.sleep(1)

            print("")
            print(f'### Launching {scenario_name} over already {active_processes} processes running ###')
            print("")
            active_processes += subprocesses_number
                
            p = mp.Process(target=play_Orchestra, kwargs=dict(scene_name=scenario_name, output_folder=custom_output_folder, plant_models=[RootBRIDGES], plant_scenarios=[scenario], 
                                                            soil_model=RhizoSoil, soil_scenario=scenario,
                                                            translator_path=openalea.rootbridges.__path__[0],
                                                            logger_class=Logger, log_settings=Logger.light_log,
                                                            scene_xrange=scene_xrange, scene_yrange=scene_yrange, sowing_density=sowing_density,
                                                            time_step=3600, n_iterations=24))
            
            p.start()
            processes.append(p)

        else:
            play_Orchestra(scene_name=scenario_name, output_folder=custom_output_folder, plant_models=[RootBRIDGES], plant_scenarios=[scenario], 
                                soil_model=RhizoSoil, soil_scenario=scenario,
                                translator_path=openalea.rootbridges.__path__[0],
                                logger_class=Logger, log_settings=Logger.heavy_log, heavy_log_period=1,
                                scene_xrange=scene_xrange, scene_yrange=scene_yrange, sowing_density=sowing_density, row_spacing=row_spacing,
                                time_step=3600, n_iterations=24*80)
            
            target_folder_key = "RootBRIDGES_0"

            analyze_data(scenarios=[scenario_name], outputs_dirpath=custom_output_folder, target_folder_key=target_folder_key,
                            inputs_dirpath="inputs",
                            on_sums=True,
                            on_performance=False,
                            animate_raw_logs=False,
                            target_properties=None
                            )
    
    # After all tasks started, join any remaining processes
    for proc in processes:
        proc.join()
