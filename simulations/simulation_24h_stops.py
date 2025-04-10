# Public packages
import os, traceback, time
import multiprocessing as mp
import numpy as np
# Model packages
from root_bridges.root_bridges_24h_stops import Model
# Utility packages
from log.logging import Logger
from analyze.analyze import analyze_data, test_output_range
from initialize.initialize import MakeScenarios as ms
from metafspm.component_factory import Choregrapher


def single_run(scenario, outputs_dirpath="outputs", simulation_length=2500, echo=True, log_settings={}, analyze=False, target_day=1):
    root_bridges = Model(time_step=3600, target_day=target_day, **scenario)

    logger = Logger(model_instance=root_bridges, components=root_bridges.components,
                    outputs_dirpath=outputs_dirpath, 
                    time_step_in_hours=1, logging_period_in_hours=24,
                    recording_shoot=False,
                    echo=echo, **log_settings)
    
    stop_file = os.path.join(outputs_dirpath + " *", "Delete_to_Stop")
    open(stop_file, "w").close()

    try:
        for _ in range(simulation_length):
            # Placed here also to capture mtg initialization
            logger()
            # logger.run_and_monitor_model_step()
            root_bridges.run()

            if not os.path.exists(stop_file):
                raise KeyboardInterrupt

    except Exception as e:
        logger.exceptions.append(traceback.format_exc())

    finally:
        logger.stop()
        if analyze:
            analyze_data(scenarios=[os.path.basename(outputs_dirpath)], outputs_dirpath=outputs_dirpath, target_properties=None, **log_settings)


def simulate_scenarios(scenarios, simulation_length=2500, echo=True, custom_prefix=None, log_settings={}):
    for scenario_name, scenario in scenarios.items():
        # Enable quick parallel testing with exact same parameters
        if custom_prefix:
            scenario_name = f"{scenario_name}_{custom_prefix}"
        
        single_run(scenario, outputs_dirpath=os.path.join("outputs", str(scenario_name)),
                                                      simulation_length=simulation_length,
                                                      echo=echo, log_settings=log_settings)
        
        # test_output_range(scenarios=[scenario_name], outputs_dirpath="outputs", test_file_dirpath="inputs/outputs_validation_root_cynaps_V0.xlsx")

        analyze_data(scenarios=[scenario_name], outputs_dirpath="outputs", inputs_dirpath="inputs",
                     on_sums=True,
                     on_performance=False,
                     animate_raw_logs=False,
                     target_properties=None
                     )


if __name__ == '__main__':
    # scenarios = ms.from_table(file_path="inputs/Scenarios_24_11_10.xlsx", which=["Rhizodep_ref"])
    scenarios = ms.from_table(file_path="inputs/Scenarios_24_11_10.xlsx", which=["RC_ref"])
    # scenarios = ms.from_table(file_path="inputs/Scenarios_24_11_10.xlsx", which=["RC_ref_low"])
    # scenarios = ms.from_table(file_path="inputs/Scenarios_24_11_10.xlsx", which=["RC_ref_high"])
    # scenarios = ms.from_table(file_path="inputs/Scenarios_24_11_10.xlsx", which=["RC_no_hair"])
    scenario_name = list(scenarios.keys())[0]
    scenario = list(scenarios.values())[0]

    # scenarios = ms.from_table(file_path="inputs/Scenarios_24_11_10.xlsx", which=["RC_debug"])
    # target_days = [ 5, 7, 10, 20, 30, 40, 50, 60]
    # target_days = np.arange(10, 61, 1)
    # target_days = [125]
    target_days = [40]

    processes = []
    max_processes = mp.cpu_count()
    for day in target_days:

        while len(processes) == max_processes:
            for proc in processes:
                if not proc.is_alive():
                    processes.remove(proc)
            time.sleep(1)

        current_scenario_name = f"{str(scenario_name)}_{day}D"

        p = mp.Process(target=single_run, kwargs=dict(scenario=scenario, 
                                                      outputs_dirpath=os.path.join("outputs", current_scenario_name),
                                                      target_day=day, simulation_length=(day + 1) * 24,
                                                      echo=True,
                                                      log_settings=Logger.light_log))
        p.start()
        processes.append(p)
2