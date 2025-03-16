# Public packages
import os, sys
# Model packages
from root_bridges.root_bridges_uncoupled import Model
# Utility packages
from log.logging import Logger
from analyze.analyze import analyze_data, test_output_range
from initialize.initialize import MakeScenarios as ms


def single_run(scenario, outputs_dirpath="outputs", simulation_length=2500, echo=True, log_settings={}, analyze=False):
    root_bridges = Model(time_step=3600, **scenario)

    logger = Logger(model_instance=root_bridges, components=root_bridges.components,
                    outputs_dirpath=outputs_dirpath, 
                    time_step_in_hours=1, logging_period_in_hours=24,
                    recording_shoot=False,
                    echo=echo, **log_settings)
    
    stop_file = os.path.join(outputs_dirpath, "Delete_to_Stop")
    open(stop_file, "w").close()

    try:
        for _ in range(simulation_length):
            # Placed here also to capture mtg initialization
            #logger()
            logger.run_and_monitor_model_step()
            #root_bridges.run()

            if not os.path.exists(stop_file):
                raise KeyboardInterrupt

    except (ZeroDivisionError, KeyboardInterrupt):
        logger.exceptions.append(sys.exc_info())

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
    scenarios = ms.from_table(file_path="inputs/Scenarios_24_11_10.xlsx", which=["RC_ref"])
    # scenarios = ms.from_table(file_path="inputs/Scenarios_24_11_10.xlsx", which=["RC_debug"])
    simulate_scenarios(scenarios, simulation_length=24*10, custom_prefix="10D", log_settings=Logger.light_log)
    