import os

# Utility packages
from analyze.analyze import analyze_data


if __name__ == '__main__':
        
    scenarios = ["RB_ref"]
    for scenario_name in scenarios:
        target_folder_key = "RootBRIDGES_0"
        output_path = os.path.join("outputs", scenario_name, target_folder_key)

        analyze_data(scenarios=scenario_name, outputs_dirpath=output_path, target_folder_key=target_folder_key,
                        inputs_dirpath="inputs",
                        on_sums=True,
                        on_performance=False,
                        animate_raw_logs=True,
                        target_properties=None)
