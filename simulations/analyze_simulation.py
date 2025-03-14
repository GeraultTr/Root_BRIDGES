import os
# Utility packages
from analyze.analyze import analyze_data, test_output_range


if __name__ == '__main__':

    # scenarios = ["Drew_1975_low", "Drew_1975_1"]
    #scenarios = ["Drew_1975_low"]
    scenarios = ["RC_ref"]

    output_path = "outputs"
    #output_path = "C:/Users/tigerault/OneDrive - agroparistech.fr/Thesis/Sujet/Modelling/saved_scenarios/05-06_hairless_tests"
    # output_path = "C:/Users/tigerault/OneDrive - agroparistech.fr/Thesis/Sujet/Modelling/saved_scenarios/01-06_ISRR 2024"

    # test_output_range(scenarios=scenarios, outputs_dirpath="outputs", test_file_dirpath="inputs/outputs_validation_root_cynaps_V0.xlsx")

    analyze_data(scenarios=scenarios, outputs_dirpath=output_path, inputs_dirpath="inputs",
                     on_sums=True,
                     on_performance=False,
                     animate_raw_logs=False,
                     target_properties=None
                     )
    # In the end put the system to sleep, Windows only
    #os.system("rundll32.exe powrprof.dll,SetSuspendState 0,1,0")
    