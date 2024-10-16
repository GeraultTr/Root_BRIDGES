import os
# Utility packages
from analyze.analyze import analyze_data


if __name__ == '__main__':

    scenarios = ["Drew_1975_low", "Drew_1975_1"]
    #scenarios = ["Drew_1975_low"]

    #output_path = "outputs"
    #output_path = "C:/Users/tigerault/OneDrive - agroparistech.fr/Thesis/Sujet/Modelling/saved_scenarios/05-06_hairless_tests"
    output_path = "C:/Users/tigerault/OneDrive - agroparistech.fr/Thesis/Sujet/Modelling/saved_scenarios/01-06_ISRR 2024"

    analyze_data(scenarios=scenarios, outputs_dirpath=output_path,
                     on_sums=False,
                     on_performance=True,
                     animate_raw_logs=True,
                     target_properties=None
                     )
    # In the end put the system to sleep, Windows only
    #os.system("rundll32.exe powrprof.dll,SetSuspendState 0,1,0")
    