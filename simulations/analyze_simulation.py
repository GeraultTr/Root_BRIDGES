import os
import numpy as np

# Utility packages
from analyze.analyze import analyze_data, test_output_range
from log.visualize import post_compress_gltf


if __name__ == '__main__':

    # scenarios = ["Drew_1975_low", "Drew_1975_1"]
    #scenarios = ["Drew_1975_low"]
    # scenarios = ["RC_ref_30D_debug"]
    # target_days = [ 5, 7, 10, 20, 30, 40, 50, 60] #, 70, 80, 90, 100]
    # target_days = np.arange(10, 61, 1)
    # target_days = [10, 20, 30, 40, 50, 60]
    target_days = [125]
    # scenarios = [f"RC_ref_{day}D" for day in target_days]
    scenarios = [f"RC_ref_{day}D" + "_images" for day in target_days]
    # scenarios = [f"RC_no_hair_{day}D" for day in target_days]
    # scenarios = [f"RC_ref_{day}D" for day in target_days] + [f"RC_no_hair_{day}D" for day in target_days]

    output_path = "outputs"
    #output_path = "C:/Users/tigerault/OneDrive - agroparistech.fr/Thesis/Sujet/Modelling/saved_scenarios/05-06_hairless_tests"
    # output_path = "C:/Users/tigerault/OneDrive - agroparistech.fr/Thesis/Sujet/Modelling/saved_scenarios/01-06_ISRR 2024"

    # test_output_range(scenarios=scenarios, outputs_dirpath="outputs", test_file_dirpath="inputs/outputs_validation_root_cynaps_V0.xlsx")

    post_compress_gltf(os.path.join(output_path, scenarios[0], "root_images"))

    analyze_data(scenarios=scenarios, outputs_dirpath=output_path, inputs_dirpath="inputs",
                     on_sums=False,
                     on_performance=False,
                     animate_raw_logs=True,
                     target_properties=None
                     )
    
    # In the end put the system to sleep, Windows only
    #os.system("rundll32.exe powrprof.dll,SetSuspendState 0,1,0")
    