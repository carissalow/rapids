import pandas as pd
from utils.utils import fetch_provider_features, run_provider_cleaning_script

sensor_data_files = dict(snakemake.input)

provider = snakemake.params["provider"]
provider_key = snakemake.params["provider_key"]
sensor_key = snakemake.params["sensor_key"]

if sensor_key == "all_cleaning_individual" or sensor_key == "all_cleaning_overall":
    # Data cleaning
    sensor_features = run_provider_cleaning_script(provider, provider_key, sensor_key, sensor_data_files)
else:
    # Extract sensor features
    del sensor_data_files["time_segments_labels"]
    time_segments_file = snakemake.input["time_segments_labels"]
    sensor_features = fetch_provider_features(provider, provider_key, sensor_key, sensor_data_files, time_segments_file)

sensor_features.to_csv(snakemake.output[0], index=False)