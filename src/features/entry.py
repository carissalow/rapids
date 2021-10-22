import pandas as pd
from utils.utils import fetch_provider_features, run_provider_cleaning_script

sensor_data_files = dict(snakemake.input)

provider = snakemake.params["provider"]
provider_key = snakemake.params["provider_key"]
sensor_key = snakemake.params["sensor_key"]

if "time_segments_labels" in sensor_data_files.keys():
    # Extract sensor features
    del sensor_data_files["time_segments_labels"]
    time_segments_file = snakemake.input["time_segments_labels"]
    sensor_features = fetch_provider_features(provider, provider_key, sensor_key, sensor_data_files, time_segments_file)
else:
    # Data cleaning
    sensor_features = run_provider_cleaning_script(provider, provider_key, sensor_key, sensor_data_files)

sensor_features.to_csv(snakemake.output[0], index=False)