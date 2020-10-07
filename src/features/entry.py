import pandas as pd
from utils.utils import fetch_provider_features

sensor_data_files = dict(snakemake.input)
del sensor_data_files["day_segments_labels"]
day_segments_file = snakemake.input["day_segments_labels"]

provider = snakemake.params["provider"]
provider_key = snakemake.params["provider_key"]
sensor_key = snakemake.params["sensor_key"]

sensor_features = fetch_provider_features(provider, provider_key, sensor_key, sensor_data_files, day_segments_file)

sensor_features.to_csv(snakemake.output[0], index=False)