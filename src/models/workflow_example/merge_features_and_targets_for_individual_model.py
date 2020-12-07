import pandas as pd
import numpy as np

index_columns = ["local_segment", "local_segment_label", "local_segment_start_datetime", "local_segment_end_datetime"]
sensor_features = pd.read_csv(snakemake.input["cleaned_sensor_features"], index_col=index_columns)
targets = pd.read_csv(snakemake.input["targets"], index_col=index_columns)

data = pd.concat([sensor_features, targets[["target"]]], axis=1, join="inner")

data.to_csv(snakemake.output[0], index=True)
