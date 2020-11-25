import pandas as pd
import numpy as np

merge_keys = ["pid", "local_segment", "local_segment_label", "local_segment_start_datetime", "local_segment_end_datetime"]
sensor_features = pd.read_csv(snakemake.input["cleaned_sensor_features"])

all_demographic_features = pd.DataFrame()
for demographic_features_path in snakemake.input["demographic_features"]:
    pid = demographic_features_path.split("/")[3]
    demographic_features = pd.read_csv(demographic_features_path)
    demographic_features = demographic_features.assign(pid=pid)
    all_demographic_features = pd.concat([all_demographic_features, demographic_features], axis=0)
    
# merge sensor features and demographic features
features = sensor_features.merge(all_demographic_features, on="pid", how="left")

all_targets = pd.DataFrame()
for targets_path in snakemake.input["targets"]:
    pid = targets_path.split("/")[3]
    targets = pd.read_csv(targets_path)
    targets = targets.assign(pid=pid)
    all_targets = pd.concat([all_targets, targets], axis=0)

# merge features and targets
data = features.merge(all_targets[["target"] + merge_keys], on=merge_keys, how="inner")

data.to_csv(snakemake.output[0], index=False)
