import pandas as pd
import numpy as np
from fitbit_step.fitbit_step_base import base_fitbit_step_features

step_data = pd.read_csv(snakemake.input["step_data"], parse_dates=["local_date_time", "local_date"])
day_segment = snakemake.params["day_segment"]
threshold_active_bout = snakemake.params["threshold_active_bout"]
include_zero_step_rows = snakemake.params["include_zero_step_rows"]
step_features = pd.DataFrame(columns=["local_date"])

requested_features = {}
requested_features["features_all_steps"] = snakemake.params["features_all_steps"]
requested_features["features_sedentary_bout"] = snakemake.params["features_sedentary_bout"]
requested_features["features_active_bout"] = snakemake.params["features_active_bout"]

step_features = step_features.merge(base_fitbit_step_features(step_data, day_segment, requested_features, threshold_active_bout, include_zero_step_rows), on="local_date", how="outer")


assert np.sum([len(x) for x in requested_features.values()]) + 1 == step_features.shape[1], "The number of features in the output dataframe (=" + str(step_features.shape[1]) + ") does not match the expected value (=" + str(np.sum([len(x) for x in requested_features.values()])) + " + 1). Verify your fitbit step feature extraction functions"

step_features.to_csv(snakemake.output[0], index=False)