import numpy as np
import pandas as pd
from accelerometer.accelerometer_base import base_accelerometer_features


acc_data = pd.read_csv(snakemake.input[0], parse_dates=["local_date_time", "local_date"])
day_segment = snakemake.params["day_segment"]

requested_features = {}
requested_features["magnitude"] = snakemake.params["magnitude"]
requested_features["exertional_activity_episode"] = snakemake.params["exertional_activity_episode"]
requested_features["nonexertional_activity_episode"] = snakemake.params["nonexertional_activity_episode"]

valid_sensed_minutes = snakemake.params["valid_sensed_minutes"]

acc_features = pd.DataFrame(columns=["local_date"])

acc_features = acc_features.merge(base_accelerometer_features(acc_data, day_segment, requested_features, valid_sensed_minutes), on="local_date", how="outer")

assert np.sum([len(x) for x in requested_features.values()]) + (1 if valid_sensed_minutes else 0) + 1 == acc_features.shape[1], "The number of features in the output dataframe (=" + str(acc_features.shape[1]) + ") does not match the expected value (=" + str(np.sum([len(x) for x in requested_features.values()]) + (1 if valid_sensed_minutes else 0)) + " + 1). Verify your accelerometer feature extraction functions"

acc_features.to_csv(snakemake.output[0], index=False)