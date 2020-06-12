import pandas as pd
from fitbit_sleep.fitbit_sleep_base import base_fitbit_sleep_features
import itertools

sleep_summary_data = pd.read_csv(snakemake.input["sleep_summary_data"])
requested_summary_features = snakemake.params["summary_features"]
requested_sleep_type = snakemake.params["sleep_types"]
day_segment = snakemake.params["day_segment"]
sleep_features = pd.DataFrame(columns=["local_date"])

sleep_features = sleep_features.merge(base_fitbit_sleep_features(sleep_summary_data, day_segment, requested_summary_features, requested_sleep_type), on="local_date", how="outer")

requested_features = ["".join(feature) for feature in itertools.product(requested_summary_features, requested_sleep_type)] if day_segment == "daily" else []

assert len(requested_features) + 1 == sleep_features.shape[1], "The number of features in the output dataframe (=" + str(sleep_features.shape[1]) + ") does not match the expected value (=" + str(len(requested_features)) + " + 1). Verify your fitbit sleep feature extraction functions"

sleep_features.to_csv(snakemake.output[0], index=False)

