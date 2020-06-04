import pandas as pd
from ar.ar_base import base_ar_features

ar_data = pd.read_csv(snakemake.input[0],parse_dates=["local_date_time"])
ar_deltas = pd.read_csv(snakemake.input[1],parse_dates=["local_start_date_time", "local_end_date_time", "local_start_date", "local_end_date"])
day_segment = snakemake.params["segment"]
requested_features = snakemake.params["features"]
ar_features = pd.DataFrame(columns=["local_date"])


ar_features = ar_features.merge(base_ar_features(ar_data, ar_deltas, day_segment, requested_features), on="local_date", how="outer")

assert len(requested_features) + 1 == ar_features.shape[1], "The number of features in the output dataframe (=" + str(ar_features.shape[1]) + ") does not match the expected value (=" + str(len(requested_features)) + " + 1). Verify your activity recognition feature extraction functions"

ar_features.to_csv(snakemake.output[0], index=False)