import pandas as pd
from fitbit_heartrate.fitbit_heartrate_base import base_fitbit_heartrate_features

heartrate_summary_data = pd.read_csv(snakemake.input["heartrate_summary_data"], index_col=["local_date"], parse_dates=["local_date"])
heartrate_intraday_data = pd.read_csv(snakemake.input["heartrate_intraday_data"], parse_dates=["local_date_time", "local_date"])
day_segment = snakemake.params["day_segment"]
requested_summary_features = snakemake.params["summary_features"]
requested_intraday_features = snakemake.params["intraday_features"]
heartrate_features = pd.DataFrame(columns=["local_date"])

heartrate_features = heartrate_features.merge(base_fitbit_heartrate_features(heartrate_summary_data, heartrate_intraday_data, day_segment, requested_summary_features, requested_intraday_features), on="local_date", how="outer")

requested_features = requested_summary_features + requested_intraday_features if day_segment == "daily" else requested_intraday_features
assert len(requested_features) + 1 == heartrate_features.shape[1], "The number of features in the output dataframe (=" + str(heartrate_features.shape[1]) + ") does not match the expected value (=" + str(len(requested_features)) + " + 1). Verify your fitbit heartrate feature extraction functions"

heartrate_features.to_csv(snakemake.output[0], index=False)
