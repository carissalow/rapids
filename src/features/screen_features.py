import pandas as pd
import itertools
from screen.screen_base import base_screen_features

screen_data = pd.read_csv(snakemake.input["screen_deltas"], parse_dates=["local_start_date_time", "local_end_date_time", "local_start_date", "local_end_date"])
phone_sensed_bins = pd.read_csv(snakemake.input["phone_sensed_bins"], parse_dates=["local_date"], index_col="local_date")
phone_sensed_bins[phone_sensed_bins > 0] = 1
day_segment = snakemake.params["day_segment"]
screen_features = pd.DataFrame(columns=["local_date"])

requested_features_deltas = ["firstuseafter" + "{0:0=2d}".format(snakemake.params["reference_hour_first_use"]) if feature_name == "firstuseafter" else feature_name for feature_name in snakemake.params["features_deltas"]]
requested_features = ["".join(feature) for feature in itertools.product(requested_features_deltas, snakemake.params["episode_types"])]

screen_features = screen_features.merge(base_screen_features(screen_data, phone_sensed_bins, day_segment, snakemake.params), on="local_date", how="outer")

assert len(requested_features) + 1 == screen_features.shape[1], "The number of features in the output dataframe (=" + str(screen_features.shape[1]) + ") does not match the expected value (=" + str(len(requested_features)) + " + 1). Verify your screen feature extraction functions"

screen_features.to_csv(snakemake.output[0], index=False)
