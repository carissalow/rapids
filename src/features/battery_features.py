import pandas as pd
from battery.battery_base import base_battery_features

battery_data = pd.read_csv(snakemake.input[0], parse_dates=["local_start_date_time", "local_end_date_time", "local_start_date", "local_end_date"])
day_segment = snakemake.params["day_segment"]
requested_features = snakemake.params["features"]
battery_features = pd.DataFrame(columns=["local_date"])

battery_features = battery_features.merge(base_battery_features(battery_data, day_segment, requested_features), on="local_date", how="outer")

assert len(requested_features) + 1 == battery_features.shape[1], "The number of features in the output dataframe (=" + str(battery_features.shape[1]) + ") does not match the expected value (=" + str(len(requested_features)) + " + 1). Verify your battery feature extraction functions"

battery_features.to_csv(snakemake.output[0], index=False)
