import pandas as pd
from light.light_base import base_light_features

light_data = pd.read_csv(snakemake.input[0], parse_dates=["local_date_time", "local_date"])
day_segment = snakemake.params["day_segment"]
metrics = snakemake.params["metrics"]
light_features = pd.DataFrame(columns=["local_date"])

light_features = light_features.merge(base_light_features(light_data, day_segment, metrics), on="local_date", how="outer")

assert len(metrics) + 1 == light_features.shape[1], "The number of features in the output dataframe (=" + str(light_features.shape[1]) + ") does not match the expected value (=" + str(len(metrics)) + " + 1). Verify your light feature extraction functions"

light_features.to_csv(snakemake.output[0], index=False)