import pandas as pd
from location_doryab.location_base import base_location_features

location_data = pd.read_csv(snakemake.input[0], parse_dates=["local_date_time", "local_date"])
day_segment = snakemake.params["day_segment"]
requested_features = snakemake.params["features"] 
location_features = pd.DataFrame(columns=["local_date"])
dbscan_eps = snakemake.params["dbscan_eps"]
dbscan_minsamples = snakemake.params["dbscan_minsamples"]
threshold_static = snakemake.params["threshold_static"]
maximum_gap_allowed = snakemake.params["maximum_gap_allowed"]
minutes_data_used = snakemake.params["minutes_data_used"]

if(minutes_data_used):
        requested_features.append("minutesdataused")

base_features = base_location_features(location_data, day_segment, requested_features, dbscan_eps, dbscan_minsamples,threshold_static,maximum_gap_allowed)

location_features = location_features.merge(base_features, on="local_date", how="outer")

assert len(requested_features) + 1 == location_features.shape[1], "The number of features in the output dataframe (=" + str(location_features.shape[1]) + ") does not match the expected value (=" + str(len(requested_features)) + " + 1). Verify your location feature extraction functions"

location_features.to_csv(snakemake.output[0], index=False)