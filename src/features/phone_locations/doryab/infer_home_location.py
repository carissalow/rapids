import warnings
import pandas as pd
from doryab_clustering import haversine, create_clustering_hyperparameters, cluster



def infer_home_location(location_data, clustering_algorithm, dbscan_eps, dbscan_minsamples, threshold_static, strategy):
    
    # Home locations are inferred based on records logged during midnight to 6am.
    if strategy == "TIME_CLUSTER":
        
        location_data_filtered = location_data[location_data["local_hour"] <= 6]
    
        if location_data_filtered.empty:
            warnings.warn("We could not infer a home location because there are no location records logged during midnight to 6am.")
            return pd.DataFrame(columns=location_data_filtered.columns.tolist() + ["distance_from_home"])
        
        hyperparameters = create_clustering_hyperparameters(clustering_algorithm, dbscan_eps, dbscan_minsamples)
        location_data_filtered = cluster(location_data_filtered, clustering_algorithm, threshold_static, **hyperparameters)

        home_location = location_data_filtered[location_data_filtered["cluster_label"] == 1][["double_latitude", "double_longitude"]].mean()
        distance_from_home = haversine(location_data["double_longitude"], location_data["double_latitude"], [home_location["double_longitude"]] * location_data.shape[0], [home_location["double_latitude"]] * location_data.shape[0])
        location_data["distance_from_home"] = distance_from_home

    else:
        raise ValueError("It only support TIME_CLUSTER strategy currently.")    
    
    return location_data



location_data = pd.read_csv(snakemake.input["sensor_input"])
accuracy_limit = snakemake.params["accuracy_limit"]
maximum_row_gap = snakemake.params["maximum_row_gap"]
maximum_row_duration = snakemake.params["maximum_row_duration"]
dbscan_eps = snakemake.params["dbscan_eps"]
dbscan_minsamples = snakemake.params["dbscan_minsamples"]
threshold_static = snakemake.params["threshold_static"]
clustering_algorithm = snakemake.params["clustering_algorithm"]
strategy = "TIME_CLUSTER"

rows_before_accuracy_filter = len(location_data)
location_data.query("accuracy < @accuracy_limit and (double_latitude != 0 or double_longitude != 0)", inplace=True)

if rows_before_accuracy_filter > 0 and len(location_data) == 0:
    warnings.warn("Cannot compute Doryab location features because there are no rows with an accuracy value lower than ACCURACY_LIMIT: {}".format(accuracy_limit))

if not location_data.timestamp.is_monotonic:
    location_data.sort_values(by=["timestamp"], inplace=True)

location_data["duration_in_seconds"] = (location_data.timestamp.diff(-1) * (-1)) / 1000
location_data.loc[location_data["duration_in_seconds"] >= maximum_row_gap, "duration_in_seconds"] = maximum_row_duration

location_data_with_home = infer_home_location(location_data, clustering_algorithm, dbscan_eps, dbscan_minsamples, threshold_static, strategy)
location_data_with_home.to_csv(snakemake.output[0], index=False)

