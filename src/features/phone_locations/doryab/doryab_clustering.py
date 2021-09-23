import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN, OPTICS



# Calculate the great-circle distance (in meters) between two points on the earth (specified in decimal degrees)
def haversine(lon1, lat1, lon2, lat2):
    # Radius of earth in kilometers. Use 3956 for miles
    r = 6371
    # Convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = np.radians([lon1, lat1, lon2, lat2])
    # Haversine formula
    distance = r * 2 * np.arcsin(np.sqrt(np.sin((lat2 - lat1) / 2.0) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2.0) ** 2)) * 1000
    return distance

# Just an approximation, but speeds up clustering by a huge amount and doesn't introduce much error over small distances
# Reference: https://jonisalonen.com/2014/computing-distance-between-coordinates-can-be-simple-and-fast/
def meters_to_degrees(distance):
    # Convert meter to nautical mile
    distance = distance / 1852
    # Convert nautical mile to degree
    distance = distance / 60
    return distance

# Relabel clusters: -1 denotes the outliers (insignificant or rarely visited locations), 1 denotes the most visited significant location, 2 denotes the 2nd most significant location,...
def label(location_data):

    # Exclude outliers (cluster_label = -1) while calculating the total duration of locations in a cluster
    label2duration = location_data[["cluster_label", "duration"]].replace(-1, np.nan).groupby("cluster_label")[["duration"]].sum().sort_values(by=["duration"], ascending=False)
    # Add the row number as the new cluster label
    label2duration["new_cluster_label"] = np.arange(len(label2duration)) + 1
    # Still use -1 to denote the outliers
    label2duration.loc[-1, "new_cluster_label"] = -1
    # Merge the new cluster label with the original location data
    location_data = location_data.merge(label2duration[["new_cluster_label"]], left_on="cluster_label", right_index=True, how="left")

    del location_data["cluster_label"]
    location_data.rename(columns={"new_cluster_label": "cluster_label"}, inplace=True)

    return location_data

def create_clustering_hyperparameters(clustering_algorithm, dbscan_eps, dbscan_minsamples):
    if clustering_algorithm == "DBSCAN":
        hyperparameters = {"eps": meters_to_degrees(dbscan_eps), "min_samples": dbscan_minsamples}
    else: # OPTICS
        hyperparameters = {"max_eps": meters_to_degrees(dbscan_eps), "min_samples": dbscan_minsamples, "metric": "euclidean", "cluster_method": "dbscan"}
    
    return hyperparameters

# Only stationary samples are clustered, hence moving samples are labeled with NA
def cluster(location_data, clustering_algorithm, **kwargs):

    if location_data.empty:
        return pd.DataFrame(columns=location_data.columns.tolist() + ["is_stationary", "cluster_label"])
    
    if "duration" not in location_data.columns:
        # Convert second to minute
        location_data = location_data.assign(duration=location_data["duration_in_seconds"] / 60)

    # Only keep stationary samples for clustering
    stationary_data = location_data[location_data["is_stationary"] == 1][["double_latitude", "double_longitude", "duration"]]

    # Remove duplicates and apply sample_weight (only available for DBSCAN currently) to reduce memory usage
    stationary_data_dedup = stationary_data.groupby(["double_latitude", "double_longitude"])[["duration"]].sum().reset_index()
    lat_lon_dedup = stationary_data_dedup[["double_latitude", "double_longitude"]].values

    if stationary_data_dedup.shape[0] < kwargs["min_samples"]:
        cluster_results = np.array([-1] * stationary_data_dedup.shape[0])
    elif clustering_algorithm == "DBSCAN":        
        clusterer = DBSCAN(**kwargs)
        cluster_results = clusterer.fit_predict(lat_lon_dedup, sample_weight=stationary_data_dedup["duration"])
    else: # OPTICS
        clusterer = OPTICS(**kwargs)
        cluster_results = clusterer.fit_predict(lat_lon_dedup)

    # Add cluster labels
    stationary_data_dedup["cluster_label"] = cluster_results
    location_data_with_labels = label(location_data.merge(stationary_data_dedup[["double_latitude", "double_longitude", "cluster_label"]], how="left", on=["double_latitude", "double_longitude"]))

    return location_data_with_labels
