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

# Add "is_stationary" column to denote whether it is stationary or not
# "distance", "duration_in_seconds" (if it does not exist), and "speed" columns are also added
def mark_as_stationary(location_data, threshold_static):

    if not location_data.index.is_monotonic:
        location_data.sort_index(inplace=True)

    # Distance in meters
    location_data["distance"] = haversine(location_data["double_longitude"], location_data["double_latitude"], location_data["double_longitude"].shift(-1), location_data["double_latitude"].shift(-1))
    # Duration in seconds
    if "duration_in_seconds" not in location_data.columns:
        location_data["duration_in_seconds"] = (location_data.timestamp.diff(-1) * (-1)) / 1000
    # Speed in km/h
    location_data["speed"] = (location_data["distance"] / location_data["duration_in_seconds"]).replace(np.inf, np.nan) * 3.6

    location_data["is_stationary"] = np.where(location_data["speed"] < threshold_static, 1, 0)

    return location_data

# Relabel clusters: -1 denotes the outliers (insignificant or rarely visited locations), 1 denotes the most visited significant location, 2 denotes the 2nd most significant location,...
def label(location_data):

    location_data["count"] = location_data.groupby(["cluster_label"], sort=False)["cluster_label"].transform("count")
    location_data.loc[location_data["cluster_label"] == -1, "count"] = -1
    
    num_clusters = location_data["cluster_label"].nunique()
    cluster_label = num_clusters - (location_data.groupby(["count", "cluster_label"], sort=True).grouper.group_info[0])
    cluster_label[cluster_label == num_clusters] = -1

    location_data["cluster_label"] = cluster_label
    del location_data["count"]

    return location_data

def create_clustering_hyperparameters(clustering_algorithm, dbscan_eps, dbscan_minsamples):
    if clustering_algorithm == "DBSCAN":
        hyperparameters = {"eps": meters_to_degrees(dbscan_eps), "min_samples": dbscan_minsamples}
    else: # OPTICS
        hyperparameters = {"max_eps": meters_to_degrees(dbscan_eps), "min_samples": dbscan_minsamples, "metric": "euclidean", "cluster_method": "dbscan"}
    
    return hyperparameters

# Only stationary samples are clustered, hence moving samples are labeled with NA
def cluster(location_data, clustering_algorithm, threshold_static, **kwargs):

    if location_data.empty:
        return pd.DataFrame(columns=location_data.columns.tolist() + ["is_stationary", "cluster_label"])

    location_data = location_data.set_index("local_date_time")

    location_data = mark_as_stationary(location_data, threshold_static)
    
    # Only keep stationary samples for clustering
    stationary_data = location_data[location_data["is_stationary"] == 1][["double_latitude", "double_longitude"]]

    # Remove duplicates and apply sample_weight (only available for DBSCAN currently) to reduce memory usage
    stationary_data_dedup = stationary_data.groupby(["double_latitude", "double_longitude"]).size().reset_index()
    lat_lon = stationary_data_dedup[["double_latitude", "double_longitude"]].values

    if stationary_data_dedup.shape[0] < kwargs["min_samples"]:
        cluster_results = np.array([-1] * stationary_data_dedup.shape[0])
    elif clustering_algorithm == "DBSCAN":
        clusterer = DBSCAN(**kwargs)
        cluster_results = clusterer.fit_predict(lat_lon, sample_weight=stationary_data_dedup[0])
    else: # OPTICS
        clusterer = OPTICS(**kwargs)
        cluster_results = clusterer.fit_predict(lat_lon)

    # Add cluster labels
    stationary_data_dedup["cluster_label"] = cluster_results
    stationary_data_with_labels = label(stationary_data.reset_index().merge(stationary_data_dedup[["double_latitude", "double_longitude", "cluster_label"]], how="left", on=["double_latitude", "double_longitude"])).set_index("local_date_time")
    location_data["cluster_label"] = stationary_data_with_labels["cluster_label"]

    return location_data
