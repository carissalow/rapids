from datetime import time
import numpy as np
import pandas as pd
from phone_locations.doryab.doryab_clustering import haversine, create_clustering_hyperparameters, cluster



def apply_cluster_strategy(location_data, time_segment, clustering_algorithm, dbscan_eps, dbscan_minsamples, cluster_on, filter_data_by_segment):

    hyperparameters = create_clustering_hyperparameters(clustering_algorithm, dbscan_eps, dbscan_minsamples)

    if cluster_on == "PARTICIPANT_DATASET":
        # clusters are created in add_doryab_extra_columns.py script
        location_data = filter_data_by_segment(location_data, time_segment)
    elif cluster_on == "TIME_SEGMENT":
        location_data = filter_data_by_segment(location_data, time_segment)
        location_data = cluster(location_data, clustering_algorithm, **hyperparameters)
    else: # TIME_SEGMENT_INSTANCE
        location_data = filter_data_by_segment(location_data, time_segment)
        location_data_clusters = pd.DataFrame()
        for segment_instance, instance_data in location_data.groupby(["local_segment"]):
            location_data_per_group = cluster(instance_data, clustering_algorithm, **hyperparameters)
            location_data_clusters = pd.concat([location_data_per_group, location_data_clusters])
        location_data = location_data_clusters

    return location_data

def variance_and_logvariance_features(location_data, location_features):
    location_data_grouped = location_data.groupby("local_segment")
    location_data["latitude_X_duration"] = location_data["double_latitude"] * location_data["duration"]
    location_data["longitude_X_duration"] = location_data["double_longitude"] * location_data["duration"]
    
    location_data["latitude_wavg"] = location_data_grouped["latitude_X_duration"].transform("sum") / location_data_grouped["duration"].transform("sum")
    location_data["longitude_wavg"] = location_data_grouped["longitude_X_duration"].transform("sum") / location_data_grouped["duration"].transform("sum")

    location_data["latitude_for_wvar"] = (location_data["double_latitude"] - location_data["latitude_wavg"]) ** 2 * location_data["duration"] * 60
    location_data["longitude_for_wvar"] = (location_data["double_longitude"] - location_data["longitude_wavg"]) ** 2 * location_data["duration"] * 60

    location_features["locationvariance"] = ((location_data_grouped["latitude_for_wvar"].sum() + location_data_grouped["longitude_for_wvar"].sum()) / (location_data_grouped["duration"].sum() * 60 - 1)).fillna(0)
    location_features["loglocationvariance"] = np.log10(location_features["locationvariance"]).replace(-np.inf, np.nan)

    return location_features

def distance_and_speed_features(moving_data):

    distance_and_speed = moving_data[["local_segment", "distance"]].groupby(["local_segment"]).sum().rename(columns={"distance": "totaldistance"})
    
    moving_data_grouped = moving_data.groupby(["local_segment"])

    moving_data["speed_X_duration"] = moving_data["speed"] * moving_data["duration"]
    distance_and_speed["avgspeed"] = moving_data_grouped["speed_X_duration"].sum() / moving_data_grouped["duration"].sum()

    moving_data["speed_wavg"] = moving_data_grouped["speed_X_duration"].transform("sum") / moving_data_grouped["duration"].transform("sum")
    moving_data["speed_for_wvar"] = (moving_data["speed"] - moving_data["speed_wavg"]) ** 2 * moving_data["duration"] * 60
    distance_and_speed["varspeed"] = moving_data_grouped["speed_for_wvar"].sum() / (moving_data_grouped["duration"].sum() * 60 - 1)
    
    return distance_and_speed

def radius_of_gyration(location_data):
    
    if location_data.empty:
        return np.nan

    # define a lambda function to compute the weighted mean for each cluster
    weighted_mean = lambda x: np.average(x, weights=location_data.loc[x.index, "duration"])
 
    # center is the centroid of the places visited during a segment instance, not the home location
    clusters = location_data.groupby(["local_segment", "cluster_label"]).agg(
        double_latitude=("double_latitude", weighted_mean),
        double_longitude=("double_longitude", weighted_mean),
        time_in_a_cluster=("duration", "sum")
    ).reset_index()

    # redefine the lambda function to compute the weighted mean across clusters
    weighted_mean = lambda x: np.average(x, weights=clusters.loc[x.index, "time_in_a_cluster"])

    clusters[["centroid_double_latitude", "centroid_double_longitude"]] = clusters.groupby(["local_segment"], sort=False)[["double_latitude", "double_longitude"]].transform(weighted_mean)
    clusters["distance_squared"] = haversine(clusters["double_longitude"], clusters["double_latitude"], clusters["centroid_double_longitude"], clusters["centroid_double_latitude"]) ** 2
    
    clusters["distance_squared_X_time_in_a_cluster"] = clusters["distance_squared"] * clusters["time_in_a_cluster"]
    rog = np.sqrt(clusters.groupby(["local_segment"])["distance_squared_X_time_in_a_cluster"].sum() / clusters.groupby(["local_segment"])["time_in_a_cluster"].sum().replace(0, np.inf))
   
    return rog

def cluster_stay(x, stay_at_clusters, cluster_n):
    topn_cluster_label = x[stay_at_clusters.loc[x.index]["cluster_label"] == cluster_n]
    time_at_topn = topn_cluster_label.iloc[0] if len(topn_cluster_label) == 1 else None
    return time_at_topn

def stay_at_topn_clusters(location_data):

    stay_at_clusters = location_data[["local_segment", "cluster_label", "duration"]].groupby(["local_segment", "cluster_label"], sort=True).sum().reset_index()

    stay_at_clusters_features = stay_at_clusters.groupby(["local_segment"]).agg(        
        timeattop1location=("duration", lambda x: cluster_stay(x, stay_at_clusters, 1)),
        timeattop2location=("duration", lambda x: cluster_stay(x, stay_at_clusters, 2)),
        timeattop3location=("duration", lambda x: cluster_stay(x, stay_at_clusters, 3)),
        maxlengthstayatclusters=("duration", "max"),
        minlengthstayatclusters=("duration", "min"),
        avglengthstayatclusters=("duration", "mean"),
        stdlengthstayatclusters=("duration", "std")
    ).fillna(0)

    return stay_at_clusters_features

def location_entropy(location_data):

    location_data = location_data.groupby(["local_segment", "cluster_label"])[["duration"]].sum().reset_index().rename(columns={"duration": "cluster_duration"})
    location_data["all_clusters_duration"] = location_data.groupby(["local_segment"])["cluster_duration"].transform("sum")
    location_data["plogp"] = (location_data["cluster_duration"] / location_data["all_clusters_duration"]).apply(lambda x: x * np.log(x))
    
    entropy = -1 * location_data.groupby(["local_segment"])[["plogp"]].sum().rename(columns={"plogp": "locationentropy"})

    entropy["num_clusters"] = location_data.groupby(["local_segment"])["cluster_label"].nunique()
    entropy["normalizedlocationentropy"] = entropy["locationentropy"] / entropy["num_clusters"]

    return entropy



def doryab_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):

    location_data = pd.read_csv(sensor_data_files["sensor_data"])
    requested_features = provider["FEATURES"]
    dbscan_eps = provider["DBSCAN_EPS"]
    dbscan_minsamples = provider["DBSCAN_MINSAMPLES"]
    cluster_on = provider["CLUSTER_ON"]
    clustering_algorithm = provider["CLUSTERING_ALGORITHM"]
    radius_from_home = provider["RADIUS_FOR_HOME"]
    
    if provider["MINUTES_DATA_USED"]:
        requested_features.append("minutesdataused")

    # name of the features this function can compute
    base_features_names = ["locationvariance","loglocationvariance","totaldistance","avgspeed","varspeed","numberofsignificantplaces","numberlocationtransitions","radiusgyration","timeattop1location","timeattop2location","timeattop3location","movingtostaticratio","outlierstimepercent","maxlengthstayatclusters","minlengthstayatclusters","avglengthstayatclusters","stdlengthstayatclusters","locationentropy","normalizedlocationentropy","minutesdataused","timeathome","homelabel"]    
    # the subset of requested features this function can compute
    features_to_compute = list(set(requested_features) & set(base_features_names))
    
    location_data = apply_cluster_strategy(location_data, time_segment, clustering_algorithm, dbscan_eps, dbscan_minsamples, cluster_on, filter_data_by_segment)

    if location_data.empty:
        return pd.DataFrame(columns=["local_segment"] + features_to_compute)
    location_features = pd.DataFrame()

    # update distance after chunk_episodes() function
    location_data["distance"] = location_data["speed"] * (location_data["duration"] / 60) * 1000 # in meters

    location_features["minutesdataused"] = location_data[["local_segment", "duration"]].groupby(["local_segment"])["duration"].sum()

    # variance features
    location_features = variance_and_logvariance_features(location_data, location_features)

    # distance and speed features
    moving_data = location_data[location_data["is_stationary"] == 0].copy()
    location_features = location_features.merge(distance_and_speed_features(moving_data), how="outer", left_index=True, right_index=True)
    location_features[["totaldistance", "avgspeed", "varspeed"]] = location_features[["totaldistance", "avgspeed", "varspeed"]].fillna(0)

    # stationary features
    stationary_data = location_data[location_data["is_stationary"] == 1].copy()
    stationary_data_without_outliers = stationary_data[stationary_data["cluster_label"] != -1]

    location_features["numberofsignificantplaces"] = stationary_data_without_outliers.groupby(["local_segment"])["cluster_label"].nunique()
    # number of location transitions: ignores transitions from moving to static and vice-versa, but counts transitions from outliers to major location clusters
    location_features["numberlocationtransitions"] = stationary_data[["local_segment", "cluster_label"]].groupby(["local_segment"])["cluster_label"].apply(lambda x: np.sum(x != x.shift()) - 1)
    location_features["radiusgyration"] = radius_of_gyration(stationary_data_without_outliers)
    
    # stay at topn clusters features
    location_features = location_features.merge(stay_at_topn_clusters(stationary_data_without_outliers), how="outer", left_index=True, right_index=True)

    # moving to static ratio
    static_time = stationary_data.groupby(["local_segment"])["duration"].sum()
    total_time = location_data.groupby(["local_segment"])["duration"].sum()
    location_features["movingtostaticratio"] = static_time / total_time

    # outliers time percent
    outliers_time = stationary_data[stationary_data["cluster_label"] == -1].groupby(["local_segment"])["duration"].sum()
    location_features["outlierstimepercent"] = (outliers_time / static_time).fillna(0)

    # entropy features
    location_features = location_features.merge(location_entropy(stationary_data_without_outliers), how="outer", left_index=True, right_index=True)

    # time at home
    if stationary_data.empty:
        location_features["timeathome"] = 0
    else:
        stationary_data["time_at_home"] = stationary_data.apply(lambda row: row["duration"] if row["distance_from_home"] <= radius_from_home else 0, axis=1)
        location_features["timeathome"] = stationary_data[["local_segment", "time_at_home"]].groupby(["local_segment"])["time_at_home"].sum()

    # home label
    location_features["homelabel"] = stationary_data[["local_segment", "home_label"]].groupby(["local_segment"]).agg(lambda x: pd.Series.mode(x)[0])

    location_features = location_features[features_to_compute].reset_index()

    return location_features
