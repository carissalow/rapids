import numpy as np
import pandas as pd
from phone_locations.doryab.doryab_clustering import haversine, create_clustering_hyperparameters, cluster



def apply_cluster_strategy(location_data, time_segment, clustering_algorithm, dbscan_eps, dbscan_minsamples, cluster_on, filter_data_by_segment):

    hyperparameters = create_clustering_hyperparameters(clustering_algorithm, dbscan_eps, dbscan_minsamples)

    if cluster_on == "PARTICIPANT_DATASET":
        # clusters are created in cluster_accross_participant_dataset.py script
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

def distance_and_speed_features(moving_data):

    distance_and_speed = moving_data[["local_segment", "distance"]].groupby(["local_segment"]).sum().rename(columns={"distance": "totaldistance"})
    
    moving_data_grouped = moving_data[["local_segment", "speed"]].groupby(["local_segment"])
    distance_and_speed["avgspeed"] = moving_data_grouped["speed"].mean()
    distance_and_speed["varspeed"] = moving_data_grouped["speed"].var()
    
    return distance_and_speed

def radius_of_gyration(location_data):
 
    # center is the centroid of the places visited during a segment instance, not the home location
    clusters = location_data.groupby(["local_segment", "cluster_label"]).agg(
        double_latitude=("double_latitude", "mean"),
        double_longitude=("double_longitude", "mean"),
        time_in_a_cluster=("duration_in_seconds", "sum")
    ).reset_index()
   
    clusters[["centroid_double_latitude", "centroid_double_longitude"]] = clusters.groupby(["local_segment"], sort=False)[["double_latitude", "double_longitude"]].transform("mean")
    clusters["distance_squared"] = haversine(clusters["double_longitude"], clusters["double_latitude"], clusters["centroid_double_longitude"], clusters["centroid_double_latitude"]) ** 2
    
    clusters["distance_squared_X_time_in_a_cluster"] = clusters["distance_squared"] * clusters["time_in_a_cluster"]
    rog = np.sqrt(clusters.groupby(["local_segment"])["distance_squared_X_time_in_a_cluster"].sum() / clusters.groupby(["local_segment"])["time_in_a_cluster"].sum().replace(0, np.inf))
   
    return rog

def cluster_stay(x, stay_at_clusters, cluster_n):
    topn_cluster_label = x[stay_at_clusters.loc[x.index]["cluster_label"] == cluster_n]
    time_at_topn = topn_cluster_label.iloc[0] if len(topn_cluster_label) == 1 else None
    return time_at_topn

def stay_at_topn_clusters(location_data):

    stay_at_clusters = location_data[["local_segment", "cluster_label", "duration_in_seconds"]].groupby(["local_segment", "cluster_label"], sort=True).sum().reset_index()
    stay_at_clusters["duration_in_minutes"] = stay_at_clusters["duration_in_seconds"] / 60

    stay_at_clusters_features = stay_at_clusters.groupby(["local_segment"]).agg(        
        timeattop1location=("duration_in_minutes", lambda x: cluster_stay(x, stay_at_clusters, 1)),
        timeattop2location=("duration_in_minutes", lambda x: cluster_stay(x, stay_at_clusters, 2)),
        timeattop3location=("duration_in_minutes", lambda x: cluster_stay(x, stay_at_clusters, 3)),
        maxlengthstayatclusters=("duration_in_minutes", "max"),
        minlengthstayatclusters=("duration_in_minutes", "min"),
        avglengthstayatclusters=("duration_in_minutes", "mean"),
        stdlengthstayatclusters=("duration_in_minutes", "std")
    ).fillna(0)

    return stay_at_clusters_features

def location_entropy(location_data):

    location_data = location_data.groupby(["local_segment", "cluster_label"])[["duration_in_seconds"]].sum().reset_index().rename(columns={"duration_in_seconds": "cluster_duration"})
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

    location_features["minutesdataused"] = location_data.drop_duplicates(subset=["local_segment", "local_date", "local_hour", "local_minute"])[["local_segment", "local_minute"]].groupby(["local_segment"])["local_minute"].count()

    # variance features
    location_features["locationvariance"] = location_data.groupby(["local_segment"])["double_latitude"].var() + location_data.groupby(["local_segment"])["double_longitude"].var()
    location_features["loglocationvariance"] = np.log10(location_features["locationvariance"]).replace(-np.inf, np.nan)

    # distance and speed features
    moving_data = location_data[location_data["is_stationary"] == 0]
    location_features = location_features.merge(distance_and_speed_features(moving_data), how="outer", left_index=True, right_index=True)

    # stationary features
    stationary_data = location_data[location_data["is_stationary"] == 1]
    stationary_data_without_outliers = stationary_data[stationary_data["cluster_label"] != -1]

    location_features["numberofsignificantplaces"] = stationary_data_without_outliers.groupby(["local_segment"])["cluster_label"].nunique()
    # number of location transitions: ignores transitions from moving to static and vice-versa, but counts transitions from outliers to major location clusters
    location_features["numberlocationtransitions"] = stationary_data[["local_segment", "cluster_label"]].groupby(["local_segment"])["cluster_label"].apply(lambda x: np.sum(x != x.shift()) - 1)
    location_features["radiusgyration"] = radius_of_gyration(stationary_data_without_outliers)
    
    # stay at topn clusters features
    location_features = location_features.merge(stay_at_topn_clusters(stationary_data_without_outliers), how="outer", left_index=True, right_index=True)

    # moving to static ratio
    static_time = stationary_data.groupby(["local_segment"])["duration_in_seconds"].sum()
    total_time = location_data.groupby(["local_segment"])["duration_in_seconds"].sum()
    location_features["movingtostaticratio"] = static_time / total_time

    # outliers time percent
    outliers_time = stationary_data[stationary_data["cluster_label"] == -1].groupby(["local_segment"])["duration_in_seconds"].sum()
    location_features["outlierstimepercent"] = outliers_time / static_time

    # entropy features
    location_features = location_features.merge(location_entropy(stationary_data_without_outliers), how="outer", left_index=True, right_index=True)

    # time at home
    location_features["timeathome"] = stationary_data[stationary_data["distance_from_home"] <= radius_from_home].groupby(["local_segment"])["duration_in_seconds"].sum() / 60

    # home label
    location_features["homelabel"] = stationary_data[["local_segment", "home_label"]].groupby(["local_segment"]).agg(lambda x: pd.Series.mode(x)[0])

    location_features = location_features[features_to_compute].reset_index()

    return location_features
