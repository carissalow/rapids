import warnings
import numpy as np
import pandas as pd
from doryab_clustering import haversine, create_clustering_hyperparameters, cluster



# Add "is_stationary" column to denote whether it is stationary or not
# "distance" and "speed" columns are also added
def mark_as_stationary(location_data, threshold_static):

    # Distance in meters
    location_data = location_data.assign(distance=haversine(location_data["double_longitude"], location_data["double_latitude"], location_data["double_longitude"].shift(-1), location_data["double_latitude"].shift(-1)))
    # Speed in km/h
    location_data.loc[:, "speed"] = (location_data["distance"] / location_data["duration_in_seconds"]).replace(np.inf, np.nan) * 3.6

    location_data.loc[:, "is_stationary"] = np.where(location_data["speed"] < threshold_static, 1, 0)

    location_data.dropna(subset=["duration_in_seconds"], inplace=True)
    return location_data

def infer_home_location(location_data, clustering_algorithm, hyperparameters, strategy, days_threshold):
    
    # Home locations are inferred based on records logged during midnight to 6am.
    # The home location is the mean coordinate of the home cluster. 
    if (strategy == "DORYAB_STRATEGY") or (strategy == "SUN_LI_VEGA_STRATEGY"):
        
        location_data_filtered = location_data[location_data["local_hour"] < 6]
    
        if location_data_filtered.empty:
            warnings.warn("We could not infer a home location because there are no location records logged during midnight to 6am.")
            return pd.DataFrame(columns=location_data_filtered.columns.tolist() + ["distance_from_home", "home_label"])
        
        location_data_filtered = cluster(location_data_filtered, clustering_algorithm, **hyperparameters)

        if strategy == "DORYAB_STRATEGY":

            # We assume the participant does not change the home location during the whole study.
            # The most common cluster of all nights are regarded as the home cluster.
            home_location = location_data_filtered[location_data_filtered["cluster_label"] == 1][["double_latitude", "double_longitude"]].mean()
            location_data["distance_from_home"] = haversine(location_data["double_longitude"], location_data["double_latitude"], [home_location["double_longitude"]] * location_data.shape[0], [home_location["double_latitude"]] * location_data.shape[0])
            location_data["home_label"] = 1

        else: # SUN_LI_VEGA_STRATEGY
            
            """
            We assume the participant might change the home location during the whole study.

            Each night will be assigned a candidate home location based on the following rules:
            if there are records within [03:30:00, 04:30:00]: (group 1)
                we choose the most common cluster during that period as the candidate of home cluster.
            elif there are records within [midnight, 03:30:00): (group 2)
                we choose the last valid cluster during that period as the candidate of home cluster.
            elif there are records within (04:30:00, 06:00:00]: (group 3)
                we choose the first valid cluster during that period as the candidate of home cluster.
            else:
                the home location is NA (missing) for that night.

            If the count of consecutive days with the same candidate home location cluster label is larger or equal to MINIMUM_DAYS_TO_DETECT_HOME_CHANGES,
            the candidate will be regarded as the home cluster; 
            otherwise, the home cluster will be the last valid day's cluster.
            (If there are no valid clusters before that day, it will be assigned the next valid day's cluster.)

            """

            # Split location data into 3 groups: [midnight, 03:30:00), [03:30:00, 04:30:00], (04:30:00, 06:00:00]
            location_data_filtered = location_data_filtered[~location_data_filtered["cluster_label"].isin([-1, np.nan])]
            location_data_filtered["group"] = location_data_filtered["local_time"].apply(lambda x: 1 if x >= "03:30:00" and x <= "04:30:00" else (2 if x < "03:30:00" else 3))
            
            # Select the smallest group number per day
            selected_groups = location_data_filtered[location_data_filtered["group"] == location_data_filtered.groupby("local_date")["group"].transform("min")][["group", "local_date", "cluster_label"]]
            
            # For group 1: [03:30:00, 04:30:00]
            group_1 = selected_groups[selected_groups["group"] == 1]
            home_clusters_group_1 = group_1.groupby(["local_date"]).agg(lambda x: pd.Series.mode(x)[0])
            # For group 2: [midnight, 03:30:00)
            group_2 = selected_groups[selected_groups["group"] == 2]
            home_clusters_group_2 = group_2.groupby(["local_date"]).last()
            # For group 3: (04:30:00, 06:00:00]
            group_3 = selected_groups[selected_groups["group"] == 3]
            home_clusters_group_3 = group_3.groupby(["local_date"]).first()
        
            home_clusters = pd.concat([home_clusters_group_1, home_clusters_group_2, home_clusters_group_3]).sort_index()
            
            # Count the consecutive days with the same candidate home location cluster label
            home_clusters["number_of_days"] = home_clusters.groupby((home_clusters["cluster_label"] != home_clusters["cluster_label"].shift(1)).cumsum())["cluster_label"].transform("count")
            # Assign the missing days with (1) the last valid day's cluster first and (2) the next valid day's cluster then
            home_clusters.loc[home_clusters["number_of_days"] < days_threshold, "cluster_label"] = np.nan
            location_data = location_data.merge(home_clusters[["cluster_label"]], left_on="local_date", right_index=True, how="left")
            location_data["cluster_label"] = location_data["cluster_label"].fillna(method="ffill").fillna(method="bfill")

            center_per_cluster = location_data_filtered.groupby(["cluster_label"])[["double_latitude", "double_longitude"]].mean().rename(columns={"double_latitude": "home_latitude", "double_longitude": "home_longitude"})
            location_data = location_data.merge(center_per_cluster, left_on="cluster_label", right_index=True, how="left")
            location_data["distance_from_home"] = haversine(location_data["double_longitude"], location_data["double_latitude"], location_data["home_longitude"], location_data["home_latitude"])    

            # reorder cluster labels
            reorder_mapping = {old_label: idx + 1 for idx, old_label in enumerate(location_data["cluster_label"].unique())}
            location_data["home_label"] = location_data["cluster_label"].map(reorder_mapping)

            location_data.drop(["cluster_label", "home_longitude", "home_latitude"], axis=1, inplace=True)

    return location_data



location_data = pd.read_csv(snakemake.input["sensor_input"])
provider = snakemake.params["provider"]

maximum_row_gap = provider["MAXIMUM_ROW_GAP"]
dbscan_eps = provider["DBSCAN_EPS"]
dbscan_minsamples = provider["DBSCAN_MINSAMPLES"]
threshold_static = provider["THRESHOLD_STATIC"]
clustering_algorithm = provider["CLUSTERING_ALGORITHM"]
cluster_on = provider["CLUSTER_ON"]
strategy = provider["INFER_HOME_LOCATION_STRATEGY"]
days_threshold = provider["MINIMUM_DAYS_TO_DETECT_HOME_CHANGES"]

if not location_data.timestamp.is_monotonic:
    location_data.sort_values(by=["timestamp"], inplace=True)

location_data["duration_in_seconds"] = -1 * location_data.timestamp.diff(-1) / 1000
location_data.loc[location_data["duration_in_seconds"] >= maximum_row_gap, "duration_in_seconds"] = np.nan

location_data = mark_as_stationary(location_data, threshold_static)

hyperparameters = create_clustering_hyperparameters(clustering_algorithm, dbscan_eps, dbscan_minsamples)
location_data_with_doryab_columns = infer_home_location(location_data, clustering_algorithm, hyperparameters, strategy, days_threshold)

selected_columns = ["local_timezone", "device_id", "start_timestamp", "end_timestamp", "provider", "double_latitude", "double_longitude", "distance", "speed", "is_stationary", "distance_from_home", "home_label"]
if cluster_on == "PARTICIPANT_DATASET":
    location_data_with_doryab_columns = cluster(location_data_with_doryab_columns, clustering_algorithm, **hyperparameters)
    selected_columns.append("cluster_label")

# Prepare for episodes
location_data_with_doryab_columns = location_data_with_doryab_columns.rename(columns={"timestamp": "start_timestamp"})
location_data_with_doryab_columns["end_timestamp"] = (location_data_with_doryab_columns["start_timestamp"] + location_data_with_doryab_columns["duration_in_seconds"] * 1000 - 1).astype(int)
location_data_with_doryab_columns[selected_columns].to_csv(snakemake.output[0], index=False)
