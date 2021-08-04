import pandas as pd
import numpy as np

def dropRowsWithCertainThreshold(data, threshold):
    data_grouped = data.groupby(["local_timezone", "local_segment", "local_date", "local_hour", "local_minute"])
    data_cleaned = data_grouped.filter(lambda x: x["timestamp"].count() > threshold)
    return data_cleaned

def getActivityEpisodes(acc_minute):    
    # rebuild local date time for resampling
    acc_minute["local_datetime"] = pd.to_datetime(acc_minute["local_date"] + \
                                        " " + acc_minute["local_hour"].apply(str) + ":" + acc_minute["local_minute"].apply(str) + ":00")

    # compute time interval between consecutive rows in minutes
    acc_minute["rows_interval"] = round(acc_minute["local_datetime"].diff().dt.total_seconds() / 60, 0)

    # put consecutive rows into the same group if (1) the interval between two rows is 1 minute and (2) have the same values of "isexertionalactivity", "local_timezone", and "local_segment"
    acc_minute["group_idx"] =  ((acc_minute[["isexertionalactivity", "local_timezone", "local_segment"]].shift() != acc_minute[["isexertionalactivity", "local_timezone", "local_segment"]]).any(axis=1) | (acc_minute["rows_interval"] != 1)).cumsum()
    
    # get activity episodes: duration column contains the number of minutes (rows) of exertional and nonexertional activity for each episode
    grouped = acc_minute.groupby("group_idx")
    activity_episodes = grouped["local_segment"].agg(duration="count")
    activity_episodes[["local_segment", "isexertionalactivity"]] = grouped[["local_segment", "isexertionalactivity"]].first()
    
    return activity_episodes

def statsFeatures(acc_data, features_to_compute, features_type, acc_features):
    if "sum" + features_type in features_to_compute:
        acc_features["sum" + features_type] = acc_data.groupby(["local_segment"])["duration"].sum()
    if "max" + features_type in features_to_compute:
        acc_features["max" + features_type] = acc_data.groupby(["local_segment"])["duration"].max()
    if "min" + features_type in features_to_compute:
        acc_features["min" + features_type] = acc_data.groupby(["local_segment"])["duration"].min()
    if "avg" + features_type in features_to_compute:
        acc_features["avg" + features_type] = acc_data.groupby(["local_segment"])["duration"].mean()
    if "median" + features_type in features_to_compute:
        acc_features["median" + features_type] = acc_data.groupby(["local_segment"])["duration"].median()
    if "std" + features_type in features_to_compute:
        acc_features["std" + features_type] = acc_data.groupby(["local_segment"])["duration"].std()

    return acc_features



def panda_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):

    acc_data = pd.read_csv(sensor_data_files["sensor_data"])
    requested_features = provider["FEATURES"]
    valid_sensed_minutes = provider["VALID_SENSED_MINUTES"]
    # name of the features this function can compute
    base_features_names_exertionalactivityepisode = ["sumdurationexertionalactivityepisode", "maxdurationexertionalactivityepisode", "mindurationexertionalactivityepisode", "avgdurationexertionalactivityepisode", "mediandurationexertionalactivityepisode", "stddurationexertionalactivityepisode"]
    base_features_names_nonexertionalactivityepisode = ["sumdurationnonexertionalactivityepisode", "maxdurationnonexertionalactivityepisode", "mindurationnonexertionalactivityepisode", "avgdurationnonexertionalactivityepisode", "mediandurationnonexertionalactivityepisode", "stddurationnonexertionalactivityepisode"]
    # the subset of requested features this function can compute
    features_to_compute_exertionalactivityepisode = list(set([x + "exertionalactivityepisode" for x in requested_features["exertional_activity_episode"]]) & set(base_features_names_exertionalactivityepisode))
    features_to_compute_nonexertionalactivityepisode = list(set([ x + "nonexertionalactivityepisode" for x in requested_features["nonexertional_activity_episode"]]) & set(base_features_names_nonexertionalactivityepisode))

    features_to_compute =  features_to_compute_exertionalactivityepisode + features_to_compute_nonexertionalactivityepisode + (["validsensedminutes"] if valid_sensed_minutes else [])

    acc_features = pd.DataFrame(columns=["local_segment"] + features_to_compute)
    if not acc_data.empty:
        acc_data = filter_data_by_segment(acc_data, time_segment)
        
        
        if not acc_data.empty:
            # drop rows where we only have one row per minute (no variance) 
            acc_data = dropRowsWithCertainThreshold(acc_data, 1)

            if not acc_data.empty:
                acc_features = pd.DataFrame()
                # check if the participant performs exertional activity for each minute
                acc_minute = pd.DataFrame()
                acc_minute["isexertionalactivity"] = (acc_data.groupby(["local_timezone", "local_segment", "local_date", "local_hour", "local_minute"])["double_values_0"].var() + acc_data.groupby(["local_timezone", "local_segment", "local_date", "local_hour", "local_minute"])["double_values_1"].var() + acc_data.groupby(["local_timezone", "local_segment", "local_date", "local_hour", "local_minute"])["double_values_2"].var()).apply(lambda x: 1 if x > 0.15 * (9.807 ** 2) else 0)
                acc_minute.reset_index(inplace=True)

                if valid_sensed_minutes:
                    acc_features["validsensedminutes"] = acc_minute.groupby(["local_segment"])["isexertionalactivity"].count()
                
                activity_episodes = getActivityEpisodes(acc_minute)
                # compute exertional episodes features
                exertionalactivity_episodes = activity_episodes[activity_episodes["isexertionalactivity"] == 1]
                acc_features = statsFeatures(exertionalactivity_episodes, features_to_compute_exertionalactivityepisode, "durationexertionalactivityepisode", acc_features)
                # compute non-exertional episodes features
                nonexertionalactivity_episodes = activity_episodes[activity_episodes["isexertionalactivity"] == 0]
                acc_features = statsFeatures(nonexertionalactivity_episodes, features_to_compute_nonexertionalactivityepisode, "durationnonexertionalactivityepisode", acc_features)

                acc_features[[colname for colname in acc_features.columns if "std" not in colname]] = acc_features[[colname for colname in acc_features.columns if "std" not in colname]].fillna(0)

            acc_features = acc_features.reset_index()

    return acc_features
