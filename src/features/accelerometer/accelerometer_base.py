import pandas as pd
import numpy as np

def getActivityEpisodes(acc_minute):    
    # rebuild local date time for resampling
    acc_minute["local_datetime"] = pd.to_datetime(acc_minute["local_date"].dt.strftime("%Y-%m-%d") + \
                                        " " + acc_minute["local_hour"].apply(str) + ":" + acc_minute["local_minute"].apply(str) + ":00")

    # resample the data into 1 minute bins, set "isexertionalactivity" column to be NA if it is missing
    resampled_acc_minute = pd.DataFrame(acc_minute.resample("1T", on="local_datetime")["isexertionalactivity"].sum(min_count=1))

    # group rows by consecutive values of "isexertionalactivity" column
    group = pd.DataFrame(resampled_acc_minute["isexertionalactivity"] != resampled_acc_minute["isexertionalactivity"].shift()).cumsum().rename(columns={"isexertionalactivity": "group_idx"})

    # combine resampled_acc_minute and group column
    resampled_acc_minute = pd.concat([resampled_acc_minute, group], axis=1)

    # drop rows where "isexertionalactivity" column is missing and reset the index
    resampled_acc_minute.dropna(subset=["isexertionalactivity"], inplace=True)
    resampled_acc_minute.reset_index(inplace=True)
    resampled_acc_minute.loc[:, "local_date"] = resampled_acc_minute["local_datetime"].dt.date

    # duration column contains the number of minutes (rows) of exertional and nonexertional activity for each episode
    activity_episode = resampled_acc_minute.groupby(["isexertionalactivity", "group_idx", "local_date"]).count().rename(columns={"local_datetime": "duration"}).reset_index()

    return activity_episode

def dropRowsWithCertainThreshold(data, threshold):
    data_grouped = data.groupby(["local_date", "local_hour", "local_minute"]).count()
    drop_dates = data_grouped[data_grouped["timestamp"] == threshold].index
    data.set_index(["local_date", "local_hour", "local_minute"], inplace = True)
    if not drop_dates.empty:
        data.drop(drop_dates, axis = 0, inplace = True)
    return data.reset_index()

def statsFeatures(acc_data, day_segment, features_to_compute, features_type, acc_features):
    if features_type == "magnitude":
        col_name = features_type
    elif features_type == "durationexertionalactivityepisode" or features_type == "durationnonexertionalactivityepisode":
        col_name = "duration"
    else:
        raise ValueError("features_type can only be one of ['magnitude', 'durationexertionalactivityepisode', 'durationnonexertionalactivityepisode'].")

    if "sum" + features_type in features_to_compute:
        acc_features["acc_" + day_segment + "_sum" + features_type] = acc_data.groupby(["local_date"])[col_name].sum()
    if "max" + features_type in features_to_compute:
        acc_features["acc_" + day_segment + "_max" + features_type] = acc_data.groupby(["local_date"])[col_name].max()
    if "min" + features_type in features_to_compute:
        acc_features["acc_" + day_segment + "_min" + features_type] = acc_data.groupby(["local_date"])[col_name].min()
    if "avg" + features_type in features_to_compute:
        acc_features["acc_" + day_segment + "_avg" + features_type] = acc_data.groupby(["local_date"])[col_name].mean()
    if "median" + features_type in features_to_compute:
        acc_features["acc_" + day_segment + "_median" + features_type] = acc_data.groupby(["local_date"])[col_name].median()
    if "std" + features_type in features_to_compute:
        acc_features["acc_" + day_segment + "_std" + features_type] = acc_data.groupby(["local_date"])[col_name].std()

    return acc_features



def base_accelerometer_features(acc_data, day_segment, requested_features, valid_sensed_minutes):
    # name of the features this function can compute
    base_features_names_magnitude = ["maxmagnitude", "minmagnitude", "avgmagnitude", "medianmagnitude", "stdmagnitude"]
    base_features_names_exertionalactivityepisode = ["sumdurationexertionalactivityepisode", "maxdurationexertionalactivityepisode", "mindurationexertionalactivityepisode", "avgdurationexertionalactivityepisode", "mediandurationexertionalactivityepisode", "stddurationexertionalactivityepisode"]
    base_features_names_nonexertionalactivityepisode = ["sumdurationnonexertionalactivityepisode", "maxdurationnonexertionalactivityepisode", "mindurationnonexertionalactivityepisode", "avgdurationnonexertionalactivityepisode", "mediandurationnonexertionalactivityepisode", "stddurationnonexertionalactivityepisode"]
    # the subset of requested features this function can compute
    features_to_compute_magnitude = list(set(requested_features["magnitude"]) & set(base_features_names_magnitude))
    features_to_compute_exertionalactivityepisode = list(set(requested_features["exertional_activity_episode"]) & set(base_features_names_exertionalactivityepisode))
    features_to_compute_nonexertionalactivityepisode = list(set(requested_features["nonexertional_activity_episode"]) & set(base_features_names_nonexertionalactivityepisode))

    features_to_compute = features_to_compute_magnitude + features_to_compute_exertionalactivityepisode + features_to_compute_nonexertionalactivityepisode + (["validsensedminutes"] if valid_sensed_minutes else [])

    acc_features = pd.DataFrame(columns=["local_date"] + ["acc_" + day_segment + "_" + x for x in features_to_compute])
    if not acc_data.empty:
        if day_segment != "daily":
            acc_data = acc_data[acc_data["local_day_segment"] == day_segment]

        if not acc_data.empty:
            acc_features = pd.DataFrame()        
            # get magnitude related features: magnitude = sqrt(x^2+y^2+z^2)
            magnitude = acc_data.apply(lambda row: np.sqrt(row["double_values_0"] ** 2 + row["double_values_1"] ** 2 + row["double_values_2"] ** 2), axis=1)
            acc_data = acc_data.assign(magnitude = magnitude.values)
            acc_features = statsFeatures(acc_data, day_segment, features_to_compute_magnitude, "magnitude", acc_features)

            
            # get extertional activity features
            # reference: https://jamanetwork.com/journals/jamasurgery/fullarticle/2753807

            # drop rows where we only have one row per minute (no variance) 
            acc_data = dropRowsWithCertainThreshold(acc_data, 1)
            if not acc_data.empty:
                # check if the participant performs exertional activity for each minute
                acc_minute = pd.DataFrame()
                acc_minute["isexertionalactivity"] = (acc_data.groupby(["local_date", "local_hour", "local_minute"])["double_values_0"].var() + acc_data.groupby(["local_date", "local_hour", "local_minute"])["double_values_1"].var() + acc_data.groupby(["local_date", "local_hour", "local_minute"])["double_values_2"].var()).apply(lambda x: 1 if x > 0.15 * (9.807 ** 2) else 0)
                acc_minute.reset_index(inplace=True)

                if valid_sensed_minutes:
                    acc_features["acc_" + day_segment + "_validsensedminutes"] = acc_minute.groupby(["local_date"])["isexertionalactivity"].count()
                
                activity_episode = getActivityEpisodes(acc_minute)
                exertionalactivity_episodes = activity_episode[activity_episode["isexertionalactivity"] == 1]
                acc_features = statsFeatures(exertionalactivity_episodes, day_segment, features_to_compute_exertionalactivityepisode, "durationexertionalactivityepisode", acc_features)

                nonexertionalactivity_episodes = activity_episode[activity_episode["isexertionalactivity"] == 0]
                acc_features = statsFeatures(nonexertionalactivity_episodes, day_segment, features_to_compute_nonexertionalactivityepisode, "durationnonexertionalactivityepisode", acc_features)

                acc_features[[colname for colname in acc_features.columns if "std" not in colname]] = acc_features[[colname for colname in acc_features.columns if "std" not in colname]].fillna(0)
            
            acc_features = acc_features.reset_index()

    return acc_features
