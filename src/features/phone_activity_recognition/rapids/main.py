import pandas as pd

def rapids_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):

    ar_episodes = pd.read_csv(sensor_data_files["sensor_episodes"])
    activity_classes = provider["ACTIVITY_CLASSES"]

    # name of the features this function can compute
    base_features_names = ["count","mostcommonactivity","countuniqueactivities","durationstationary","durationmobile","durationvehicle"]
    # the subset of requested features this function can compute
    requested_features = provider["FEATURES"]
    features_to_compute = list(set(requested_features) & set(base_features_names))

    ar_features = pd.DataFrame(columns=["local_segment"] + features_to_compute)
    if not ar_episodes.empty:
        ar_episodes = filter_data_by_segment(ar_episodes, time_segment)

        if not ar_episodes.empty:
            ar_features = pd.DataFrame()

            if "count" in features_to_compute:
                ar_features["count"] = ar_episodes.groupby(["local_segment"]).count()["episode_id"]
            if "mostcommonactivity" in features_to_compute:
                ar_features["mostcommonactivity"] = ar_episodes.groupby(["local_segment"])["activity_type"].agg(lambda x: pd.Series.mode(x)[0])
            if "countuniqueactivities" in features_to_compute:
                ar_features["countuniqueactivities"] = ar_episodes.groupby(["local_segment"])["activity_type"].nunique()
            
            # duration features    
            for column, activity_labels in activity_classes.items():
                if "duration" + column.lower() in features_to_compute:
                    filtered_data = ar_episodes[ar_episodes["activity_name"].isin(pd.Series(activity_labels))]
                    if not filtered_data.empty:
                        ar_features["duration" + column.lower()] = ar_episodes[ar_episodes["activity_name"].isin(pd.Series(activity_labels))].groupby(["local_segment"])["duration"].sum()
                    else:
                        ar_features["duration" + column.lower()] = 0

            ar_features.index.names = ["local_segment"]
            ar_features = ar_features.reset_index()
    
    ar_features.fillna(value={"count": 0, "countuniqueactivities": 0, "durationstationary": 0, "durationmobile": 0, "durationvehicle": 0}, inplace=True)

    return ar_features
