import pandas as pd
import numpy as np

def rapids_features(sensor_data_files, day_segment, provider, filter_data_by_segment, *args, **kwargs):

    chunk_episodes = kwargs["chunk_episodes"]

    ar_episodes = pd.read_csv(sensor_data_files["sensor_episodes"])
    activity_classes = provider["ACTIVITY_CLASSES"]

    # name of the features this function can compute
    base_features_names = ["count","mostcommonactivity","countuniqueactivities","durationstationary","durationmobile","durationvehicle"]
    # the subset of requested features this function can compute
    requested_features = provider["FEATURES"]
    features_to_compute = list(set(requested_features) & set(base_features_names))

    ar_features = pd.DataFrame(columns=["local_segment"] + ["ar_rapids_" + x for x in features_to_compute])
    if not ar_episodes.empty:
        ar_episodes = filter_data_by_segment(ar_episodes, day_segment)

        if not ar_episodes.empty:
            # chunk episodes
            ar_episodes = chunk_episodes(ar_episodes)

        if not ar_episodes.empty:
            ar_features = pd.DataFrame()

            if "count" in features_to_compute:
                ar_features["ar_rapids_count"] = ar_episodes.groupby(["local_segment"]).count()["episode_id"]
            if "mostcommonactivity" in features_to_compute:
                ar_features["ar_rapids_mostcommonactivity"] = ar_episodes.groupby(["local_segment"])["activity_type"].agg(lambda x: pd.Series.mode(x)[0])
            if "countuniqueactivities" in features_to_compute:
                ar_features["ar_rapids_countuniqueactivities"] = ar_episodes.groupby(["local_segment"])["activity_type"].nunique()
            
            # duration features    
            for column, activity_labels in activity_classes.items():
                if "duration" + column.lower() in features_to_compute:
                    filtered_data = ar_episodes[ar_episodes["activity_name"].isin(pd.Series(activity_labels))]
                    if not filtered_data.empty:
                        ar_features["ar_rapids_duration_" + column] = ar_episodes[ar_episodes["activity_name"].isin(pd.Series(activity_labels))].groupby(["local_segment"])["duration"].sum().fillna(0)
                    else:
                        ar_features["ar_rapids_duration_" + column] = 0

            ar_features.index.names = ["local_segment"]
            ar_features = ar_features.reset_index()
    
    return ar_features






























    """

    if not ar_data.empty:
        ar_data = filter_data_by_segment(ar_data, day_segment)

        if not ar_data.empty:
            # chunk_episodes
            ar_data = chunk_episodes(ar_data)

        if not ar_data.empty:
        
            ar_data["episode_id"] = ((ar_data.ar_status != ar_data.ar_status.shift()) | (ar_data.start_timestamp - ar_data.end_timestamp.shift() > 1)).cumsum()
            grouped = ar_data.groupby(by=["local_segment", "episode_id", "ar_status"])
            ar_episodes= grouped[["duration"]].sum()
            ar_episodes["ar_diff"] = grouped["ar_level"].first() - grouped["ar_level"].last()
            ar_episodes["ar_consumption_rate"] = ar_episodes["ar_diff"] / ar_episodes["duration"]
            ar_episodes.reset_index(inplace=True)

            # for discharge episodes
            ar_discharge_episodes = ar_episodes[(ar_episodes["ar_status"] == 3) | (ar_episodes["ar_status"] == 4)]
            ar_discharge_features = pd.DataFrame()
            if "countdischarge" in features_to_compute:
                ar_discharge_features["ar_rapids_countdischarge"] = ar_discharge_episodes.groupby(["local_segment"])["episode_id"].count()
            if "sumdurationdischarge" in features_to_compute:
                ar_discharge_features["ar_rapids_sumdurationdischarge"] = ar_discharge_episodes.groupby(["local_segment"])["duration"].sum()
            if "avgconsumptionrate" in features_to_compute:
                ar_discharge_features["ar_rapids_avgconsumptionrate"] = ar_discharge_episodes.groupby(["local_segment"])["ar_consumption_rate"].mean()
            if "maxconsumptionrate" in features_to_compute:
                ar_discharge_features["ar_rapids_maxconsumptionrate"] = ar_discharge_episodes.groupby(["local_segment"])["ar_consumption_rate"].max()

            # for charge episodes
            ar_charge_episodes = ar_episodes[(ar_episodes["ar_status"] == 2) | (ar_episodes["ar_status"] == 5)]
            ar_charge_features = pd.DataFrame()
            if "countcharge" in features_to_compute:
                ar_charge_features["ar_rapids_countcharge"] = ar_charge_episodes.groupby(["local_segment"])["episode_id"].count()
            if "sumdurationcharge" in features_to_compute:
                ar_charge_features["ar_rapids_sumdurationcharge"] = ar_charge_episodes.groupby(["local_segment"])["duration"].sum()
            
            # combine discharge features and charge features; fill the missing values with ZERO
            ar_features = pd.concat([ar_discharge_features, ar_charge_features], axis=1, sort=True).fillna(0)

            ar_features.index.rename("local_segment", inplace=True)
            ar_features = ar_features.reset_index()

    return ar_features
    """