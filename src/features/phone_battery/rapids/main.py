import pandas as pd
from datetime import datetime, timedelta, time

def rapids_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):
    
    battery_data = pd.read_csv(sensor_data_files["sensor_episodes"])

    # name of the features this function can compute
    base_features_names = ["countdischarge", "sumdurationdischarge", "countcharge", "sumdurationcharge", "avgconsumptionrate", "maxconsumptionrate"]
    # the subset of requested features this function can compute
    requested_features = provider["FEATURES"]
    features_to_compute = list(set(requested_features) & set(base_features_names))

    battery_features = pd.DataFrame(columns=["local_segment"] + features_to_compute)
    if not battery_data.empty:
        battery_data = filter_data_by_segment(battery_data, time_segment)

        if not battery_data.empty:
        
            battery_data["episode_id"] = ((battery_data.battery_status != battery_data.battery_status.shift()) | (battery_data.start_timestamp - battery_data.end_timestamp.shift() > 1)).cumsum()
            grouped = battery_data.groupby(by=["local_segment", "episode_id", "battery_status"])
            battery_episodes= grouped[["duration"]].sum()
            battery_episodes["battery_diff"] = grouped["battery_level"].first() - grouped["battery_level"].last()
            battery_episodes["battery_consumption_rate"] = battery_episodes["battery_diff"] / battery_episodes["duration"]
            battery_episodes.reset_index(inplace=True)

            # for discharge episodes
            battery_discharge_episodes = battery_episodes[(battery_episodes["battery_status"] == 3) | (battery_episodes["battery_status"] == 4)]
            battery_discharge_episodes = battery_discharge_episodes[battery_discharge_episodes['battery_consumption_rate'] !=0 ]
            battery_discharge_features = pd.DataFrame()
            if "countdischarge" in features_to_compute:
                battery_discharge_features["countdischarge"] = battery_discharge_episodes.groupby(["local_segment"])["episode_id"].count()
            if "sumdurationdischarge" in features_to_compute:
                battery_discharge_features["sumdurationdischarge"] = battery_discharge_episodes.groupby(["local_segment"])["duration"].sum()
            if "avgconsumptionrate" in features_to_compute:
                battery_discharge_features["avgconsumptionrate"] = battery_discharge_episodes.groupby(["local_segment"])["battery_consumption_rate"].mean()
            if "maxconsumptionrate" in features_to_compute:
                battery_discharge_features["maxconsumptionrate"] = battery_discharge_episodes.groupby(["local_segment"])["battery_consumption_rate"].max()

            # for charge episodes
            battery_charge_episodes = battery_episodes[(battery_episodes["battery_status"] == 2) | (battery_episodes["battery_status"] == 5)]
            battery_charge_episodes = battery_charge_episodes[battery_charge_episodes['battery_consumption_rate'] !=0 ]
            battery_charge_features = pd.DataFrame()
            if "countcharge" in features_to_compute:
                battery_charge_features["countcharge"] = battery_charge_episodes.groupby(["local_segment"])["episode_id"].count()
            if "sumdurationcharge" in features_to_compute:
                battery_charge_features["sumdurationcharge"] = battery_charge_episodes.groupby(["local_segment"])["duration"].sum()
            
            # combine discharge features and charge features; fill the missing values with ZERO
            battery_features = pd.concat([battery_discharge_features, battery_charge_features], axis=1, sort=True).fillna(0)

            battery_features.index.rename("local_segment", inplace=True)
            battery_features = battery_features.reset_index()

    return battery_features