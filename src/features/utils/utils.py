rapids_log_tag =  "RAPIDS:"

def filter_data_by_segment(data, day_segment):
    date_regex = "[0-9]{4}[\-|\/][0-9]{2}[\-|\/][0-9]{2}"
    hour_regex = "[0-9]{2}:[0-9]{2}:[0-9]{2}"
    segment_regex = "\[({}#{}#{}#{}#{})\]".format(day_segment, date_regex, hour_regex, date_regex, hour_regex)
    data["local_segment"] = data["assigned_segments"].str.extract(segment_regex, expand=True)
    return(data.dropna(subset = ["local_segment"]))

def chunk_episodes(sensor_episodes):
    import pytz, copy
    import pandas as pd
    from datetime import datetime

    # avoid warning messages: SettingWithCopyWarning
    sensor_episodes = sensor_episodes.copy()

    # convert string to datetime with local timezone
    sensor_episodes["start_datetime"] = pd.to_datetime(sensor_episodes["local_segment"].str[-39:-20], format="%Y-%m-%d#%H:%M:%S")
    sensor_episodes["start_datetime"] = pd.concat([data["start_datetime"].dt.tz_localize(tz) for tz, data in sensor_episodes.groupby("local_timezone")])

    sensor_episodes["end_datetime"] = pd.to_datetime(sensor_episodes["local_segment"].str[-19:], format="%Y-%m-%d#%H:%M:%S")
    sensor_episodes["end_datetime"] = pd.concat([data["end_datetime"].dt.tz_localize(tz) for tz, data in sensor_episodes.groupby("local_timezone")])

    # unix timestamp in milliseconds
    sensor_episodes["start_timestamp"] = sensor_episodes["start_datetime"].apply(lambda dt: dt.timestamp() * 1000)
    sensor_episodes["end_timestamp"] = sensor_episodes["end_datetime"].apply(lambda dt: dt.timestamp() * 1000)

    # compute chunked timestamp
    sensor_episodes["chunked_start_timestamp"] = sensor_episodes[["timestamp", "start_timestamp"]].max(axis=1)

    sensor_episodes["timestamp_plus_duration"] = sensor_episodes["timestamp"] + sensor_episodes["duration"] * 1000 * 60
    sensor_episodes["chunked_end_timestamp"] = sensor_episodes[["timestamp_plus_duration", "end_timestamp"]].min(axis=1)

    # time_diff: intersection of current row and segment
    sensor_episodes["time_diff"] = (sensor_episodes["chunked_end_timestamp"] - sensor_episodes["chunked_start_timestamp"]) / (1000 * 60)

    # compute chunked datetime
    sensor_episodes["chunked_start_datetime"] = pd.to_datetime(sensor_episodes["chunked_start_timestamp"], unit="ms", utc=True)
    sensor_episodes["chunked_start_datetime"] = pd.concat([data["chunked_start_datetime"].dt.tz_convert(tz) for tz, data in sensor_episodes.groupby("local_timezone")])

    sensor_episodes["chunked_end_datetime"] = pd.to_datetime(sensor_episodes["chunked_end_timestamp"], unit="ms", utc=True)
    sensor_episodes["chunked_end_datetime"] = pd.concat([data["chunked_end_datetime"].dt.tz_convert(tz) for tz, data in sensor_episodes.groupby("local_timezone")])

    # merge episodes
    sensor_episodes_grouped = sensor_episodes.groupby(["episode_id", "episode", "screen_sequence"])
    merged_sensor_episodes = sensor_episodes_grouped[["time_diff"]].sum()
    merged_sensor_episodes["local_segment"] = sensor_episodes_grouped["local_segment"].first()

    merged_sensor_episodes["start_timestamp"] = sensor_episodes_grouped["chunked_start_timestamp"].first()
    merged_sensor_episodes["end_timestamp"] = sensor_episodes_grouped["chunked_end_timestamp"].last()

    merged_sensor_episodes["local_start_date_time"] = sensor_episodes_grouped["chunked_start_datetime"].first()
    merged_sensor_episodes["local_end_date_time"] = sensor_episodes_grouped["chunked_end_datetime"].last()

    merged_sensor_episodes.reset_index(inplace=True)

    return merged_sensor_episodes

def fetch_provider_features(provider, provider_key, config_key, sensor_data_file, day_segments_file):
    import pandas as pd
    from importlib import import_module, util

    sensor_features = pd.DataFrame(columns=["local_segment"])
    sensor_data = pd.read_csv(sensor_data_file)
    day_segments_labels = pd.read_csv(day_segments_file, header=0)
    if "FEATURES" not in provider:
        raise ValueError("Provider config[{}][PROVIDERS][{}] is missing a FEATURES attribute in config.yaml".format(config_key.upper(), provider_key))

    if provider["COMPUTE"] == True:
            code_path = provider["SRC_FOLDER"] + ".main"
            feature_module = import_module(code_path)
            feature_function = getattr(feature_module,  provider["SRC_FOLDER"] + "_features")
            
            for day_segment in day_segments_labels["label"]:
                    print("{} Processing {} {} {}".format(rapids_log_tag, config_key, provider_key, day_segment))
                    features = feature_function(sensor_data, day_segment, provider, filter_data_by_segment=filter_data_by_segment, chunk_episodes=chunk_episodes)
                    sensor_features = sensor_features.merge(features, how="outer")
    else:
            for feature in provider["FEATURES"]:
                    sensor_features[feature] = None
    segment_colums = pd.DataFrame()
    split_segemnt_columns = sensor_features["local_segment"].str.split(pat="#", expand=True)
    new_segment_columns = split_segemnt_columns if split_segemnt_columns.shape[1] == 5 else pd.DataFrame(columns=["local_segment_label", "local_start_date", "local_start_time", "local_end_date", "local_end_time"])
    segment_colums[["local_segment_label", "local_start_date", "local_start_time", "local_end_date", "local_end_time"]] = new_segment_columns
    for i in range(segment_colums.shape[1]):
            sensor_features.insert(1 + i, segment_colums.columns[i], segment_colums[segment_colums.columns[i]])

    return sensor_features
