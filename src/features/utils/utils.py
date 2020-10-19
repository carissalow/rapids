rapids_log_tag =  "RAPIDS:"

def filter_data_by_segment(data, day_segment):
    datetime_regex = "[0-9]{4}[\-|\/][0-9]{2}[\-|\/][0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}"
    timestamps_regex = "[0-9]{13}"
    segment_regex = "\[({}#{},{};{},{})\]".format(day_segment, datetime_regex, datetime_regex, timestamps_regex, timestamps_regex)
    data["local_segment"] = data["assigned_segments"].str.extract(segment_regex, expand=True)
    data = data.drop(columns=["assigned_segments"])
    data = data.dropna(subset = ["local_segment"])
    if(data.shape[0] == 0): # there are no rows belonging to day_segment
        data["timestamps_segment"] = None
    else:
        data[["local_segment","timestamps_segment"]] = data["local_segment"].str.split(pat =";",n=1, expand=True)
    return(data)

# Each minute could fall into two segments.
# Firstly, we generate two rows for each resampled minute via resample_episodes rule:
# the first row's timestamp column is the start_timestamp, while the second row's timestamp column is the end_timestamp.
# Then, we check if the segments of start_timestamp are the same as the segments of end_timestamp:
# if they are the same (only fall into one segment), we will discard the second row;
# otherwise (fall into two segments), we will keep both.
def chunk_episodes(sensor_episodes):
    import copy
    import pandas as pd

    # Deduplicate episodes
    # Drop rows where segments of start_timestamp and end_timestamp are the same
    sensor_episodes = sensor_episodes.drop_duplicates(subset=["start_timestamp", "end_timestamp", "local_segment"], keep="first")

    # Delete useless columns
    for drop_col in ["utc_date_time", "local_date_time", "local_date", "local_time", "local_hour", "local_minute"]:
        del sensor_episodes[drop_col]
    
    # Avoid SettingWithCopyWarning
    sensor_episodes = sensor_episodes.copy()

    # Unix timestamp for current segment in milliseconds
    sensor_episodes[["segment_start_timestamp", "segment_end_timestamp"]] = sensor_episodes["timestamps_segment"].str.split(",", expand=True)

    # Compute chunked timestamp
    sensor_episodes["chunked_start_timestamp"] = sensor_episodes[["start_timestamp", "segment_start_timestamp"]].max(axis=1)
    sensor_episodes["chunked_end_timestamp"] = sensor_episodes[["end_timestamp", "segment_end_timestamp"]].min(axis=1)

    # Compute duration: intersection of current row and segment
    sensor_episodes["duration"] = (sensor_episodes["chunked_end_timestamp"] - sensor_episodes["chunked_start_timestamp"]) / (1000 * 60)

    # Compute chunked datetime
    sensor_episodes["chunked_start_datetime"] = pd.to_datetime(sensor_episodes["chunked_start_timestamp"], unit="ms", utc=True)
    sensor_episodes["chunked_start_datetime"] = pd.concat([data["chunked_start_datetime"].dt.tz_convert(tz) for tz, data in sensor_episodes.groupby("local_timezone")])

    sensor_episodes["chunked_end_datetime"] = pd.to_datetime(sensor_episodes["chunked_end_timestamp"], unit="ms", utc=True)
    sensor_episodes["chunked_end_datetime"] = pd.concat([data["chunked_end_datetime"].dt.tz_convert(tz) for tz, data in sensor_episodes.groupby("local_timezone")])

    # Merge episodes
    cols_for_groupby = [col for col in sensor_episodes.columns if col not in ["local_timezone", "timestamps_segment", "timestamp", "assigned_segments", "start_datetime", "end_datetime", "start_timestamp", "end_timestamp", "duration", "segment_start_timestamp", "segment_end_timestamp", "chunked_start_timestamp", "chunked_end_timestamp", "chunked_start_datetime", "chunked_end_datetime"]]

    sensor_episodes_grouped = sensor_episodes.groupby(by=cols_for_groupby)
    merged_sensor_episodes = sensor_episodes_grouped[["duration"]].sum()

    merged_sensor_episodes["start_timestamp"] = sensor_episodes_grouped["chunked_start_timestamp"].first()
    merged_sensor_episodes["end_timestamp"] = sensor_episodes_grouped["chunked_end_timestamp"].last()

    merged_sensor_episodes["local_start_date_time"] = sensor_episodes_grouped["chunked_start_datetime"].first()
    merged_sensor_episodes["local_end_date_time"] = sensor_episodes_grouped["chunked_end_datetime"].last()

    merged_sensor_episodes.reset_index(inplace=True)

    return merged_sensor_episodes

def fetch_provider_features(provider, provider_key, sensor_key, sensor_data_files, day_segments_file):
    import pandas as pd
    from importlib import import_module, util

    sensor_features = pd.DataFrame(columns=["local_segment"])
    day_segments_labels = pd.read_csv(day_segments_file, header=0)
    if "FEATURES" not in provider:
        raise ValueError("Provider config[{}][PROVIDERS][{}] is missing a FEATURES attribute in config.yaml".format(sensor_key.upper(), provider_key.upper()))

    if provider["COMPUTE"] == True:

            code_path =  sensor_key + "." + provider["SRC_FOLDER"] + ".main"
            feature_module = import_module(code_path)
            feature_function = getattr(feature_module,  provider["SRC_FOLDER"] + "_features")
            
            for day_segment in day_segments_labels["label"]:
                    print("{} Processing {} {} {}".format(rapids_log_tag, sensor_key, provider_key, day_segment))
                    features = feature_function(sensor_data_files, day_segment, provider, filter_data_by_segment=filter_data_by_segment, chunk_episodes=chunk_episodes)
                    sensor_features = sensor_features.merge(features, how="outer")
    else:
            for feature in provider["FEATURES"]:
                    sensor_features[feature] = None
    segment_colums = pd.DataFrame()
    split_segemnt_columns = sensor_features["local_segment"].str.split(pat="(.*)#(.*),(.*)", expand=True)
    new_segment_columns = split_segemnt_columns.iloc[:,1:4] if split_segemnt_columns.shape[1] == 5 else pd.DataFrame(columns=["local_segment_label", "local_segment_start_datetime","local_segment_end_datetime"])
    segment_colums[["local_segment_label", "local_segment_start_datetime", "local_segment_end_datetime"]] = new_segment_columns
    for i in range(segment_colums.shape[1]):
            sensor_features.insert(1 + i, segment_colums.columns[i], segment_colums[segment_colums.columns[i]])

    return sensor_features
