rapids_log_tag =  "RAPIDS:"

import os, sys
import importlib

def import_path(path):
    module_name = os.path.basename(path).replace('-', '_')
    spec = importlib.util.spec_from_loader(
        module_name,
        importlib.machinery.SourceFileLoader(module_name, path)
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sys.modules[module_name] = module
    return module

def filter_data_by_segment(data, time_segment):
    data.dropna(subset=["assigned_segments"], inplace=True)
    if(data.shape[0] == 0): # data is empty
        data["local_segment"] = data["timestamps_segment"] = None
        return data

    datetime_regex = "[0-9]{4}[\-|\/][0-9]{2}[\-|\/][0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}"
    timestamps_regex = "[0-9]{13}"
    segment_regex = "\[({}#{},{};{},{})\]".format(time_segment, datetime_regex, datetime_regex, timestamps_regex, timestamps_regex)
    data["local_segment"] = data["assigned_segments"].str.extract(segment_regex, expand=True)
    data = data.drop(columns=["assigned_segments"])
    data = data.dropna(subset = ["local_segment"])
    if(data.shape[0] == 0): # there are no rows belonging to time_segment after droping na
        data["timestamps_segment"] = None
    else:
        data[["local_segment","timestamps_segment"]] = data["local_segment"].str.split(pat =";",n=1, expand=True)
    
    # chunk episodes
    if (not data.empty) and ("start_timestamp" in data.columns) and ("end_timestamp" in data.columns):
        data = chunk_episodes(data)
    
    return data

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
    for drop_col in ["local_date_time", "local_date", "local_time", "local_hour", "local_minute"]:
        del sensor_episodes[drop_col]
    
    # Avoid SettingWithCopyWarning
    sensor_episodes = sensor_episodes.copy()

    # Unix timestamp for current segment in milliseconds
    sensor_episodes[["segment_start_timestamp", "segment_end_timestamp"]] = sensor_episodes["timestamps_segment"].str.split(",", expand=True).astype(int)

    # Compute chunked timestamp
    sensor_episodes["chunked_start_timestamp"] = sensor_episodes[["start_timestamp", "segment_start_timestamp"]].max(axis=1)
    sensor_episodes["chunked_end_timestamp"] = sensor_episodes[["end_timestamp", "segment_end_timestamp"]].min(axis=1)

    # Compute duration: intersection of current row and segment
    sensor_episodes["duration"] = (sensor_episodes["chunked_end_timestamp"] - sensor_episodes["chunked_start_timestamp"]) / (1000 * 60)

    # Merge episodes
    cols_for_groupby = [col for col in sensor_episodes.columns if col not in ["timestamps_segment", "timestamp", "assigned_segments", "start_datetime", "end_datetime", "start_timestamp", "end_timestamp", "duration", "chunked_start_timestamp", "chunked_end_timestamp"]]

    sensor_episodes_grouped = sensor_episodes.groupby(by=cols_for_groupby, sort=False, dropna=False)
    merged_sensor_episodes = sensor_episodes_grouped[["duration"]].sum()

    merged_sensor_episodes["start_timestamp"] = sensor_episodes_grouped["chunked_start_timestamp"].first()
    merged_sensor_episodes["end_timestamp"] = sensor_episodes_grouped["chunked_end_timestamp"].last()

    merged_sensor_episodes.reset_index(inplace=True)

    # Compute datetime
    merged_sensor_episodes["local_start_date_time"] = pd.to_datetime(merged_sensor_episodes["start_timestamp"], unit="ms", utc=True)
    merged_sensor_episodes["local_start_date_time"] = pd.concat([data["local_start_date_time"].dt.tz_convert(tz) for tz, data in merged_sensor_episodes.groupby("local_timezone")]).apply(lambda x: x.tz_localize(None).replace(microsecond=0))

    merged_sensor_episodes["local_end_date_time"] = pd.to_datetime(merged_sensor_episodes["end_timestamp"], unit="ms", utc=True)
    merged_sensor_episodes["local_end_date_time"] = pd.concat([data["local_end_date_time"].dt.tz_convert(tz) for tz, data in merged_sensor_episodes.groupby("local_timezone")]).apply(lambda x: x.tz_localize(None).replace(microsecond=0))

    return merged_sensor_episodes

def fetch_provider_features(provider, provider_key, sensor_key, sensor_data_files, time_segments_file):
    import pandas as pd
    from importlib import import_module, util

    sensor_features = pd.DataFrame(columns=["local_segment"])
    time_segments_labels = pd.read_csv(time_segments_file, header=0)
    if "FEATURES" not in provider:
        raise ValueError("Provider config[{}][PROVIDERS][{}] is missing a FEATURES attribute in config.yaml".format(sensor_key.upper(), provider_key.upper()))

    if provider["COMPUTE"] == True:

        feature_module = import_path(provider["SRC_SCRIPT"])
        feature_function = getattr(feature_module,  provider_key.lower() + "_features")
        
        if time_segments_labels["label"].empty:
            time_segments_labels["label"] = [""]
        for time_segment in time_segments_labels["label"]:
            print("{} Processing {} {} {}".format(rapids_log_tag, sensor_key, provider_key, time_segment))
            features = feature_function(sensor_data_files, time_segment, provider, filter_data_by_segment=filter_data_by_segment, chunk_episodes=chunk_episodes)
            if not "local_segment" in features.columns:
                raise ValueError("The dataframe returned by the " + sensor_key + " provider '" + provider_key + "' is missing the 'local_segment' column added by the 'filter_data_by_segment()' function. Check the provider script is using such function and is not removing 'local_segment' by accident (" + provider["SRC_SCRIPT"] + ")\n  The 'local_segment' column is used to index a provider's features (each row corresponds to a different time segment instance (e.g. 2020-01-01, 2020-01-02, 2020-01-03, etc.)")
            features.columns = ["{}{}".format("" if col.startswith("local_segment") else (sensor_key + "_"+ provider_key + "_"), col) for col in features.columns]
            sensor_features = pd.concat([sensor_features, features], axis=0, sort=False)
    else:
        for feature in provider["FEATURES"]:
            sensor_features[feature] = None
    segment_colums = pd.DataFrame()
    sensor_features['local_segment'] = sensor_features['local_segment'].str.replace(r'_RR\d+SS', '')
    split_segemnt_columns = sensor_features["local_segment"].str.split(pat="(.*)#(.*),(.*)", expand=True)
    new_segment_columns = split_segemnt_columns.iloc[:,1:4] if split_segemnt_columns.shape[1] == 5 else pd.DataFrame(columns=["local_segment_label", "local_segment_start_datetime","local_segment_end_datetime"])
    segment_colums[["local_segment_label", "local_segment_start_datetime", "local_segment_end_datetime"]] = new_segment_columns
    for i in range(segment_colums.shape[1]):
            sensor_features.insert(1 + i, segment_colums.columns[i], segment_colums[segment_colums.columns[i]])

    return sensor_features

def run_provider_cleaning_script(provider, provider_key, sensor_key, sensor_data_files):
    from importlib import import_module, util
    print("{} Processing {} {}".format(rapids_log_tag, sensor_key, provider_key))

    cleaning_module = import_path(provider["SRC_SCRIPT"])
    cleaning_function = getattr(cleaning_module,  provider_key.lower() + "_cleaning")
    sensor_features = cleaning_function(sensor_data_files, provider)

    return sensor_features
