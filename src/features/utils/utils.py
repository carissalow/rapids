rapids_log_tag =  "RAPIDS:"

def filter_data_by_segment(data, day_segment):
    date_regex = "[0-9]{4}[\-|\/][0-9]{2}[\-|\/][0-9]{2}"
    hour_regex = "[0-9]{2}:[0-9]{2}:[0-9]{2}"
    segment_regex = "\[({}#{}#{}#{}#{})\]".format(day_segment, date_regex, hour_regex, date_regex, hour_regex)
    data["local_segment"] = data["assigned_segments"].str.extract(segment_regex, expand=True)
    return(data.dropna(subset = ["local_segment"]))

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
                    print("---")
                    features = feature_function(sensor_data, day_segment, provider, filter_data_by_segment=filter_data_by_segment)
                    print("2")
                    sensor_features = sensor_features.merge(features, how="outer")
    else:
            for feature in provider["FEATURES"]:
                    sensor_features[feature] = None
    print("3")
    segment_colums = pd.DataFrame()
    split_segemnt_columns = sensor_features["local_segment"].str.split(pat="#", expand=True)
    new_segment_columns = split_segemnt_columns if split_segemnt_columns.shape[1] == 5 else pd.DataFrame(columns=["local_segment_label", "local_start_date", "local_start_time", "local_end_date", "local_end_time"])
    segment_colums[["local_segment_label", "local_start_date", "local_start_time", "local_end_date", "local_end_time"]] = new_segment_columns
    for i in range(segment_colums.shape[1]):
            sensor_features.insert(1 + i, segment_colums.columns[i], segment_colums[segment_colums.columns[i]])

    return sensor_features