import pandas as pd
import numpy as np

def rapids_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):

    light_data = pd.read_csv(sensor_data_files["sensor_data"])
    requested_features = provider["FEATURES"]
    # name of the features this function can compute
    base_features_names = ["count", "maxlux", "minlux", "avglux", "medianlux", "stdlux"]
    # the subset of requested features this function can compute
    features_to_compute = list(set(requested_features) & set(base_features_names))

    light_features = pd.DataFrame(columns=["local_segment"] + features_to_compute)
    if not light_data.empty:
        light_data = filter_data_by_segment(light_data, time_segment)
        
        if not light_data.empty:
            light_features = pd.DataFrame()
            if "count" in features_to_compute:
                light_features["count"] = light_data.groupby(["local_segment"]).count()["timestamp"]
            
            # get light ambient luminance related features
            if "maxlux" in features_to_compute:
                light_features["maxlux"] = light_data.groupby(["local_segment"])["double_light_lux"].max()
            if "minlux" in features_to_compute:
                light_features["minlux"] = light_data.groupby(["local_segment"])["double_light_lux"].min()
            if "avglux" in features_to_compute:
                light_features["avglux"] = light_data.groupby(["local_segment"])["double_light_lux"].mean()
            if "medianlux" in features_to_compute:
                light_features["medianlux"] = light_data.groupby(["local_segment"])["double_light_lux"].median()
            if "stdlux" in features_to_compute:
                light_features["stdlux"] = light_data.groupby(["local_segment"])["double_light_lux"].std()
            
            light_features = light_features.reset_index()

    return light_features