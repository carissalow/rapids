import pandas as pd

def base_light_features(light_data, day_segment, requested_features):
    # name of the features this function can compute
    base_features_names = ["count", "maxlux", "minlux", "avglux", "medianlux", "stdlux"]
    # the subset of requested features this function can compute
    features_to_compute = list(set(requested_features) & set(base_features_names))

    if light_data.empty:
        light_features = pd.DataFrame(columns=["local_date"] + ["light_" + day_segment + "_" + x for x in features_to_compute])
    else:
        if day_segment != "daily":
            light_data =light_data[light_data["local_day_segment"] == day_segment]
        
        if not light_data.empty:
            light_features = pd.DataFrame()
            if "count" in features_to_compute:
                light_features["light_" + day_segment + "_count"] = light_data.groupby(["local_date"]).count()["timestamp"]
            
            # get light ambient luminance related features
            if "maxlux" in features_to_compute:
                light_features["light_" + day_segment + "_maxlux"] = light_data.groupby(["local_date"])["double_light_lux"].max()
            if "minlux" in features_to_compute:
                light_features["light_" + day_segment + "_minlux"] = light_data.groupby(["local_date"])["double_light_lux"].min()
            if "avglux" in features_to_compute:
                light_features["light_" + day_segment + "_avglux"] = light_data.groupby(["local_date"])["double_light_lux"].mean()
            if "medianlux" in features_to_compute:
                light_features["light_" + day_segment + "_medianlux"] = light_data.groupby(["local_date"])["double_light_lux"].median()
            if "stdlux" in features_to_compute:
                light_features["light_" + day_segment + "_stdlux"] = light_data.groupby(["local_date"])["double_light_lux"].std()
            
            light_features = light_features.reset_index()

    return light_features