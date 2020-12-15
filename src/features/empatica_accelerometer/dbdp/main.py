import pandas as pd
import numpy as np

def dbdp_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):

    sensor_data = pd.read_csv(sensor_data_files["sensor_data"])
    requested_features = provider["FEATURES"]
    # name of the features this function can compute
    base_features_names = [] # ["maxmagnitude", "minmagnitude", "avgmagnitude", "medianmagnitude", "stdmagnitude"]
    # the subset of requested features this function can compute
    features_to_compute = list(set(requested_features) & set(base_features_names))

    features = pd.DataFrame(columns=["local_segment"] + features_to_compute)
    if not sensor_data.empty:
        sensor_data = filter_data_by_segment(sensor_data, time_segment)
        
        if not sensor_data.empty:
            features = pd.DataFrame()
            

    return features