import pandas as pd
import numpy as np

def rapids_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):

    acc_data = pd.read_csv(sensor_data_files["sensor_data"])
    requested_features = provider["FEATURES"]
    # name of the features this function can compute
    base_features_names = ["maxmagnitude", "minmagnitude", "avgmagnitude", "medianmagnitude", "stdmagnitude"]
    # the subset of requested features this function can compute
    features_to_compute = list(set(requested_features) & set(base_features_names))

    acc_features = pd.DataFrame(columns=["local_segment"] + features_to_compute)
    if not acc_data.empty:
        acc_data = filter_data_by_segment(acc_data, time_segment)
        
        if not acc_data.empty:
            acc_features = pd.DataFrame()
            # get magnitude related features: magnitude = sqrt(x^2+y^2+z^2)
            magnitude = acc_data.apply(lambda row: np.sqrt(row["double_values_0"] ** 2 + row["double_values_1"] ** 2 + row["double_values_2"] ** 2), axis=1)
            acc_data = acc_data.assign(magnitude = magnitude.values)
            
            if "maxmagnitude" in features_to_compute:
                acc_features["maxmagnitude"] = acc_data.groupby(["local_segment"])["magnitude"].max()
            if "minmagnitude" in features_to_compute:
                acc_features["minmagnitude"] = acc_data.groupby(["local_segment"])["magnitude"].min()
            if "avgmagnitude" in features_to_compute:
                acc_features["avgmagnitude"] = acc_data.groupby(["local_segment"])["magnitude"].mean()
            if "medianmagnitude" in features_to_compute:
                acc_features["medianmagnitude"] = acc_data.groupby(["local_segment"])["magnitude"].median()
            if "stdmagnitude" in features_to_compute:
                acc_features["stdmagnitude"] = acc_data.groupby(["local_segment"])["magnitude"].std()
            
            acc_features = acc_features.reset_index()

    return acc_features
