import pandas as pd
from scipy.stats import entropy

def statsFeatures(heartrate_data, features, features_type, heartrate_features):

    if features_type == "hr":
        col_name = "heartrate"
    elif features_type == "restinghr":
        col_name = "heartrate_daily_restinghr"
    elif features_type == "caloriesoutofrange":
        col_name = "heartrate_daily_caloriesoutofrange"
    elif features_type == "caloriesfatburn":
        col_name = "heartrate_daily_caloriesfatburn"
    elif features_type == "caloriescardio":
        col_name = "heartrate_daily_caloriescardio"
    elif features_type == "caloriespeak":
        col_name = "heartrate_daily_caloriespeak"
    else:
        raise ValueError("features_type can only be one of ['hr', 'restinghr', 'caloriesoutofrange', 'caloriesfatburn', 'caloriescardio', 'caloriespeak'].")

    if "sum" + features_type in features:
        heartrate_features["sum" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].sum()
    if "max" + features_type in features:
        heartrate_features["max" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].max()
    if "min" + features_type in features:
        heartrate_features["min" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].min()
    if "avg" + features_type in features:
        heartrate_features["avg" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].mean()
    if "median" + features_type in features:
        heartrate_features["median" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].median()
    if "mode" + features_type in features:
        heartrate_features["mode" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].agg(lambda x: pd.Series.mode(x)[0])
    if "std" + features_type in features:
        heartrate_features["std" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].std()
    if "diffmaxmode" + features_type in features:
        heartrate_features["diffmaxmode" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].max() - heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].agg(lambda x: pd.Series.mode(x)[0])
    if "diffminmode" + features_type in features:
        heartrate_features["diffminmode" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].agg(lambda x: pd.Series.mode(x)[0]) - heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].min()
    if "entropy" + features_type in features:
        heartrate_features["entropy" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].agg(entropy)
    
    return heartrate_features

def extractHRFeaturesFromIntradayData(heartrate_intraday_data, features, time_segment, filter_data_by_segment):
    heartrate_intraday_features = pd.DataFrame(columns=["local_segment"] + features)
    if not heartrate_intraday_data.empty:
        num_rows_per_minute = heartrate_intraday_data.groupby(["local_date", "local_hour", "local_minute"]).count().mean()["device_id"]
        heartrate_intraday_data = filter_data_by_segment(heartrate_intraday_data, time_segment)

        if not heartrate_intraday_data.empty:
            heartrate_intraday_features = pd.DataFrame()
        
            # get stats of heartrate
            heartrate_intraday_features = statsFeatures(heartrate_intraday_data, features, "hr", heartrate_intraday_features)

            # get number of minutes in each heart rate zone
            for feature_name in list(set(["minutesonoutofrangezone", "minutesonfatburnzone", "minutesoncardiozone", "minutesonpeakzone"]) & set(features)):
                heartrate_zone = heartrate_intraday_data[heartrate_intraday_data["heartrate_zone"] == feature_name[9:-4]]
                heartrate_intraday_features[feature_name] = heartrate_zone.groupby(["local_segment"])["device_id"].count() / num_rows_per_minute
                heartrate_intraday_features.fillna(value={feature_name: 0}, inplace=True)
            
            heartrate_intraday_features.reset_index(inplace=True)

    return heartrate_intraday_features


def rapids_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):

    heartrate_intraday_data = pd.read_csv(sensor_data_files["sensor_data"])

    requested_intraday_features = provider["FEATURES"]
    # name of the features this function can compute
    base_intraday_features_names = ["maxhr", "minhr", "avghr", "medianhr", "modehr", "stdhr", "diffmaxmodehr", "diffminmodehr", "entropyhr", "minutesonoutofrangezone", "minutesonfatburnzone", "minutesoncardiozone", "minutesonpeakzone"]
    # the subset of requested features this function can compute
    intraday_features_to_compute = list(set(requested_intraday_features) & set(base_intraday_features_names))
    
    # extract features from intraday data
    heartrate_intraday_features = extractHRFeaturesFromIntradayData(heartrate_intraday_data, intraday_features_to_compute, time_segment, filter_data_by_segment)
    
    return heartrate_intraday_features
