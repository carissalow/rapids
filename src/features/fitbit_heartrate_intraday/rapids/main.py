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

    if "intradaysum" + features_type in features:
        heartrate_features["heartrate_rapids_intradaysum" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].sum()
    if "intradaymax" + features_type in features:
        heartrate_features["heartrate_rapids_intradaymax" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].max()
    if "intradaymin" + features_type in features:
        heartrate_features["heartrate_rapids_intradaymin" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].min()
    if "intradayavg" + features_type in features:
        heartrate_features["heartrate_rapids_intradayavg" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].mean()
    if "intradaymedian" + features_type in features:
        heartrate_features["heartrate_rapids_intradaymedian" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].median()
    if "intradaymode" + features_type in features:
        heartrate_features["heartrate_rapids_intradaymode" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].agg(lambda x: pd.Series.mode(x)[0])
    if "intradaystd" + features_type in features:
        heartrate_features["heartrate_rapids_intradaystd" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].std()
    if "intradaydiffmaxmode" + features_type in features:
        heartrate_features["heartrate_rapids_intradaydiffmaxmode" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].max() - heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].agg(lambda x: pd.Series.mode(x)[0])
    if "intradaydiffminmode" + features_type in features:
        heartrate_features["heartrate_rapids_intradaydiffminmode" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].agg(lambda x: pd.Series.mode(x)[0]) - heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].min()
    if "intradayentropy" + features_type in features:
        heartrate_features["heartrate_rapids_intradayentropy" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].agg(entropy)
    
    return heartrate_features

def extractHRFeaturesFromIntradayData(heartrate_intraday_data, features, day_segment, filter_data_by_segment):
    heartrate_intraday_features = pd.DataFrame(columns=["local_segment"] + ["heartrate_rapids_" + x for x in features])
    if not heartrate_intraday_data.empty:
        num_rows_per_minute = heartrate_intraday_data.groupby(["local_date", "local_hour", "local_minute"]).count().mean()["device_id"]
        heartrate_intraday_data = filter_data_by_segment(heartrate_intraday_data, day_segment)

        if not heartrate_intraday_data.empty:
            heartrate_intraday_features = pd.DataFrame()
        
            # get stats of heartrate
            heartrate_intraday_features = statsFeatures(heartrate_intraday_data, features, "hr", heartrate_intraday_features)

            # get number of minutes in each heart rate zone
            for feature_name in list(set(["intradayminutesonoutofrangezone", "intradayminutesonfatburnzone", "intradayminutesoncardiozone", "intradayminutesonpeakzone"]) & set(features)):
                heartrate_zone = heartrate_intraday_data[heartrate_intraday_data["heartrate_zone"] == feature_name[17:-4]]
                heartrate_intraday_features["heartrate_rapids_" + feature_name] = heartrate_zone.groupby(["local_segment"])["device_id"].count() / num_rows_per_minute
                heartrate_intraday_features.fillna(value={"heartrate_rapids_" + feature_name: 0}, inplace=True)
        heartrate_intraday_features.reset_index(inplace=True)

    return heartrate_intraday_features


def rapids_features(sensor_data_files, day_segment, provider, filter_data_by_segment, *args, **kwargs):

    heartrate_intraday_data = pd.read_csv(sensor_data_files["sensor_data"])

    requested_intraday_features = ["intraday" + x for x in provider["FEATURES"]]
    # name of the features this function can compute
    base_intraday_features_names = ["intradaymaxhr", "intradayminhr", "intradayavghr", "intradaymedianhr", "intradaymodehr", "intradaystdhr", "intradaydiffmaxmodehr", "intradaydiffminmodehr", "intradayentropyhr", "intradayminutesonoutofrangezone", "intradayminutesonfatburnzone", "intradayminutesoncardiozone", "intradayminutesonpeakzone"]
    # the subset of requested features this function can compute
    intraday_features_to_compute = list(set(requested_intraday_features) & set(base_intraday_features_names))
    
    # extract features from intraday data
    heartrate_intraday_features = extractHRFeaturesFromIntradayData(heartrate_intraday_data, intraday_features_to_compute, day_segment, filter_data_by_segment)
    
    return heartrate_intraday_features
