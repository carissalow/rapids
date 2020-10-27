import pandas as pd
from scipy.stats import entropy

def extractHRFeaturesFromSummaryData(heartrate_summary_data, summary_features):
    heartrate_summary_data.set_index("local_segment", inplace=True)
    heartrate_summary_features = pd.DataFrame()
    if "restinghr" in summary_features:
        heartrate_summary_features["heartrate_rapids_restinghr"] = heartrate_summary_data["heartrate_daily_restinghr"]
    # calories features might be inaccurate: they depend on users' fitbit profile (weight, height, etc.)
    if "caloriesoutofrange" in summary_features:
        heartrate_summary_features["heartrate_rapids_caloriesoutofrange"] = heartrate_summary_data["heartrate_daily_caloriesoutofrange"]
    if "caloriesfatburn" in summary_features:
        heartrate_summary_features["heartrate_rapids_caloriesfatburn"] = heartrate_summary_data["heartrate_daily_caloriesfatburn"]
    if "caloriescardio" in summary_features:
        heartrate_summary_features["heartrate_rapids_caloriescardio"] = heartrate_summary_data["heartrate_daily_caloriescardio"]
    if "caloriespeak" in summary_features:
        heartrate_summary_features["heartrate_rapids_caloriespeak"] = heartrate_summary_data["heartrate_daily_caloriespeak"]
    heartrate_summary_features.reset_index(inplace=True)
    
    return heartrate_summary_features

def extractHRFeaturesFromIntradayData(heartrate_intraday_data, features, day_segment, filter_data_by_segment):
    heartrate_intraday_features = pd.DataFrame(columns=["local_segment"] + ["heartrate_rapids_" + x for x in features])
    if not heartrate_intraday_data.empty:
        num_rows_per_minute = heartrate_intraday_data.groupby(["local_date", "local_hour", "local_minute"]).count().mean()["device_id"]
        heartrate_intraday_data = filter_data_by_segment(heartrate_intraday_data, day_segment)

        if not heartrate_intraday_data.empty:
            heartrate_intraday_features = pd.DataFrame()
        
            # get stats of heartrate
            if "maxhr" in features:
                heartrate_intraday_features["heartrate_rapids_maxhr"] = heartrate_intraday_data[["local_segment", "heartrate"]].groupby(["local_segment"])["heartrate"].max()
            if "minhr" in features:
                heartrate_intraday_features["heartrate_rapids_minhr"] = heartrate_intraday_data[["local_segment", "heartrate"]].groupby(["local_segment"])["heartrate"].min()
            if "avghr" in features:
                heartrate_intraday_features["heartrate_rapids_avghr"] = heartrate_intraday_data[["local_segment", "heartrate"]].groupby(["local_segment"])["heartrate"].mean()
            if "medianhr" in features:
                heartrate_intraday_features["heartrate_rapids_medianhr"] = heartrate_intraday_data[["local_segment", "heartrate"]].groupby(["local_segment"])["heartrate"].median()
            if "modehr" in features:
                heartrate_intraday_features["heartrate_rapids_modehr"] = heartrate_intraday_data[["local_segment", "heartrate"]].groupby(["local_segment"])["heartrate"].agg(lambda x: pd.Series.mode(x)[0])
            if "stdhr" in features:
                heartrate_intraday_features["heartrate_rapids_stdhr"] = heartrate_intraday_data[["local_segment", "heartrate"]].groupby(["local_segment"])["heartrate"].std()
            if "diffmaxmodehr" in features:
                heartrate_intraday_features["heartrate_rapids_diffmaxmodehr"] = heartrate_intraday_data[["local_segment", "heartrate"]].groupby(["local_segment"])["heartrate"].max() - heartrate_intraday_data[["local_segment", "heartrate"]].groupby(["local_segment"])["heartrate"].agg(lambda x: pd.Series.mode(x)[0])
            if "diffminmodehr" in features:
                heartrate_intraday_features["heartrate_rapids_diffminmodehr"] = heartrate_intraday_data[["local_segment", "heartrate"]].groupby(["local_segment"])["heartrate"].agg(lambda x: pd.Series.mode(x)[0]) - heartrate_intraday_data[["local_segment", "heartrate"]].groupby(["local_segment"])["heartrate"].min()
            if "entropyhr" in features:
                heartrate_intraday_features["heartrate_rapids_entropyhr"] = heartrate_intraday_data[["local_segment", "heartrate"]].groupby(["local_segment"])["heartrate"].agg(entropy)

            # get number of minutes in each heart rate zone
            for feature_name in list(set(["minutesonoutofrangezone", "minutesonfatburnzone", "minutesoncardiozone", "minutesonpeakzone"]) & set(features)):
                heartrate_zone = heartrate_intraday_data[heartrate_intraday_data["heartrate_zone"] == feature_name[9:-4]]
                heartrate_intraday_features["heartrate_rapids_" + feature_name] = heartrate_zone.groupby(["local_segment"])["device_id"].count() / num_rows_per_minute
                heartrate_intraday_features.fillna(value={"heartrate_rapids_" + feature_name: 0}, inplace=True)
        heartrate_intraday_features.reset_index(inplace=True)

    return heartrate_intraday_features


def rapids_features(sensor_data_files, day_segment, provider, filter_data_by_segment, *args, **kwargs):

    heartrate_summary_data = pd.read_csv(sensor_data_files["sensor_data"][0])
    heartrate_intraday_data = pd.read_csv(sensor_data_files["sensor_data"][1])

    requested_summary_features = provider["FEATURES"]["SUMMARY"]
    requested_intraday_features = provider["FEATURES"]["INTRADAY"]
    # name of the features this function can compute
    base_summary_features_names = ["restinghr", "caloriesoutofrange", "caloriesfatburn", "caloriescardio", "caloriespeak"]
    base_intraday_features_names = ["maxhr", "minhr", "avghr", "medianhr", "modehr", "stdhr", "diffmaxmodehr", "diffminmodehr", "entropyhr", "minutesonoutofrangezone", "minutesonfatburnzone", "minutesoncardiozone", "minutesonpeakzone"]
    # the subset of requested features this function can compute
    summary_features_to_compute = list(set(requested_summary_features) & set(base_summary_features_names))
    intraday_features_to_compute = list(set(requested_intraday_features) & set(base_intraday_features_names))

    heartrate_intraday_features = extractHRFeaturesFromIntradayData(heartrate_intraday_data, intraday_features_to_compute, day_segment, filter_data_by_segment)
    if not heartrate_summary_data.empty and day_segment == "daily" and summary_features_to_compute != []:
        # filter by segment and skipping any non-daily segment
        heartrate_summary_data = filter_data_by_segment(heartrate_summary_data, "daily")
        
        datetime_start_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2} 00:00:00"
        datetime_end_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2} 23:59:59"

        segment_regex = "daily#{},{}".format(datetime_start_regex, datetime_end_regex)
        heartrate_summary_data = heartrate_summary_data[heartrate_summary_data["local_segment"].str.match(segment_regex)]

        # extract daily features from summary data
        heartrate_summary_features = extractHRFeaturesFromSummaryData(heartrate_summary_data, summary_features_to_compute)

        # merge summary features and intraday features
        heartrate_features = heartrate_intraday_features.merge(heartrate_summary_features, on=["local_segment"], how="outer")
    else:
        heartrate_features = heartrate_intraday_features
    
    return heartrate_features
