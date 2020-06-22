import pandas as pd
from scipy.stats import entropy

def extractHRFeaturesFromSummaryData(heartrate_summary_data, summary_features):
    heartrate_summary_features = pd.DataFrame()
    if "restinghr" in summary_features:
        heartrate_summary_features["heartrate_daily_restinghr"] = heartrate_summary_data["heartrate_daily_restinghr"]
    # calories features might be inaccurate: they depend on users' fitbit profile (weight, height, etc.)
    if "caloriesoutofrange" in summary_features:
        heartrate_summary_features["heartrate_daily_caloriesoutofrange"] = heartrate_summary_data["heartrate_daily_caloriesoutofrange"]
    if "caloriesfatburn" in summary_features:
        heartrate_summary_features["heartrate_daily_caloriesfatburn"] = heartrate_summary_data["heartrate_daily_caloriesfatburn"]
    if "caloriescardio" in summary_features:
        heartrate_summary_features["heartrate_daily_caloriescardio"] = heartrate_summary_data["heartrate_daily_caloriescardio"]
    if "caloriespeak" in summary_features:
        heartrate_summary_features["heartrate_daily_caloriespeak"] = heartrate_summary_data["heartrate_daily_caloriespeak"]
    heartrate_summary_features.reset_index(inplace=True)
    
    return heartrate_summary_features

def extractHRFeaturesFromIntradayData(heartrate_intraday_data, features, day_segment):
    heartrate_intraday_features = pd.DataFrame(columns=["local_date"] + ["heartrate_" + day_segment + "_" + x for x in features])
    if not heartrate_intraday_data.empty:
        device_id = heartrate_intraday_data["device_id"][0]
        num_rows_per_minute = heartrate_intraday_data.groupby(["local_date", "local_hour", "local_minute"]).count().mean()["device_id"]
        if day_segment != "daily":
            heartrate_intraday_data = heartrate_intraday_data[heartrate_intraday_data["local_day_segment"] == day_segment]

        if not heartrate_intraday_data.empty:
            heartrate_intraday_features = pd.DataFrame()
        
            # get stats of heartrate
            if "maxhr" in features:
                heartrate_intraday_features["heartrate_" + day_segment + "_maxhr"] = heartrate_intraday_data[["local_date", "heartrate"]].groupby(["local_date"])["heartrate"].max()
            if "minhr" in features:
                heartrate_intraday_features["heartrate_" + day_segment + "_minhr"] = heartrate_intraday_data[["local_date", "heartrate"]].groupby(["local_date"])["heartrate"].min()
            if "avghr" in features:
                heartrate_intraday_features["heartrate_" + day_segment + "_avghr"] = heartrate_intraday_data[["local_date", "heartrate"]].groupby(["local_date"])["heartrate"].mean()
            if "medianhr" in features:
                heartrate_intraday_features["heartrate_" + day_segment + "_medianhr"] = heartrate_intraday_data[["local_date", "heartrate"]].groupby(["local_date"])["heartrate"].median()
            if "modehr" in features:
                heartrate_intraday_features["heartrate_" + day_segment + "_modehr"] = heartrate_intraday_data[["local_date", "heartrate"]].groupby(["local_date"])["heartrate"].agg(lambda x: pd.Series.mode(x)[0])
            if "stdhr" in features:
                heartrate_intraday_features["heartrate_" + day_segment + "_stdhr"] = heartrate_intraday_data[["local_date", "heartrate"]].groupby(["local_date"])["heartrate"].std()
            if "diffmaxmodehr" in features:
                heartrate_intraday_features["heartrate_" + day_segment + "_diffmaxmodehr"] = heartrate_intraday_data[["local_date", "heartrate"]].groupby(["local_date"])["heartrate"].max() - heartrate_intraday_data[["local_date", "heartrate"]].groupby(["local_date"])["heartrate"].agg(lambda x: pd.Series.mode(x)[0])
            if "diffminmodehr" in features:
                heartrate_intraday_features["heartrate_" + day_segment + "_diffminmodehr"] = heartrate_intraday_data[["local_date", "heartrate"]].groupby(["local_date"])["heartrate"].agg(lambda x: pd.Series.mode(x)[0]) - heartrate_intraday_data[["local_date", "heartrate"]].groupby(["local_date"])["heartrate"].min()
            if "entropyhr" in features:
                heartrate_intraday_features["heartrate_" + day_segment + "_entropyhr"] = heartrate_intraday_data[["local_date", "heartrate"]].groupby(["local_date"])["heartrate"].agg(entropy)

            # get number of minutes in each heart rate zone
            for feature_name in list(set(["minutesonoutofrangezone", "minutesonfatburnzone", "minutesoncardiozone", "minutesonpeakzone"]) & set(features)):
                heartrate_zone = heartrate_intraday_data[heartrate_intraday_data["heartrate_zone"] == feature_name[9:-4]]
                heartrate_intraday_features["heartrate_" + day_segment + "_" + feature_name] = heartrate_zone.groupby(["local_date"])["device_id"].count() / num_rows_per_minute
                heartrate_intraday_features.fillna(value={"heartrate_" + day_segment + "_" + feature_name: 0}, inplace=True)
        heartrate_intraday_features.reset_index(inplace=True)

    return heartrate_intraday_features

def base_fitbit_heartrate_features(heartrate_summary_data, heartrate_intraday_data, day_segment, requested_summary_features, requested_intraday_features):
    # name of the features this function can compute
    base_summary_features_names = ["restinghr", "caloriesoutofrange", "caloriesfatburn", "caloriescardio", "caloriespeak"]
    base_intraday_features_names = ["maxhr", "minhr", "avghr", "medianhr", "modehr", "stdhr", "diffmaxmodehr", "diffminmodehr", "entropyhr", "minutesonoutofrangezone", "minutesonfatburnzone", "minutesoncardiozone", "minutesonpeakzone"]
    # the subset of requested features this function can compute
    summary_features_to_compute = list(set(requested_summary_features) & set(base_summary_features_names))
    intraday_features_to_compute = list(set(requested_intraday_features) & set(base_intraday_features_names))

    heartrate_intraday_features = extractHRFeaturesFromIntradayData(heartrate_intraday_data, intraday_features_to_compute, day_segment)
    if not heartrate_summary_data.empty and day_segment == "daily" and summary_features_to_compute != []:
        heartrate_summary_features = extractHRFeaturesFromSummaryData(heartrate_summary_data, summary_features_to_compute)
        heartrate_features = heartrate_intraday_features.merge(heartrate_summary_features, on=["local_date"], how="outer")
    else:
        heartrate_features = heartrate_intraday_features

    return heartrate_features
