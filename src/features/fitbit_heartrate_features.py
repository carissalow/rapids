import pandas as pd
import numpy as np
from scipy.stats import entropy
import json


def extractHRFeaturesFromSummaryData(heartrate_summary_data, daily_features_from_summary_data):
    heartrate_summary_features = pd.DataFrame()
    if "restinghr" in daily_features_from_summary_data:
        heartrate_summary_features["heartrate_daily_restinghr"] = heartrate_summary_data["heartrate_daily_restinghr"]
    # calories features might be inaccurate: they depend on users' fitbit profile (weight, height, etc.)
    if "caloriesoutofrange" in daily_features_from_summary_data:
        heartrate_summary_features["heartrate_daily_caloriesoutofrange"] = heartrate_summary_data["heartrate_daily_caloriesoutofrange"]
    if "caloriesfatburn" in daily_features_from_summary_data:
        heartrate_summary_features["heartrate_daily_caloriesfatburn"] = heartrate_summary_data["heartrate_daily_caloriesfatburn"]
    if "caloriescardio" in daily_features_from_summary_data:
        heartrate_summary_features["heartrate_daily_caloriescardio"] = heartrate_summary_data["heartrate_daily_caloriescardio"]
    if "caloriespeak" in daily_features_from_summary_data:
        heartrate_summary_features["heartrate_daily_caloriespeak"] = heartrate_summary_data["heartrate_daily_caloriespeak"]
    heartrate_summary_features.reset_index(inplace=True)
    
    return heartrate_summary_features

def extractHRFeaturesFromIntradayData(heartrate_intraday_data, features):
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
            for feature_name in list(set(["lengthoutofrange", "lengthfatburn", "lengthcardio", "lengthpeak"]) & set(features)):
                heartrate_zone = heartrate_intraday_data[heartrate_intraday_data["heartrate_zone"] == feature_name[6:]]
                heartrate_intraday_features["heartrate_" + day_segment + "_" + feature_name] = heartrate_zone.groupby(["local_date"])["device_id"].count() / num_rows_per_minute
                heartrate_intraday_features.fillna(value={"heartrate_" + day_segment + "_" + feature_name: 0}, inplace=True)
        heartrate_intraday_features.reset_index(inplace=True)

    return heartrate_intraday_features


heartrate_summary_data = pd.read_csv(snakemake.input["heartrate_summary_data"], index_col=["local_date"], parse_dates=["local_date"])
heartrate_intraday_data = pd.read_csv(snakemake.input["heartrate_intraday_data"], parse_dates=["local_date_time", "local_date"])
day_segment = snakemake.params["day_segment"]
features = snakemake.params["features"]
daily_features_from_summary_data = snakemake.params["daily_features_from_summary_data"]

heartrate_intraday_features = extractHRFeaturesFromIntradayData(heartrate_intraday_data, features)
if not heartrate_summary_data.empty and day_segment == "daily" and daily_features_from_summary_data != []:
    heartrate_summary_features = extractHRFeaturesFromSummaryData(heartrate_summary_data, daily_features_from_summary_data)
    heartrate_features = heartrate_intraday_features.merge(heartrate_summary_features, on=["local_date"], how="outer")
else:
    heartrate_features = heartrate_intraday_features

heartrate_features.to_csv(snakemake.output[0], index=False)
