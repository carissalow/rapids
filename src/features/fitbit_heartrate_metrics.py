import pandas as pd
import numpy as np
from scipy.stats import entropy
import json


heartrate_data = pd.read_csv(snakemake.input[0], parse_dates=["local_date_time", "local_date"])
day_segment = snakemake.params["day_segment"]
metrics = snakemake.params["metrics"]


heartrate_features = pd.DataFrame(columns=["local_date"] + ["heartrate_" + day_segment + "_" + x for x in metrics])
if not heartrate_data.empty:
    device_id = heartrate_data["device_id"][0]
    num_rows_per_minute = heartrate_data.groupby(["local_date", "local_hour", "local_minute"]).count().mean()["device_id"]
    if day_segment != "daily":
        heartrate_data =heartrate_data[heartrate_data["local_day_segment"] == day_segment]
    
    if not heartrate_data.empty:
        heartrate_features = pd.DataFrame()
     
        # get stats of heartrate
        if "maxhr" in metrics:
            heartrate_features["heartrate_" + day_segment + "_maxhr"] = heartrate_data.groupby(["local_date"])["heartrate"].max()
        if "minhr" in metrics:
            heartrate_features["heartrate_" + day_segment + "_minhr"] = heartrate_data.groupby(["local_date"])["heartrate"].min()
        if "avghr" in metrics:
            heartrate_features["heartrate_" + day_segment + "_avghr"] = heartrate_data.groupby(["local_date"])["heartrate"].mean()
        if "medianhr" in metrics:
            heartrate_features["heartrate_" + day_segment + "_medianhr"] = heartrate_data.groupby(["local_date"])["heartrate"].median()
        if "modehr" in metrics:
            heartrate_features["heartrate_" + day_segment + "_modehr"] = heartrate_data.groupby(["local_date"])["heartrate"].agg(pd.Series.mode)
        if "stdhr" in metrics:
            heartrate_features["heartrate_" + day_segment + "_stdhr"] = heartrate_data.groupby(["local_date"])["heartrate"].std()
        if "diffmaxmodehr" in metrics:
            heartrate_features["heartrate_" + day_segment + "_diffmaxmodehr"] = heartrate_data.groupby(["local_date"])["heartrate"].max() - heartrate_data.groupby(["local_date"])["heartrate"].agg(pd.Series.mode)
        if "diffminmodehr" in metrics:
            heartrate_features["heartrate_" + day_segment + "_diffminmodehr"] = heartrate_data.groupby(["local_date"])["heartrate"].agg(pd.Series.mode) - heartrate_data.groupby(["local_date"])["heartrate"].min()
        if "entropyhr" in metrics:
            heartrate_features["heartrate_" + day_segment + "_entropyhr"] = heartrate_data.groupby(["local_date"])["heartrate"].agg(entropy)

        # get number of minutes in each heart rate zone
        for feature_name in list(set(["lengthoutofrange", "lengthfatburn", "lengthcardio", "lengthpeak"]) & set(metrics)):
            heartrate_zone = heartrate_data[heartrate_data["heartrate_zone"] == feature_name[6:]]
            heartrate_features["heartrate_" + day_segment + "_" + feature_name] = heartrate_zone.groupby(["local_date"])["device_id"].count() / num_rows_per_minute

        heartrate_features = heartrate_features.reset_index()

heartrate_features.to_csv(snakemake.output[0], index=False)
