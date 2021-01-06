import json
import numpy as np
import pandas as pd
from datetime import datetime


CALORIES_INTRADAY_COLUMNS = ("device_id",
                                "level", "mets", "value",
                                "local_date_time", "timestamp")

def parseCaloriesData(calories_data):
    if calories_data.empty:
        return pd.DataFrame(), pd.DataFrame(columns=CALORIES_INTRADAY_COLUMNS)
    device_id = calories_data["device_id"].iloc[0]
    records_intraday = []
    # Parse JSON into individual records
    for record in calories_data.fitbit_data:
        record = json.loads(record)  # Parse text into JSON
        curr_date = datetime.strptime(
            record["activities-calories"][0]["dateTime"], "%Y-%m-%d")
        dataset = record["activities-calories-intraday"]["dataset"]
        for data in dataset:
            d_time = datetime.strptime(data["time"], '%H:%M:%S').time()
            d_datetime = datetime.combine(curr_date, d_time)

            row_intraday = (device_id,
                            data["level"], data["mets"], data["value"],
                            d_datetime, 0)

            records_intraday.append(row_intraday)

    return pd.DataFrame(data=[], columns=["local_date_time", "timestamp"]), pd.DataFrame(data=records_intraday, columns=CALORIES_INTRADAY_COLUMNS)

table_format = snakemake.params["table_format"]
timezone = snakemake.params["timezone"]

if table_format == "JSON":
    json_raw = pd.read_csv(snakemake.input[0])
    summary, intraday = parseCaloriesData(json_raw)
elif table_format == "CSV":
    summary = pd.read_csv(snakemake.input[0], parse_dates=["local_date_time"], date_parser=lambda col: pd.to_datetime(col).tz_localize(None))
    intraday = pd.read_csv(snakemake.input[1], parse_dates=["local_date_time"], date_parser=lambda col: pd.to_datetime(col).tz_localize(None))

#    if not pd.isnull(local_start_date) and not pd.isnull(local_end_date):

if summary.shape[0] > 0:
    summary["timestamp"] = summary["local_date_time"].dt.tz_localize(timezone, ambiguous=False, nonexistent="NaT").dropna().astype(np.int64) // 10**6
    summary.dropna(subset=['timestamp'], inplace=True)
if intraday.shape[0] > 0:
    intraday["timestamp"] = intraday["local_date_time"].dt.tz_localize(timezone, ambiguous=False, nonexistent="NaT").dropna().astype(np.int64) // 10**6
    intraday.dropna(subset=['timestamp'], inplace=True)

summary.to_csv(snakemake.output["summary_data"], index=False)
intraday.to_csv(snakemake.output["intraday_data"], index=False)