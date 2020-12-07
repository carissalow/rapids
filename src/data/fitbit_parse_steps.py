import json, yaml
import pandas as pd
import numpy as np
from datetime import datetime, timezone
from math import trunc

STEPS_COLUMNS = ("device_id", "steps", "local_date_time", "timestamp")


def parseStepsData(steps_data, fitbit_data_type):
    if steps_data.empty:
        return pd.DataFrame(), pd.DataFrame(columns=STEPS_INTRADAY_COLUMNS)
    device_id = steps_data["device_id"].iloc[0]
    records_summary, records_intraday = [], []
    
    # Parse JSON into individual records
    for record in steps_data.fitbit_data:
        record = json.loads(record)  # Parse text into JSON
        curr_date = datetime.strptime(record["activities-steps"][0]["dateTime"], "%Y-%m-%d")
        
        # Parse summary data
        if fitbit_data_type == "summary":

            row_summary = (device_id,
                record["activities-steps"][0]["value"],
                curr_date,
                0)
            
            records_summary.append(row_summary)

        # Parse intraday data
        if fitbit_data_type == "intraday":
            dataset = record["activities-steps-intraday"]["dataset"]
            for data in dataset:
                d_time = datetime.strptime(data["time"], '%H:%M:%S').time()
                d_datetime = datetime.combine(curr_date, d_time)

                row_intraday = (device_id,
                    data["value"],
                    d_datetime,
                    0)

                records_intraday.append(row_intraday)
    
    if fitbit_data_type == "summary":
        parsed_data = pd.DataFrame(data=records_summary, columns=STEPS_COLUMNS)
    elif fitbit_data_type == "intraday":
        parsed_data = pd.DataFrame(data=records_intraday, columns=STEPS_COLUMNS)
    else:
        raise ValueError("fitbit_data_type can only be one of ['summary', 'intraday'].")

    return parsed_data



timezone = snakemake.params["timezone"]
column_format = snakemake.params["column_format"]
fitbit_data_type = snakemake.params["fitbit_data_type"]

with open(snakemake.input["participant_file"], "r", encoding="utf-8") as f:
    participant_file = yaml.safe_load(f)
local_start_date = pd.Timestamp(participant_file["FITBIT"]["START_DATE"])
local_end_date = pd.Timestamp(participant_file["FITBIT"]["END_DATE"]) + pd.DateOffset(1)

if column_format == "JSON":
    json_raw = pd.read_csv(snakemake.input["raw_data"])
    parsed_data = parseStepsData(json_raw, fitbit_data_type)
elif column_format == "PLAIN_TEXT":
    parsed_data = pd.read_csv(snakemake.input["raw_data"], parse_dates=["local_date_time"], date_parser=lambda col: pd.to_datetime(col).tz_localize(None))
else:
    raise ValueError("column_format can only be one of ['JSON', 'PLAIN_TEXT'].")

# Only keep dates in the range of [local_start_date, local_end_date)
parsed_data = parsed_data.loc[(parsed_data["local_date_time"] >= local_start_date) & (parsed_data["local_date_time"] < local_end_date)]

if parsed_data.shape[0] > 0:
    parsed_data["timestamp"] = parsed_data["local_date_time"].dt.tz_localize(timezone).astype(np.int64) // 10**6

parsed_data.to_csv(snakemake.output[0], index=False)
