import json
import pandas as pd
from datetime import datetime

CALORIES_INTRADAY_COLUMNS = ("device_id", "level", "mets", "value", "local_date_time", "timestamp")

def parseCaloriesData(calories_data):
    if calories_data.empty:
        return pd.DataFrame(columns=CALORIES_INTRADAY_COLUMNS)
    device_id = calories_data["device_id"].iloc[0]
    records_intraday = []

    # Parse JSON into individual records
    for record in calories_data.json_fitbit_column:
        record = json.loads(record)  # Parse text into JSON
        if "activities-calories" in record and "activities-calories-intraday" in record:
            curr_date = datetime.strptime(record["activities-calories"][0]["dateTime"], "%Y-%m-%d")
            dataset = record["activities-calories-intraday"]["dataset"]
            for data in dataset:
                d_time = datetime.strptime(data["time"], '%H:%M:%S').time()
                d_datetime = datetime.combine(curr_date, d_time).strftime("%Y-%m-%d %H:%M:%S")
                row_intraday = (device_id, data["level"], data["mets"], data["value"], d_datetime, 0)
                records_intraday.append(row_intraday)

    return pd.DataFrame(data=records_intraday, columns=CALORIES_INTRADAY_COLUMNS)

def main(json_raw, stream_parameters):
    parsed_data = parseCaloriesData(json_raw)
    parsed_data["mets"] = parsed_data["mets"] / 10
    return parsed_data
