import json
import pandas as pd
from datetime import datetime

STEPS_COLUMNS = ("device_id", "steps", "local_date_time", "timestamp")


def parseStepsData(steps_data):
    if steps_data.empty:
        return pd.DataFrame(columns=STEPS_COLUMNS)

    device_id = steps_data["device_id"].iloc[0]
    records = []
    
    # Parse JSON into individual records
    for record in steps_data.fitbit_data:
        record = json.loads(record)  # Parse text into JSON
        if "activities-steps" in record.keys():
            curr_date = datetime.strptime(record["activities-steps"][0]["dateTime"], "%Y-%m-%d")

            row_summary = (device_id,
                record["activities-steps"][0]["value"],
                curr_date,
                0)

            records.append(row_summary)
    
    parsed_data = pd.DataFrame(data=records, columns=STEPS_COLUMNS)

    return parsed_data


def main(json_raw, stream_parameters):
    parsed_data = parseStepsData(json_raw)
    parsed_data["timestamp"] = 0 # this column is added at readable_datetime.R because we neeed to take into account multiple timezones
    parsed_data['local_date_time'] = parsed_data['local_date_time'].dt.strftime('%Y-%m-%d %H:%M:%S')
    return(parsed_data)
