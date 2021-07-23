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
    for record in steps_data.json_fitbit_column:
        record = json.loads(record)  # Parse text into JSON
        if "activities-steps" in record.keys():
            curr_date = datetime.strptime(record["activities-steps"][0]["dateTime"], "%Y-%m-%d")

            # Parse intraday data
            if "activities-steps-intraday" in record.keys():
                dataset = record["activities-steps-intraday"]["dataset"]
                for data in dataset:
                    d_time = datetime.strptime(data["time"], '%H:%M:%S').time()
                    d_datetime = datetime.combine(curr_date, d_time).strftime("%Y-%m-%d %H:%M:%S")

                    row_intraday = (device_id,
                        data["value"],
                        d_datetime,
                        0)

                    records.append(row_intraday)
    
    parsed_data = pd.DataFrame(data=records, columns=STEPS_COLUMNS)

    return parsed_data


def main(json_raw, stream_parameters):
    parsed_data = parseStepsData(json_raw)
    return parsed_data
