import json, yaml
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import dateutil.parser

SLEEP_SUMMARY_COLUMNS = ("device_id", "efficiency",
                                "minutes_after_wakeup", "minutes_asleep", "minutes_awake", "minutes_to_fall_asleep", "minutes_in_bed",
                                "is_main_sleep", "type",
                                "local_start_date_time", "local_end_date_time",
                                "timestamp")


# Parse one record for sleep API version 1.2
def parseOneSleepRecord(record, device_id, d_is_main_sleep, records_summary, episode_type):
    
    sleep_record_type = episode_type

    d_start_datetime = datetime.strptime(record["startTime"][:18], "%Y-%m-%dT%H:%M:%S")
    d_end_datetime = datetime.strptime(record["endTime"][:18], "%Y-%m-%dT%H:%M:%S")
    # Summary data
    row_summary = (device_id, record["efficiency"],
                    record["minutesAfterWakeup"], record["minutesAsleep"], record["minutesAwake"], record["minutesToFallAsleep"], record["timeInBed"],
                    d_is_main_sleep, sleep_record_type,
                    d_start_datetime, d_end_datetime,
                    0)
        
    records_summary.append(row_summary)
    
    return records_summary
    


def parseSleepData(sleep_data):
    if sleep_data.empty:
        return pd.DataFrame(columns=SLEEP_SUMMARY_COLUMNS)

    device_id = sleep_data["device_id"].iloc[0]
    records_summary = []
    # Parse JSON into individual records
    for multi_record in sleep_data.json_fitbit_column:
        sleep_record = json.loads(multi_record)
        if "sleep" in sleep_record:
            for record in sleep_record["sleep"]:
                # Whether the sleep episode is nap (0) or main sleep (1)
                d_is_main_sleep = 1 if record["isMainSleep"] else 0

                # For sleep API version 1
                if "awakeCount" in record:
                    records_summary = parseOneSleepRecord(record, device_id, d_is_main_sleep, records_summary, "classic")
                # For sleep API version 1.2
                else:
                    records_summary = parseOneSleepRecord(record, device_id, d_is_main_sleep, records_summary, record['type'])

    parsed_data = pd.DataFrame(data=records_summary, columns=SLEEP_SUMMARY_COLUMNS)

    return parsed_data


def main(json_raw, stream_parameters):
    parsed_data = parseSleepData(json_raw)

    if pd.api.types.is_datetime64_any_dtype( parsed_data['local_start_date_time']):
        parsed_data['local_start_date_time'] = parsed_data['local_start_date_time'].dt.strftime('%Y-%m-%d %H:%M:%S')
    if pd.api.types.is_datetime64_any_dtype( parsed_data['local_end_date_time']):
        parsed_data['local_end_date_time'] = parsed_data['local_end_date_time'].dt.strftime('%Y-%m-%d %H:%M:%S')

    return parsed_data
