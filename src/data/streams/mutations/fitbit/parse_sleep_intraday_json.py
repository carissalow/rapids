import json
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import dateutil.parser

SLEEP_CODE_TO_LEVEL = {
    1: "asleep", 
    2: "restless", 
    3: "awake"
}

SLEEP_INTRADAY_COLUMNS = [
    "device_id",
    "type_episode_id",
    "duration",
    "level",            # for "classic" type, one of {"awake", "restless", "asleep"}; for "stages" type, one of {"wake", "deep", "light", "rem"}
    "is_main_sleep",    # one of {0, 1} where 0: nap, 1: main sleep
    "type",             # one of {"classic", "stages"}
    "local_date_time",
    "timestamp"
]

def resampleLevelsData(levels_data, element):
    if not element in ["data", "shortData"]:
        raise ValueError(f"Unrecognized levels element: {element}")

    window_length = 30
    resampled_data = []

    for data in levels_data[element]:
        for i in range(data["seconds"] // window_length):
            row = {
                "dateTime": dateutil.parser.parse(data["dateTime"]) + timedelta(seconds=i*window_length), 
                "level": data["level"]
            }
            resampled_data.append(row)

    resampled_data = pd.DataFrame(resampled_data, columns=["dateTime", "level"]).set_index("dateTime")
    return resampled_data


def mergeLongAndShortLevelsData(levels_data):
    long_data = resampleLevelsData(levels_data, "data")

    if "shortData" in levels_data:
        short_data = resampleLevelsData(levels_data, "shortData")
        long_data["level"] = np.where(long_data.index.isin(short_data.index), "wake", long_data["level"])

    long_data.reset_index(inplace=True)
    return long_data.values.tolist()


def parseOneRecordForV1(record, is_main_sleep, records_intraday, type_episode_id):
    sleep_record_type = "classic"
    duration = 60
    timestamp = 0

    d_start_datetime = datetime.strptime(record["startTime"][:18], "%Y-%m-%dT%H:%M:%S")
    d_end_datetime = datetime.strptime(record["endTime"][:18], "%Y-%m-%dT%H:%M:%S")

    start_date = d_start_datetime.date()
    end_date = d_end_datetime.date()
    is_before_midnight = True
    curr_date = start_date

    for data in record["minuteData"]:
        # For overnight episodes, use end_date once we are over midnight
        d_time = datetime.strptime(data["dateTime"], "%H:%M:%S").time()
        if is_before_midnight and d_time.hour == 0:
            curr_date = end_date
        
        d_datetime = datetime.combine(curr_date, d_time).strftime("%Y-%m-%d %H:%M:%S")
        d_original_level = SLEEP_CODE_TO_LEVEL[int(data["value"])]

        row_intraday = {
            "type_episode_id": type_episode_id,
            "duration": duration,
            "level": d_original_level,           
            "is_main_sleep": is_main_sleep,   
            "type": sleep_record_type,             
            "local_date_time": d_datetime,
            "timestamp": timestamp
        }
        records_intraday.append(row_intraday)

    return records_intraday


def parseOneRecordForV12(record, is_main_sleep, records_intraday, type_episode_id):
    stages_duration = 30
    timestamp = 0

    sleep_record_type = record["type"]

    if sleep_record_type == "classic":
        for data in record["levels"]["data"]:
            d_datetime = data["dateTime"][:19].replace("T", " ")
            row_intraday = {
                "type_episode_id": type_episode_id,
                "duration": data["seconds"],
                "level": data["level"],           
                "is_main_sleep": is_main_sleep,   
                "type": sleep_record_type,             
                "local_date_time": d_datetime,
                "timestamp": timestamp
            }
            records_intraday.append(row_intraday)
    
    elif sleep_record_type == "stages":
        for data in mergeLongAndShortLevelsData(record["levels"]):
            d_datetime = data[0].strftime("%Y-%m-%d %H:%M:%S")
            row_intraday = {
                "type_episode_id": type_episode_id,
                "duration": stages_duration,
                "level": data[1],           
                "is_main_sleep": is_main_sleep,   
                "type": sleep_record_type,             
                "local_date_time": d_datetime,
                "timestamp": timestamp
            }
            records_intraday.append(row_intraday)
    
    else:
        raise ValueError(f"Unrecognized sleep record type: {sleep_record_type}")

    return records_intraday
    

def parseSleepData(sleep_data):
    if sleep_data.empty:
        return pd.DataFrame(columns=SLEEP_INTRADAY_COLUMNS)

    device_id = sleep_data["device_id"].iloc[0]
    records_intraday = []
    type_episode_id = 0

    for multi_record in sleep_data.json_fitbit_column:
        sleep_record = json.loads(multi_record)
        if "sleep" in sleep_record:
            for record in sleep_record["sleep"]:
                is_main_sleep = 1 if record["isMainSleep"] else 0

                if "awakeCount" in record: 
                    records_intraday = parseOneRecordForV1(record, is_main_sleep, records_intraday, type_episode_id)
                else:
                    records_intraday = parseOneRecordForV12(record, is_main_sleep, records_intraday, type_episode_id)
                
                type_episode_id += 1

    parsed_data = pd.DataFrame(data=records_intraday)
    parsed_data.insert(0, column="device_id", value=device_id)
    
    return parsed_data


def main(json_raw, stream_parameters):
    parsed_data = parseSleepData(json_raw)
    return parsed_data
