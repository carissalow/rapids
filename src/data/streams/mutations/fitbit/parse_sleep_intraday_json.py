import json
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import dateutil.parser

SLEEP_CODE2LEVEL = ["asleep", "restless", "awake"]

SLEEP_INTRADAY_COLUMNS = ("device_id",
                            "type_episode_id",
                            "duration",
                            # For "classic" type, original_level is one of {"awake", "restless", "asleep"}
                            # For "stages" type, original_level is one of {"wake", "deep", "light", "rem"}
                            "level",
                            # one of {0, 1} where 0: nap, 1: main sleep
                            "is_main_sleep", 
                            # one of {"classic", "stages"}
                            "type",
                            "local_date_time",
                            "timestamp")

def mergeLongAndShortData(data_intraday):
    long_data = pd.DataFrame(columns=["dateTime", "level"])
    short_data = pd.DataFrame(columns=["dateTime", "level"])

    window_length = 30

    for data in data_intraday["data"]:
        counter = 0
        for times in range(data["seconds"] // window_length):
            row = {"dateTime": dateutil.parser.parse(data["dateTime"])+timedelta(seconds=counter*window_length), "level": data["level"]}
            long_data = long_data.append(row, ignore_index = True)
            counter = counter + 1

    for data in data_intraday["shortData"]:
        counter = 0
        for times in range(data["seconds"] // window_length):
            row = {"dateTime": dateutil.parser.parse(data["dateTime"])+timedelta(seconds=counter*window_length), "level": data["level"]}
            short_data = short_data.append(row, ignore_index = True)
            counter = counter + 1
    long_data.set_index("dateTime",inplace=True)
    short_data.set_index("dateTime",inplace=True)
    long_data["level"] = np.where(long_data.index.isin(short_data.index) == True, "wake", long_data["level"])
    
    long_data.reset_index(inplace=True)
    
    return long_data.values.tolist()

# Parse one record for sleep API version 1
def parseOneRecordForV1(record, device_id, d_is_main_sleep, records_intraday, type_episode_id):

    sleep_record_type = "classic"

    d_start_datetime = datetime.strptime(record["startTime"][:18], "%Y-%m-%dT%H:%M:%S")
    d_end_datetime = datetime.strptime(record["endTime"][:18], "%Y-%m-%dT%H:%M:%S")

    # Intraday data
    start_date = d_start_datetime.date()
    end_date = d_end_datetime.date()
    is_before_midnight = True
    curr_date = start_date
    for data in record["minuteData"]:
        # For overnight episodes, use end_date once we are over midnight
        d_time = datetime.strptime(data["dateTime"], '%H:%M:%S').time()
        if is_before_midnight and d_time.hour == 0:
            curr_date = end_date
        d_datetime = datetime.combine(curr_date, d_time).strftime("%Y-%m-%d %H:%M:%S")

        # API 1.2 stores original_level as strings, so we convert original_levels of API 1 to strings too
        # (1: "asleep", 2: "restless", 3: "awake")
        d_original_level = SLEEP_CODE2LEVEL[int(data["value"])-1]


        row_intraday = (device_id, type_episode_id, 60,
                        d_original_level, d_is_main_sleep, sleep_record_type,
                        d_datetime, 0)

        records_intraday.append(row_intraday)

    return records_intraday

# Parse one record for sleep API version 1.2
def parseOneRecordForV12(record, device_id, d_is_main_sleep, records_intraday, type_episode_id):
    
    sleep_record_type = record['type']

    if sleep_record_type == "classic":
        for data in record["levels"]["data"]:
            d_datetime = data["dateTime"][:19].replace("T", " ")

            row_intraday = (device_id, type_episode_id, data["seconds"],
                data["level"], d_is_main_sleep, sleep_record_type,
                d_datetime, 0)
            records_intraday.append(row_intraday)
    else:
        # For sleep type "stages"
        for data in mergeLongAndShortData(record["levels"]):
            d_datetime = data[0].strftime("%Y-%m-%d %H:%M:%S")
            row_intraday = (device_id, type_episode_id, 30,
                data[1], d_is_main_sleep, sleep_record_type,
                d_datetime, 0)

            records_intraday.append(row_intraday)
    
    return records_intraday
    


def parseSleepData(sleep_data):
    if sleep_data.empty:
        return pd.DataFrame(columns=SLEEP_INTRADAY_COLUMNS)
    device_id = sleep_data["device_id"].iloc[0]
    records_intraday = []
    type_episode_id = 0
    # Parse JSON into individual records
    for multi_record in sleep_data.json_fitbit_column:
        sleep_record = json.loads(multi_record)
        if "sleep" in sleep_record:
            for record in json.loads(multi_record)["sleep"]:
                # Whether the sleep episode is nap (0) or main sleep (1)
                d_is_main_sleep = 1 if record["isMainSleep"] else 0

                # For sleep API version 1
                if "awakeCount" in record:
                    records_intraday = parseOneRecordForV1(record, device_id, d_is_main_sleep, records_intraday, type_episode_id)
                # For sleep API version 1.2
                else:
                    records_intraday = parseOneRecordForV12(record, device_id, d_is_main_sleep, records_intraday, type_episode_id)
                
                type_episode_id = type_episode_id + 1

    parsed_data = pd.DataFrame(data=records_intraday, columns=SLEEP_INTRADAY_COLUMNS)

    return parsed_data



def main(json_raw, stream_parameters):
    parsed_data = parseSleepData(json_raw)
    return parsed_data
