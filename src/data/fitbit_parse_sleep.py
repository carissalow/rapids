import json, yaml
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import dateutil.parser

SLEEP_CODE2LEVEL = ["asleep", "restless", "awake"]


SLEEP_SUMMARY_COLUMNS_V1_2 = ("device_id", "efficiency",
                                "minutes_after_wakeup", "minutes_asleep", "minutes_awake", "minutes_to_fall_asleep", "minutes_in_bed",
                                "is_main_sleep", "type",
                                "local_start_date_time", "local_end_date_time",
                                "timestamp")
SLEEP_SUMMARY_COLUMNS_V1 = SLEEP_SUMMARY_COLUMNS_V1_2 + ("count_awake", "duration_awake", "count_awakenings", "count_restless", "duration_restless")

SLEEP_INTRADAY_COLUMNS = (# Extract "type_episode_id" field based on summary data: start from 0
                            "type_episode_id",
                            "duration",
                            # For "classic" type, original_level is one of {"awake", "restless", "asleep"}
                            # For "stages" type, original_level is one of {"wake", "deep", "light", "rem"}
                            "level",
                            # For "classic" type, unified_level is one of {0, 1} where 0: awake {"awake" + "restless"}, 1: asleep {"asleep"}
                            # For "stages" type, unified_level is one of {0, 1} where 0: awake {"wake"}, 1: asleep {"deep" + "light" + "rem"}
                            "unified_level",
                            # One of {0, 1} where 0: nap, 1: main sleep
                            "is_main_sleep", 
                            # One of {"classic", "stages"}
                            "type",
                            "local_date_time",
                            "start_timestamp",
                            "end_timestamp")


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
def parseOneRecordForV1(record, device_id, type_episode_id, d_is_main_sleep, records_summary, records_intraday, fitbit_data_type):

    sleep_record_type = "classic"

    d_start_datetime = datetime.strptime(record["startTime"][:18], "%Y-%m-%dT%H:%M:%S")
    d_end_datetime = datetime.strptime(record["endTime"][:18], "%Y-%m-%dT%H:%M:%S")

    # Summary data
    if fitbit_data_type == "summary":
        row_summary = (device_id, record["efficiency"],
                        record["minutesAfterWakeup"], record["minutesAsleep"], record["minutesAwake"], record["minutesToFallAsleep"], record["timeInBed"],
                        d_is_main_sleep, sleep_record_type,
                        d_start_datetime, d_end_datetime,
                        0,
                        record["awakeCount"], record["awakeDuration"], record["awakeningsCount"],
                        record["restlessCount"], record["restlessDuration"])
        
        records_summary.append(row_summary)

    # Intraday data
    if fitbit_data_type == "intraday":
        start_date = d_start_datetime.date()
        end_date = d_end_datetime.date()
        is_before_midnight = True
        curr_date = start_date
        for data in record["minuteData"]:
            # For overnight episodes, use end_date once we are over midnight
            d_time = datetime.strptime(data["dateTime"], '%H:%M:%S').time()
            if is_before_midnight and d_time.hour == 0:
                curr_date = end_date
            d_datetime = datetime.combine(curr_date, d_time)

            # API 1.2 stores original_level as strings, so we convert original_levels of API 1 to strings too
            # (1: "asleep", 2: "restless", 3: "awake")
            d_original_level = SLEEP_CODE2LEVEL[int(data["value"])-1]


            row_intraday = (type_episode_id, 60,
                            d_original_level, -1, d_is_main_sleep, sleep_record_type,
                            d_datetime, 0, 0)

            records_intraday.append(row_intraday)

    return records_summary, records_intraday

# Parse one record for sleep API version 1.2
def parseOneRecordForV12(record, device_id, type_episode_id, d_is_main_sleep, records_summary, records_intraday, fitbit_data_type):
    
    sleep_record_type = record['type']

    d_start_datetime = datetime.strptime(record["startTime"][:18], "%Y-%m-%dT%H:%M:%S")
    d_end_datetime = datetime.strptime(record["endTime"][:18], "%Y-%m-%dT%H:%M:%S")

    # Summary data
    if fitbit_data_type == "summary":
        row_summary = (device_id, record["efficiency"],
                        record["minutesAfterWakeup"], record["minutesAsleep"], record["minutesAwake"], record["minutesToFallAsleep"], record["timeInBed"],
                        d_is_main_sleep, sleep_record_type,
                        d_start_datetime, d_end_datetime,
                        0)
        
        records_summary.append(row_summary)
    
    # Intraday data
    if fitbit_data_type == "intraday":
        if sleep_record_type == "classic":
            for data in record["levels"]["data"]:
                d_datetime = dateutil.parser.parse(data["dateTime"])

                row_intraday = (type_episode_id, data["seconds"],
                    data["level"], -1, d_is_main_sleep, sleep_record_type,
                    d_datetime, 0, 0)
                records_intraday.append(row_intraday)
        else:
            # For sleep type "stages"
            for data in mergeLongAndShortData(record["levels"]):
                row_intraday = (type_episode_id, 30,
                    data[1], -1, d_is_main_sleep, sleep_record_type,
                    data[0], 0, 0)

                records_intraday.append(row_intraday)
    
    return records_summary, records_intraday

def parseSleepData(sleep_data, fitbit_data_type):
    SLEEP_SUMMARY_COLUMNS = SLEEP_SUMMARY_COLUMNS_V1_2
    if sleep_data.empty:
        if fitbit_data_type == "summary":
            return pd.DataFrame(columns=SLEEP_SUMMARY_COLUMNS)
        elif fitbit_data_type == "intraday":
            return pd.DataFrame(columns=SLEEP_INTRADAY_COLUMNS)
    device_id = sleep_data["device_id"].iloc[0]
    records_summary, records_intraday = [], []
    type_episode_id = 0
    # Parse JSON into individual records
    for multi_record in sleep_data.fitbit_data:
        for record in json.loads(multi_record)["sleep"]:
            # Whether the sleep episode is nap (0) or main sleep (1)
            d_is_main_sleep = 1 if record["isMainSleep"] else 0

            # For sleep API version 1
            if "awakeCount" in record:
                SLEEP_SUMMARY_COLUMNS = SLEEP_SUMMARY_COLUMNS_V1
                records_summary, records_intraday = parseOneRecordForV1(record, device_id, type_episode_id, d_is_main_sleep, records_summary, records_intraday, fitbit_data_type)
            # For sleep API version 1.2
            else:
                SLEEP_SUMMARY_COLUMNS = SLEEP_SUMMARY_COLUMNS_V1_2
                records_summary, records_intraday = parseOneRecordForV12(record, device_id, type_episode_id, d_is_main_sleep, records_summary, records_intraday, fitbit_data_type)
            
            type_episode_id = type_episode_id + 1

    if fitbit_data_type == "summary":
        parsed_data = pd.DataFrame(data=records_summary, columns=SLEEP_SUMMARY_COLUMNS)
    elif fitbit_data_type == "intraday":
        parsed_data = pd.DataFrame(data=records_intraday, columns=SLEEP_INTRADAY_COLUMNS)

    return parsed_data

def mergeSleepEpisodes(sleep_data, cols_for_groupby):
    sleep_episodes = pd.DataFrame(columns=["type_episode_id", "level_episode_id", "level", "unified_level", "is_main_sleep", "type", "start_timestamp", "end_timestamp"])
    if not sleep_data.empty:
        sleep_data = sleep_data.groupby(by=cols_for_groupby)
        sleep_episodes = sleep_data[["start_timestamp"]].first()
        sleep_episodes["end_timestamp"] = sleep_data["end_timestamp"].last()
    
        sleep_episodes.reset_index(inplace=True, drop=False)

    return sleep_episodes



timezone = snakemake.params["timezone"]
column_format = snakemake.params["column_format"]
fitbit_data_type = snakemake.params["fitbit_data_type"]
sleep_episode_timestamp = snakemake.params["sleep_episode_timestamp"]

with open(snakemake.input["participant_file"], "r", encoding="utf-8") as f:
    participant_file = yaml.safe_load(f)
local_start_date = pd.Timestamp(participant_file["FITBIT"]["START_DATE"])
local_end_date = pd.Timestamp(participant_file["FITBIT"]["END_DATE"]) + pd.DateOffset(1)

if column_format == "JSON":
    json_raw = pd.read_csv(snakemake.input["raw_data"])
    parsed_data = parseSleepData(json_raw, fitbit_data_type)
elif column_format == "PLAIN_TEXT":
    if fitbit_data_type == "summary":
        parsed_data = pd.read_csv(snakemake.input["raw_data"], parse_dates=["local_start_date_time", "local_end_date_time"], date_parser=lambda col: pd.to_datetime(col).tz_localize(None))
    elif fitbit_data_type == "intraday":
        parsed_data = pd.read_csv(snakemake.input["raw_data"], parse_dates=["local_date_time"], date_parser=lambda col: pd.to_datetime(col).tz_localize(None))
else:
    raise ValueError("column_format can only be one of ['JSON', 'PLAIN_TEXT'].")

# Drop duplicates
parsed_data.drop_duplicates(inplace=True)

if parsed_data.shape[0] > 0 and fitbit_data_type == "summary":
    if sleep_episode_timestamp != "start" and sleep_episode_timestamp != "end":
        raise ValueError("SLEEP_EPISODE_TIMESTAMP can only be one of ['start', 'end'].")
    # Column name to be considered as the event datetime
    datetime_column = "local_" + sleep_episode_timestamp + "_date_time"
    
    if not pd.isnull(local_start_date) and not pd.isnull(local_end_date):
        parsed_data = parsed_data.loc[(parsed_data[datetime_column] >= local_start_date) & (parsed_data[datetime_column] < local_end_date)]
    
    # Sort by "local_start_date_time" column
    parsed_data.sort_values(by="local_start_date_time", ascending=True, inplace=True)

    parsed_data["timestamp"] = parsed_data[datetime_column].dt.tz_localize(timezone, ambiguous=False, nonexistent="NaT").dropna().astype(np.int64) // 10**6
    parsed_data.dropna(subset=['timestamp'], inplace=True)
    parsed_data.drop(["local_start_date_time", "local_end_date_time"], axis = 1, inplace=True)

if parsed_data.shape[0] > 0 and fitbit_data_type == "intraday":
    if not pd.isnull(local_start_date) and not pd.isnull(local_end_date):
        parsed_data = parsed_data.loc[(parsed_data["local_date_time"] >= local_start_date) & (parsed_data["local_date_time"] < local_end_date)]
    
    # Sort by "local_date_time" column
    parsed_data.sort_values(by="local_date_time", ascending=True, inplace=True)

    parsed_data["start_timestamp"] = parsed_data["local_date_time"].dt.tz_localize(timezone, ambiguous=False, nonexistent="NaT").dropna().astype(np.int64) // 10**6
    parsed_data.dropna(subset=['start_timestamp'], inplace=True)
    parsed_data["end_timestamp"] = parsed_data["start_timestamp"] + ((parsed_data["duration"] - 1) * 1000) + 999
    parsed_data["unified_level"] = np.where(parsed_data["level"].isin(["awake", "restless", "wake"]), 0, 1)
    
    # Put consecutive rows with the same "level" field together and merge episodes
    parsed_data.insert(2, "level_episode_id", (parsed_data[["type_episode_id", "level"]] != parsed_data[["type_episode_id", "level"]].shift()).any(axis=1).cumsum())
    parsed_data = mergeSleepEpisodes(parsed_data, ["type_episode_id", "level_episode_id", "level", "unified_level", "is_main_sleep", "type"])


parsed_data.to_csv(snakemake.output[0], index=False)
