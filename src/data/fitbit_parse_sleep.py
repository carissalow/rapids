import json
import pandas as pd
from datetime import datetime
import numpy as np
import dateutil.parser
from datetime import timedelta

SLEEP_CODE2LEVEL = ["asleep", "restless", "awake"]


SLEEP_SUMMARY_COLUMNS_V1_2 = ("device_id", "efficiency",
                                "minutes_after_wakeup", "minutes_asleep", "minutes_awake", "minutes_to_fall_asleep", "minutes_in_bed",
                                "is_main_sleep", "type",
                                "local_start_date_time", "local_end_date_time",
                                "start_timestamp", "end_timestamp")
SLEEP_SUMMARY_COLUMNS_V1 = SLEEP_SUMMARY_COLUMNS_V1_2 + ("count_awake", "duration_awake", "count_awakenings", "count_restless", "duration_restless")

SLEEP_INTRADAY_COLUMNS = ("device_id",
                            # For "classic" type, original_level is one of {"awake", "restless", "asleep"}
                            # For "stages" type, original_level is one of {"wake", "deep", "light", "rem"}
                            "level",
                            # For "classic" type, unified_level is one of {0, 1} where 0: awake {"awake" + "restless"}, 1: asleep {"asleep"}
                            # For "stages" type, unified_level is one of {0, 1} where 0: awake {"wake"}, 1: asleep {"deep" + "light" + "rem"}
                            "unified_level",
                            # one of {0, 1} where 0: nap, 1: main sleep
                            "is_main_sleep", 
                            # one of {"classic", "stages"}
                            "type",
                            "local_date_time",
                            "timestamp")

def mergeLongAndShortData(data_summary):
    longData = pd.DataFrame(columns=['dateTime', 'level', 'seconds'])
    shortData = pd.DataFrame(columns=['dateTime','level', 'seconds'])

    windowLength = 30

    for data in data_summary['data']:
        origEntry = data
        counter = 0
        numberOfSplits = origEntry['seconds']//windowLength
        for times in range(numberOfSplits):
            newRow = {'dateTime':dateutil.parser.parse(origEntry['dateTime'])+timedelta(seconds=counter*windowLength),'level':origEntry['level'],'seconds':windowLength}
            longData = longData.append(newRow, ignore_index = True)
            counter = counter + 1

    for data in data_summary['shortData']:
        origEntry = data
        counter = 0
        numberOfSplits = origEntry['seconds']//windowLength
        for times in range(numberOfSplits):
            newRow = {'dateTime':dateutil.parser.parse(origEntry['dateTime'])+timedelta(seconds=counter*windowLength),'level':origEntry['level'],'seconds':windowLength}
            shortData = shortData.append(newRow,ignore_index = True)
            counter = counter + 1
    longData.set_index('dateTime',inplace=True)
    shortData.set_index('dateTime',inplace=True)
    longData['level'] = np.where(longData.index.isin(shortData.index) == True,'wake',longData['level'])
    
    longData.reset_index(inplace=True)
    
    return longData.values.tolist()

def classicData1min(data_summary):
    dataList = list()
    for data in data_summary['data']:
        origEntry = data
        counter = 0
        timeDuration = 60
        numberOfSplits = origEntry['seconds']//timeDuration
        for times in range(numberOfSplits):
            newRow = {'dateTime':dateutil.parser.parse(origEntry['dateTime'])+timedelta(seconds=counter*timeDuration),'level':origEntry['level'],'seconds':timeDuration}
            dataList.append(newRow)
            counter = counter + 1
    # print(dataList)
    return dataList
# Parse one record for sleep API version 1
def parseOneRecordForV1(record, device_id, d_is_main_sleep, records_summary, records_intraday):

    # Summary data
    sleep_record_type = "classic"

    d_start_datetime = datetime.strptime(record["startTime"][:18], "%Y-%m-%dT%H:%M:%S")
    d_end_datetime = datetime.strptime(record["endTime"][:18], "%Y-%m-%dT%H:%M:%S")

    row_summary = (device_id, record["efficiency"],
                    record["minutesAfterWakeup"], record["minutesAsleep"], record["minutesAwake"], record["minutesToFallAsleep"], record["timeInBed"],
                    d_is_main_sleep, sleep_record_type,
                    d_start_datetime, d_end_datetime,
                    d_start_datetime.date(), d_end_datetime.date(),
                    0,0,
                    record["awakeCount"], record["awakeDuration"], record["awakeningsCount"],
                    record["restlessCount"], record["restlessDuration"])
    
    records_summary.append(row_summary)

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
        d_datetime = datetime.combine(curr_date, d_time)

        # API 1.2 stores original_level as strings, so we convert original_levels of API 1 to strings too
        # (1: "asleep", 2: "restless", 3: "awake")
        d_original_level = SLEEP_CODE2LEVEL[int(data["value"])-1]


        row_intraday = (device_id,
                        d_original_level, -1, d_is_main_sleep, sleep_record_type,
                        d_datetime, 0)

        records_intraday.append(row_intraday)

    return records_summary, records_intraday

# Parse one record for sleep API version 1.2
def parseOneRecordForV12(record, device_id, d_is_main_sleep, records_summary, records_intraday):
    
    # Summary data
    sleep_record_type = record['type']

    d_start_datetime = datetime.strptime(record["startTime"][:18], "%Y-%m-%dT%H:%M:%S")
    d_end_datetime = datetime.strptime(record["endTime"][:18], "%Y-%m-%dT%H:%M:%S")

    row_summary = (device_id, record["efficiency"],
                    record["minutesAfterWakeup"], record["minutesAsleep"], record["minutesAwake"], record["minutesToFallAsleep"], record["timeInBed"],
                    d_is_main_sleep, sleep_record_type,
                    d_start_datetime, d_end_datetime,
                    0,0)
    
    records_summary.append(row_summary)
    if sleep_record_type == 'classic':
            # Intraday data
            start_date = d_start_datetime.date()
            end_date = d_end_datetime.date()
            is_before_midnight = True
            curr_date = start_date
            data_summary = record['levels']
            dataSplitted = classicData1min(data_summary)  ##Calling the function to split the data in regular 60 seconds interval
            for data in dataSplitted:
                # For overnight episodes, use end_date once we are over midnight
                d_time = data["dateTime"].time()
                if is_before_midnight and d_time.hour == 0:
                    curr_date = end_date
                d_datetime = datetime.combine(curr_date, d_time)

                d_original_level = data["level"]

                row_intraday = (device_id,
                    d_original_level, -1, d_is_main_sleep, sleep_record_type,
                    d_datetime, 0)
                records_intraday.append(row_intraday)
    else:
            ## for sleep type "stages"
            start_date = d_start_datetime.date()
            end_date = d_end_datetime.date()
            is_before_midnight = True
            curr_date = start_date
            data_summary = record['levels']
            dataList = mergeLongAndShortData(data_summary)
            for data in dataList:

                d_time = data[0].time()
                if is_before_midnight and d_time.hour == 0:
                    curr_date = end_date
                d_datetime = datetime.combine(curr_date, d_time)

                d_original_level = data[1]

                row_intraday = (device_id,
                    d_original_level, -1, d_is_main_sleep, sleep_record_type,
                    d_datetime, 0)

                records_intraday.append(row_intraday)
    
    return records_summary, records_intraday
    


def parseSleepData(sleep_data):
    SLEEP_SUMMARY_COLUMNS = SLEEP_SUMMARY_COLUMNS_V1_2
    if sleep_data.empty:
        return pd.DataFrame(columns=SLEEP_SUMMARY_COLUMNS), pd.DataFrame(columns=SLEEP_INTRADAY_COLUMNS)
    device_id = sleep_data["device_id"].iloc[0]
    records_summary, records_intraday = [], []
    # Parse JSON into individual records
    for multi_record in sleep_data.fitbit_data:
        for record in json.loads(multi_record)["sleep"]:
            # Whether the sleep episode is nap (0) or main sleep (1)
            d_is_main_sleep = 1 if record["isMainSleep"] else 0

            # For sleep API version 1
            if "awakeCount" in record:
                SLEEP_SUMMARY_COLUMNS = SLEEP_SUMMARY_COLUMNS_V1
                records_summary, records_intraday = parseOneRecordForV1(record, device_id, d_is_main_sleep, records_summary, records_intraday)
            # For sleep API version 1.2
            else:
                SLEEP_SUMMARY_COLUMNS = SLEEP_SUMMARY_COLUMNS_V1_2
                records_summary, records_intraday = parseOneRecordForV12(record, device_id, d_is_main_sleep, records_summary, records_intraday)

    return pd.DataFrame(data=records_summary, columns=SLEEP_SUMMARY_COLUMNS), pd.DataFrame(data=records_intraday, columns=SLEEP_INTRADAY_COLUMNS)

table_format = snakemake.params["table_format"]
timezone = snakemake.params["timezone"]

if table_format == "JSON":
    json_raw = pd.read_csv(snakemake.input[0])
    summary, intraday = parseSleepData(json_raw)
elif table_format == "CSV":
    summary = pd.read_csv(snakemake.input[0])
    intraday = pd.read_csv(snakemake.input[1])

if summary.shape[0] > 0:
    summary["start_timestamp"] = summary["local_start_date_time"].dt.tz_localize(timezone).astype(np.int64) // 10**6
    summary["end_timestamp"] = summary["local_end_date_time"].dt.tz_localize(timezone).astype(np.int64) // 10**6
if intraday.shape[0] > 0:
    intraday["timestamp"] = intraday["local_date_time"].dt.tz_localize(timezone).astype(np.int64) // 10**6

# Unifying level
intraday["unified_level"] = np.where(intraday["level"].isin(["awake", "wake", "restless"]), 0, 1)

summary.to_csv(snakemake.output["summary_data"], index=False)
intraday.to_csv(snakemake.output["intraday_data"], index=False)