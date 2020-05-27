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
                                "local_start_date", "local_end_date",
                                "local_start_day_segment", "local_end_day_segment")
SLEEP_SUMMARY_COLUMNS_V1 = SLEEP_SUMMARY_COLUMNS_V1_2 + ("count_awake", "duration_awake", "count_awakenings", "count_restless", "duration_restless")

SLEEP_INTRADAY_COLUMNS = ("device_id",
                            # For "classic" type, original_level is one of {"awake", "restless", "asleep"}
                            # For "stages" type, original_level is one of {"wake", "deep", "light", "rem"}
                            "original_level",
                            # For "classic" type, unified_level is one of {0, 1} where 0: awake {"awake" + "restless"}, 1: asleep {"asleep"}
                            # For "stages" type, unified_level is one of {0, 1} where 0: awake {"wake"}, 1: asleep {"deep" + "light" + "rem"}
                            "unified_level",
                            # one of {0, 1} where 0: nap, 1: main sleep
                            "is_main_sleep", 
                            # one of {"classic", "stages"}
                            "type",
                            "local_date_time", "local_date", "local_month", "local_day",
                            "local_day_of_week", "local_time", "local_hour", "local_minute",
                            "local_day_segment")

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

# Parse one record for sleep API version 1
def parseOneRecordForV1(record, device_id, d_is_main_sleep, records_summary, records_intraday, HOUR2EPOCH):

    # Summary data
    sleep_record_type = "classic"

    d_start_datetime = datetime.strptime(record["startTime"][:18], "%Y-%m-%dT%H:%M:%S")
    d_end_datetime = datetime.strptime(record["endTime"][:18], "%Y-%m-%dT%H:%M:%S")

    row_summary = (device_id, record["efficiency"],
                    record["minutesAfterWakeup"], record["minutesAsleep"], record["minutesAwake"], record["minutesToFallAsleep"], record["timeInBed"],
                    d_is_main_sleep, sleep_record_type,
                    d_start_datetime, d_end_datetime,
                    d_start_datetime.date(), d_end_datetime.date(),
                    HOUR2EPOCH[d_start_datetime.hour], HOUR2EPOCH[d_end_datetime.hour],
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

        # unified_level summarises original_level (we came up with this classification)
        # 0 is awake, 1 is asleep
        # {"awake" + "restless"} are set to 0 and {"asleep"} is set to 1
        d_unified_level = 0 if d_original_level == "awake" or d_original_level == "restless" else 1

        row_intraday = (device_id,
                        d_original_level, d_unified_level, d_is_main_sleep, sleep_record_type,
                        d_datetime, d_datetime.date(), d_datetime.month, d_datetime.day,
                        d_datetime.weekday(), d_datetime.time(), d_datetime.hour, d_datetime.minute,
                        HOUR2EPOCH[d_datetime.hour])

        records_intraday.append(row_intraday)

    return records_summary, records_intraday

# Parse one record for sleep API version 1.2
def parseOneRecordForV12(record, device_id, d_is_main_sleep, records_summary, records_intraday, HOUR2EPOCH):
    
    # Summary data
    sleep_record_type = record['type']

    d_start_datetime = datetime.strptime(record["startTime"][:18], "%Y-%m-%dT%H:%M:%S")
    d_end_datetime = datetime.strptime(record["endTime"][:18], "%Y-%m-%dT%H:%M:%S")

    row_summary = (device_id, record["efficiency"],
                    record["minutesAfterWakeup"], record["minutesAsleep"], record["minutesAwake"], record["minutesToFallAsleep"], record["timeInBed"],
                    d_is_main_sleep, sleep_record_type,
                    d_start_datetime, d_end_datetime,
                    d_start_datetime.date(), d_end_datetime.date(),
                    HOUR2EPOCH[d_start_datetime.hour], HOUR2EPOCH[d_end_datetime.hour])
    
    records_summary.append(row_summary)
    if sleep_record_type == 'classic':
            # Intraday data
            start_date = d_start_datetime.date()
            end_date = d_end_datetime.date()
            is_before_midnight = True
            curr_date = start_date
            data_summary = record['levels']
            for data in data_summary['data']:
                # For overnight episodes, use end_date once we are over midnight
                d_time = dateutil.parser.parse(data["dateTime"]).time()
                if is_before_midnight and d_time.hour == 0:
                    curr_date = end_date
                d_datetime = datetime.combine(curr_date, d_time)

                d_original_level = data["level"]

                d_unified_level = 0 if d_original_level == "awake" or d_original_level == "restless" else 1

                row_intraday = (device_id,
                    d_original_level, d_unified_level, d_is_main_sleep, sleep_record_type,
                    d_datetime, d_datetime.date(), d_datetime.month, d_datetime.day,
                    d_datetime.weekday(), d_datetime.time(), d_datetime.hour, d_datetime.minute,
                    HOUR2EPOCH[d_datetime.hour])
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

                d_unified_level = 1 if d_original_level == "deep" or d_original_level == "light" or d_original_level == "rem" else 0

                row_intraday = (device_id,
                    d_original_level, d_unified_level, d_is_main_sleep, sleep_record_type,
                    d_datetime, d_datetime.date(), d_datetime.month, d_datetime.day,
                    d_datetime.weekday(), d_datetime.time(), d_datetime.hour, d_datetime.minute,
                    HOUR2EPOCH[d_datetime.hour])

                records_intraday.append(row_intraday)
    
    return records_summary, records_intraday
    


def parseSleepData(sleep_data, HOUR2EPOCH):
    if sleep_data.empty:
        return pd.DataFrame(columns=SLEEP_SUMMARY_COLUMNS_V1), pd.DataFrame(columns=SLEEP_INTRADAY_COLUMNS)
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
                records_summary, records_intraday = parseOneRecordForV1(record, device_id, d_is_main_sleep, records_summary, records_intraday, HOUR2EPOCH)
            # For sleep API version 1.2
            else:
                SLEEP_SUMMARY_COLUMNS = SLEEP_SUMMARY_COLUMNS_V1_2
                records_summary, records_intraday = parseOneRecordForV12(record, device_id, d_is_main_sleep, records_summary, records_intraday, HOUR2EPOCH)

    return pd.DataFrame(data=records_summary, columns=SLEEP_SUMMARY_COLUMNS), pd.DataFrame(data=records_intraday, columns=SLEEP_INTRADAY_COLUMNS)
