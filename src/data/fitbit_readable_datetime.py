import pandas as pd
import pytz, json
from datetime import datetime


NIGHT = "night"
MORNING = "morning"
AFTERNOON = "afternoon"
EVENING = "evening"
HOUR2EPOCH = [NIGHT] * 6 + [MORNING] * 6 + [AFTERNOON] * 6 + [EVENING] * 6


HR_COLUMNS = ("device_id",
                "heartrate", 
                "local_date_time", "local_date", "local_month", "local_day",
                "local_day_of_week", "local_time", "local_hour", "local_minute", 
                "local_day_segment")

SLEEP_COLUMNS = ("device_id",
                    "sleep", # 1: "asleep", 2: "restless", or 3: "awake"
                    "local_date_time", "local_date", "local_month", "local_day",
                    "local_day_of_week", "local_time", "local_hour", "local_minute",
                    "local_day_segment")

STEPS_COLUMNS = ("device_id",
                    "steps",
                    "local_date_time", "local_date", "local_month", "local_day",
                    "local_day_of_week", "local_time", "local_hour", "local_minute",
                    "local_day_segment")

def drop_duplicates(data, local_timezone):
    """
    Data is pulled in intraday manner. Since data will be duplicated until the
    last record from that day, first sort by time, then drop all but
    the last record for each day. Drop duplicates based on aware timestamp.
    """
    local_date_col = data["timestamp"].apply(lambda ts: str(datetime.fromtimestamp(ts/1000, tz=local_timezone).date()))
    data = data.assign(local_date=local_date_col.values)
    data.sort_values(by="timestamp", ascending=True, inplace=True)
    data.drop_duplicates(subset="local_date", keep="last", inplace=True)

    return data

def parse_steps_data(steps_data):
    device_id = steps_data["device_id"].iloc[0]
    records = []
    # Parse JSON into individual records
    for record in steps_data.fitbit_data:
        record = json.loads(record)  # Parse text into JSON
        curr_date = datetime.strptime(
            record["activities-steps"][0]["dateTime"], "%Y-%m-%d")
        dataset = record["activities-steps-intraday"]["dataset"]
        for data in dataset:
            d_time = datetime.strptime(data["time"], '%H:%M:%S').time()
            d_datetime = datetime.combine(curr_date, d_time)

            row = (device_id,
                data["value"],
                d_datetime,
                d_datetime.date(),
                d_datetime.month,
                d_datetime.day,
                d_datetime.weekday(),
                d_datetime.time(),
                d_datetime.hour,
                d_datetime.minute,
                HOUR2EPOCH[d_datetime.hour])

            records.append(row)

    return pd.DataFrame(data=records, columns=STEPS_COLUMNS)

def parse_sleep_data(sleep_data):
    device_id = sleep_data["device_id"].iloc[0]
    records = []
    # Parse JSON into individual records
    for multi_record in sleep_data.fitbit_data:
        for record in json.loads(multi_record)["sleep"]:

            # Compute date when sleep episodes span two days
            start_date = datetime.strptime(record["startTime"][:10], "%Y-%m-%d")
            end_date = datetime.strptime(record["endTime"][:10], "%Y-%m-%d")
            flag = 1 if start_date == end_date else 0
            for data in record["minuteData"]:
                d_time = datetime.strptime(data["dateTime"], '%H:%M:%S').time()
                if not flag and not d_time.hour:
                    flag = 1
                curr_date = end_date if flag else start_date
                d_datetime = datetime.combine(curr_date, d_time)

                row = (device_id,
                    data["value"],
                    d_datetime,
                    d_datetime.date(),
                    d_datetime.month,
                    d_datetime.day,
                    d_datetime.weekday(),
                    d_datetime.time(),
                    d_datetime.hour,
                    d_datetime.minute,
                    HOUR2EPOCH[d_datetime.hour])

                records.append(row)

    return pd.DataFrame(data=records, columns=SLEEP_COLUMNS)

def parse_heartrate_data(heartrate_data):
    device_id = heartrate_data["device_id"].iloc[0]
    records = []
    # Parse JSON into individual records
    for record in heartrate_data.fitbit_data:
        record = json.loads(record)  # Parse text into JSON
        curr_date = datetime.strptime(record["activities-heart"][0]["dateTime"], "%Y-%m-%d")
        dataset = record["activities-heart-intraday"]["dataset"]
        for data in dataset:
            d_time = datetime.strptime(data["time"], '%H:%M:%S').time()
            d_datetime = datetime.combine(curr_date, d_time)

            row = (device_id,
                data["value"],
                d_datetime,
                d_datetime.date(),
                d_datetime.month,
                d_datetime.day,
                d_datetime.weekday(),
                d_datetime.time(),
                d_datetime.hour,
                d_datetime.minute,
                HOUR2EPOCH[d_datetime.hour])

            records.append(row)

    return pd.DataFrame(data=records, columns=HR_COLUMNS)


fitbit_data = pd.read_csv(snakemake.input[0])
local_timezone = pytz.timezone(snakemake.params["local_timezone"])
sensor = snakemake.params["fitbit_sensor"]

data = fitbit_data[fitbit_data["fitbit_data_type"] == sensor]
data = drop_duplicates(data, local_timezone)

if sensor == "heartrate":
    data_preprocesed = parse_heartrate_data(data)
elif sensor == "sleep":
    data_preprocesed = parse_sleep_data(data)
elif sensor == "steps":
    data_preprocesed = parse_steps_data(data)

data_preprocesed.to_csv(snakemake.output[0], index=False)
