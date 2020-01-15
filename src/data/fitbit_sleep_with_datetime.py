import pandas as pd
import pytz, json
from datetime import datetime



NIGHT = "night"
MORNING = "morning"
AFTERNOON = "afternoon"
EVENING = "evening"
HOUR2EPOCH = [NIGHT] * 6 + [MORNING] * 6 + [AFTERNOON] * 6 + [EVENING] * 6


SLEEP_COLUMNS = ("device_id",
                 "sleep", # 1: "asleep", 2: "restless", or 3: "awake"
                 "local_date_time",
                 "local_date",
                 "local_month",
                 "local_day",
                 "local_day_of_week",
                 "local_time",
                 "local_hour",
                 "local_minute",
                 "local_day_segment")

fitbit_data = pd.read_csv(snakemake.input[0])
sleep_data = fitbit_data[fitbit_data["fitbit_data_type"] == "sleep"]

local_timezone = pytz.timezone(snakemake.params["local_timezone"])


"""
Data is pulled in intraday manner. Since data will be duplicated until the
last record from that day, first sort by time, then drop all but
the last record for each day. Drop duplicates based on aware timestamp.
"""
local_date_col = sleep_data["timestamp"].apply(lambda ts: str(datetime.fromtimestamp(ts/1000, tz=local_timezone).date()))
sleep_data = sleep_data.assign(local_date=local_date_col.values)
sleep_data.sort_values(by="timestamp", ascending=True, inplace=True)
sleep_data.drop_duplicates(subset="local_date", keep="last", inplace=True)

device_id = sleep_data["device_id"].iloc[0]
records = []
# Parse JSON into individual records
for multi_record in sleep_data.fitbit_data:
    for record in json.loads(multi_record)["sleep"]:
        start_date = datetime.strptime(record["startTime"][:10], "%Y-%m-%d")
        end_date = datetime.strptime(record["endTime"][:10], "%Y-%m-%d")
        flag = 1 if start_date == end_date else 0
        for data in record["minuteData"]:
            d_time = datetime.strptime(data["dateTime"], '%H:%M:%S').time()
            if not flag and not d_time.hour:
                flag = 1
            curr_date = end_date if flag else start_date
            d_datetime = datetime.combine(curr_date, d_time)

            # Create tuple of parsed data
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

            # Append the data to a list
            records.append(row)

# Create a new DataFrame from the list of tuples.
sleep_preprocessed = pd.DataFrame(data=records, columns=SLEEP_COLUMNS)

sleep_preprocessed.to_csv(snakemake.output[0], index=False)
