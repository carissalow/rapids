import pandas as pd
import pytz, json
from datetime import datetime


NIGHT = "night"
MORNING = "morning"
AFTERNOON = "afternoon"
EVENING = "evening"
HOUR2EPOCH = [NIGHT] * 6 + [MORNING] * 6 + [AFTERNOON] * 6 + [EVENING] * 6


STEPS_COLUMNS = ("device_id",
                 "steps",
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
steps_data = fitbit_data[fitbit_data["fitbit_data_type"] == "steps"]

local_timezone = pytz.timezone(snakemake.params["local_timezone"])


"""
Data is pulled in intraday manner. Since data will be duplicated until the
last record from that day, first sort by time, then drop all but
the last record for each day. Drop duplicates based on aware timestamp.
"""
local_date_col = steps_data["timestamp"].apply(lambda ts: str(datetime.fromtimestamp(ts/1000, tz=local_timezone).date()))
steps_data = steps_data.assign(local_date=local_date_col.values)
steps_data.sort_values(by="timestamp", ascending=True, inplace=True)
steps_data.drop_duplicates(subset="local_date", keep="last", inplace=True)

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
steps_preprocessed = pd.DataFrame(data=records, columns=STEPS_COLUMNS)

steps_preprocessed.to_csv(snakemake.output[0], index=False)
