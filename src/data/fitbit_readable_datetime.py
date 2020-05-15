import pandas as pd
import pytz, json
from datetime import datetime
from fitbit_parse_sensors.fitbit_parse_heartrate import parseHeartrateData
from fitbit_parse_sensors.fitbit_parse_sleep import parseSleepData
from fitbit_parse_sensors.fitbit_parse_steps import parseStepsData
from fitbit_parse_sensors.fitbit_parse_calories import parseCaloriesData


NIGHT = "night"
MORNING = "morning"
AFTERNOON = "afternoon"
EVENING = "evening"
HOUR2EPOCH = [NIGHT] * 6 + [MORNING] * 6 + [AFTERNOON] * 6 + [EVENING] * 6


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



fitbit_data = pd.read_csv(snakemake.input[0])
local_timezone = pytz.timezone(snakemake.params["local_timezone"])
sensor = snakemake.params["fitbit_sensor"]

data = fitbit_data[fitbit_data["fitbit_data_type"] == sensor]
data = drop_duplicates(data, local_timezone)

if sensor == "heartrate":
    summary_data, intraday_data = parseHeartrateData(data, HOUR2EPOCH)
elif sensor == "sleep":
    summary_data, intraday_data = parseSleepData(data, HOUR2EPOCH)
elif sensor == "steps":
    summary_data, intraday_data = parseStepsData(data, HOUR2EPOCH)
elif sensor == "calories":
    summary_data, intraday_data = parseCaloriesData(data, HOUR2EPOCH)
else:
    raise ValueError("Please check the FITBIT_SENSORS list in config.yaml file.")

# Summary data will be empty for steps and calories as it is not provided by Fitbit's API
summary_data.to_csv(snakemake.output["summary_data"], index=False)
intraday_data.to_csv(snakemake.output["intraday_data"], index=False)
