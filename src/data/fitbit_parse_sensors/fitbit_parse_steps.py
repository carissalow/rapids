import json
import pandas as pd
from datetime import datetime

STEPS_INTRADAY_COLUMNS = ("device_id",
                            "steps",
                            "local_date_time", "local_date", "local_month", "local_day",
                            "local_day_of_week", "local_time", "local_hour", "local_minute",
                            "local_day_segment")


def parseStepsData(steps_data, HOUR2EPOCH):
    if steps_data.empty:
        return pd.DataFrame(), pd.DataFrame(columns=STEPS_COLUMNS)
    device_id = steps_data["device_id"].iloc[0]
    records_intraday = []
    # Parse JSON into individual records
    for record in steps_data.fitbit_data:
        record = json.loads(record)  # Parse text into JSON
        curr_date = datetime.strptime(
            record["activities-steps"][0]["dateTime"], "%Y-%m-%d")
        dataset = record["activities-steps-intraday"]["dataset"]
        for data in dataset:
            d_time = datetime.strptime(data["time"], '%H:%M:%S').time()
            d_datetime = datetime.combine(curr_date, d_time)

            row_intraday = (device_id,
                data["value"],
                d_datetime, d_datetime.date(), d_datetime.month, d_datetime.day,
                d_datetime.weekday(), d_datetime.time(), d_datetime.hour, d_datetime.minute,
                HOUR2EPOCH[d_datetime.hour])

            records_intraday.append(row_intraday)

    return pd.DataFrame(), pd.DataFrame(data=records_intraday, columns=STEPS_INTRADAY_COLUMNS)
