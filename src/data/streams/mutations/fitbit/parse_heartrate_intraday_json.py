import json
import pandas as pd
from datetime import datetime


HR_INTRADAY_COLUMNS = ("device_id",
                        "heartrate", 
                        "heartrate_zone",
                        "local_date_time",
                        "timestamp")

def parseHeartrateZones(heartrate_data):
    # Get the range of heartrate zones: outofrange, fatburn, cardio, peak
    # refer to: https://help.fitbit.com/articles/en_US/Help_article/1565

    heartrate_fitbit_data = heartrate_data["activities-heart"][0]
    # API Version X: not sure the exact version
    if "heartRateZones" in heartrate_fitbit_data:
        heartrate_zones = heartrate_fitbit_data["heartRateZones"]
    # API VERSION Y: not sure the exact version
    elif "value" in heartrate_fitbit_data:
        heartrate_zones = heartrate_fitbit_data["value"]["heartRateZones"]
    else:
        raise ValueError("Heartrate zone are stored in an unknown format, this could mean Fitbit's heartrate API changed")
    
    heartrate_zones_range = {}
    for hrzone in heartrate_zones:
        heartrate_zones_range[hrzone["name"].lower().replace(" ", "")] = [hrzone["min"], hrzone["max"]]
    return heartrate_zones_range


def parseHeartrateIntradayData(records_intraday, dataset, device_id, curr_date, heartrate_zones_range):
    for data in dataset:
        d_time = datetime.strptime(data["time"], '%H:%M:%S').time()
        d_datetime = datetime.combine(curr_date, d_time).strftime("%Y-%m-%d %H:%M:%S")
        d_hr =  data["value"]

        # Get heartrate zone by range: min <= heartrate < max
        d_hrzone = None
        for hrzone, hrrange in heartrate_zones_range.items():
            if d_hr >= hrrange[0] and d_hr <= hrrange[1]:
                d_hrzone = hrzone
                break

        row_intraday = (device_id,
                        d_hr, d_hrzone,
                        d_datetime,
                        0)

        records_intraday.append(row_intraday)
    return records_intraday



def parseHeartrateData(heartrate_data):
    if heartrate_data.empty:
        return pd.DataFrame(columns=HR_INTRADAY_COLUMNS)

    device_id = heartrate_data["device_id"].iloc[0]
    records_intraday = []


    # Parse JSON into individual records
    for record in heartrate_data.json_fitbit_column:
        record = json.loads(record)  # Parse text into JSON
        if "activities-heart" in record:
            heartrate_zones_range = parseHeartrateZones(record)
            curr_date = datetime.strptime(record["activities-heart"][0]["dateTime"], "%Y-%m-%d")

            if "activities-heart-intraday" in record:
                dataset = record["activities-heart-intraday"]["dataset"]
                records_intraday = parseHeartrateIntradayData(records_intraday, dataset, device_id, curr_date, heartrate_zones_range)
    
    parsed_data = pd.DataFrame(data=records_intraday, columns=HR_INTRADAY_COLUMNS)
    return parsed_data
    


def main(json_raw, stream_parameters):
    parsed_data = parseHeartrateData(json_raw)
    return parsed_data
