import json
import pandas as pd


HR_SUMMARY_COLUMNS = ("device_id",
                        "local_date_time",
                        "timestamp",
                        "heartrate_daily_restinghr",
                        "heartrate_daily_caloriesoutofrange",
                        "heartrate_daily_caloriesfatburn",
                        "heartrate_daily_caloriescardio",
                        "heartrate_daily_caloriespeak")


def parseHeartrateSummaryData(record_summary, device_id, curr_date):
    # API Version X: not sure the exact version
    if "heartRateZones" in record_summary:
        heartrate_zones = record_summary["heartRateZones"]
        d_resting_heartrate = record_summary["value"] if "value" in record_summary else None
    # API VERSION Y: not sure the exact version
    elif "value" in record_summary:
        heartrate_zones = record_summary["value"]["heartRateZones"]
        d_resting_heartrate = record_summary["value"]["restingHeartRate"] if "restingHeartRate" in record_summary["value"] else None
    else:
        ValueError("Heartrate zone are stored in an unkown format, this could mean Fitbit's heartrate API changed")
    
    if "caloriesOut" in heartrate_zones[0]:
        d_calories_outofrange = heartrate_zones[0]["caloriesOut"]
        d_calories_fatburn = heartrate_zones[1]["caloriesOut"]
        d_calories_cardio = heartrate_zones[2]["caloriesOut"]
        d_calories_peak = heartrate_zones[3]["caloriesOut"]
    else:
        d_calories_outofrange, d_calories_fatburn, d_calories_cardio, d_calories_peak = None, None, None, None
    
    row_summary = (device_id,
                    curr_date,
                    0,
                    d_resting_heartrate,
                    d_calories_outofrange,
                    d_calories_fatburn,
                    d_calories_cardio,
                    d_calories_peak)
    return row_summary

def parseHeartrateData(heartrate_data):
    if heartrate_data.empty:
        return pd.DataFrame(columns=HR_SUMMARY_COLUMNS)

    device_id = heartrate_data["device_id"].iloc[0]
    records_summary = []


    # Parse JSON into individual records
    for record in heartrate_data.json_fitbit_column:
        record = json.loads(record)  # Parse text into JSON
        if "activities-heart" in record:
            curr_date = record["activities-heart"][0]["dateTime"] + " 00:00:00"

            record_summary = record["activities-heart"][0]
            row_summary = parseHeartrateSummaryData(record_summary, device_id, curr_date)
            records_summary.append(row_summary)
    parsed_data = pd.DataFrame(data=records_summary, columns=HR_SUMMARY_COLUMNS)

    return parsed_data
    

def main(json_raw, stream_parameters):
    parsed_data = parseHeartrateData(json_raw)
    return parsed_data
