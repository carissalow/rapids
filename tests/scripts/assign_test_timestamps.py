import numpy as np
import pandas as pd
import argparse, glob

DayOfWeek2Date = {"Fri": ["2020-03-06", "2020-10-30"],
                  "Sat": ["2020-03-07", "2020-10-31"],
                  "Sun": ["2020-03-08", "2020-11-01"],
                  "Mon": ["2020-03-09", "2020-11-02"]}

def assign_test_timestamps(file_path):

    columns_to_delete = ["test_time", "day_of_week", "time"]

    data = pd.read_csv(file_path)
    data[["day_of_week", "time"]] = data["test_time"].str.split(pat=" ", n=1, expand=True)

    data_with_timestamps = pd.DataFrame()

    # 0 is for March and 1 is for Nov
    for i in [0, 1]:
        data["local_date_time"] = pd.to_datetime(data.apply(lambda row: DayOfWeek2Date[row["day_of_week"]][i] + " " + row["time"], axis=1))
        data_with_timestamps = pd.concat([data_with_timestamps, data], axis=0)
    
    if "fitbit" in file_path:
        data_with_timestamps.insert(0, "timestamp", 0)
        data_with_timestamps["local_date_time"] = data_with_timestamps["local_date_time"].dt.strftime('%Y-%m-%d %H:%M:%S')
    else:
        # Convert local_date_time with timezone to timestamp
        data_with_timestamps.insert(0, "timestamp", data_with_timestamps["local_date_time"].dt.tz_localize(tz="America/New_York").astype(np.int64) // 10**6)
        columns_to_delete.append("local_date_time")
    
    # Discard useless columns
    for col in columns_to_delete:
        del data_with_timestamps[col]
    
    return data_with_timestamps



parser = argparse.ArgumentParser()
parser.add_argument("-f", "--files", nargs="+", help="Assign timestamps to the selected files, it could be a single file name or multiple file names separated by whitespace(s) (e.g. phone_battery_raw.csv)")
parser.add_argument("-a", "--all", action="store_true", help="Assign timestamps to all files under the tests/data/manual/aware_csv folder")

args = parser.parse_args()
if args.all:
    for file_path in glob.glob("tests/data/manual/aware_csv/*"):
        data_with_timestamps = assign_test_timestamps(file_path)
        data_with_timestamps.to_csv(file_path.replace("manual", "external"), index=False)
        print(file_path + " was processed.")

if args.files:
    for file_name in args.files:
        file_path = "tests/data/manual/aware_csv/" + file_name
        data_with_timestamps = assign_test_timestamps(file_path)
        data_with_timestamps.to_csv(file_path.replace("manual", "external"), index=False)
        print(file_path + " was processed.")
