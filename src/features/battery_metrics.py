import pandas as pd
import datetime


battery_data = pd.read_csv(snakemake.input[0])
if battery_data.empty:
    battery_features = pd.DataFrame(columns=["battery_diff", "time_diff", "battery_decrease_times","battery_consumption_rate", "local_date"])
else:
    for col in ["local_start_date_time", "local_end_date_time", "local_start_date", "local_end_date"]:
        battery_data[col] = pd.to_datetime(battery_data[col])

    # split the row into 2 rows when local_start_date + 1 = local_end_date
    battery_data_overnight = battery_data[battery_data["local_start_date"] + datetime.timedelta(days=1) == battery_data["local_end_date"]]
    if not battery_data_overnight.empty:
        battery_data_overnight_first, battery_data_overnight_second = pd.DataFrame(columns=battery_data.columns), pd.DataFrame(columns=battery_data.columns)
        battery_data_overnight_first["total_battery_diff"], battery_data_overnight_second["total_battery_diff"] = battery_data_overnight["battery_diff"], battery_data_overnight["battery_diff"]
        battery_data_overnight_first["total_time_diff"], battery_data_overnight_second["total_time_diff"] = battery_data_overnight["time_diff"], battery_data_overnight["time_diff"]
        # let start = start OR end = end, then fill the left col with the start+23:59:59.
        battery_data_overnight_first["local_start_date_time"] = battery_data_overnight["local_start_date_time"]
        battery_data_overnight_first["local_end_date_time"] = battery_data_overnight["local_start_date"].apply(lambda x: datetime.datetime.combine(x, datetime.time(23,59,59)))
        battery_data_overnight_first["local_start_date"] = battery_data_overnight["local_start_date"]
        battery_data_overnight_second["local_end_date_time"] = battery_data_overnight["local_end_date_time"]
        battery_data_overnight_second["local_start_date_time"] = battery_data_overnight["local_start_date"].apply(lambda x: datetime.datetime.combine(x, datetime.time(23,59,59)))
        battery_data_overnight_second["local_start_date"] = battery_data_overnight["local_end_date"]
        battery_data_overnight = pd.concat([battery_data_overnight_first, battery_data_overnight_second])
        # calculate battery_diff and time_diff
        battery_data_overnight["time_diff"] = (battery_data_overnight["local_end_date_time"]-battery_data_overnight["local_start_date_time"]).apply(lambda x: x.total_seconds()/3600)
        battery_data_overnight["battery_diff"] = battery_data_overnight["total_battery_diff"]*(battery_data_overnight["time_diff"]/battery_data_overnight["total_time_diff"])
        del battery_data_overnight["total_battery_diff"], battery_data_overnight["total_time_diff"]

    # filter out the rows when local_start_date + 1 < local_end_date
    battery_data = battery_data[battery_data["local_start_date"] == battery_data["local_end_date"]]

    # combine
    battery_data = pd.concat([battery_data, battery_data_overnight])

    # split into decrease table and charge table
    battery_data_decrease = battery_data[battery_data["battery_diff"] > 0]
    battery_data_charge = battery_data[battery_data["battery_diff"] <= 0]

    # for battery_data_decrease:
    battery_decrease_count = battery_data_decrease.groupby(["local_start_date"])["local_start_date"].count()
    battery_data_decrease = battery_data_decrease.groupby(["local_start_date"]).sum()
    battery_data_decrease["battery_decrease_count"] = battery_decrease_count
    battery_data_decrease["battery_decrease_duration"] = battery_data_decrease["time_diff"]
    battery_data_decrease["battery_consumption_rate"] = battery_data_decrease["battery_diff"]/battery_data_decrease["time_diff"]
    del battery_data_decrease["battery_diff"], battery_data_decrease["time_diff"]

    # for battery_data_charge:
    battery_charge_count = battery_data_charge.groupby(["local_start_date"])["local_start_date"].count()
    battery_data_charge = battery_data_charge.groupby(["local_start_date"]).sum()
    battery_data_charge["battery_charge_count"] = battery_charge_count
    battery_data_charge["battery_charge_duration"] = battery_data_charge["time_diff"]
    del battery_data_charge["battery_diff"], battery_data_charge["time_diff"]

    # combine decrease features and charge features
    battery_features = pd.concat([battery_data_decrease, battery_data_charge], axis=1, sort=True)
    battery_features["local_date"] = battery_features.index
    battery_features.reset_index(inplace=True, drop=True)
battery_features.to_csv(snakemake.output[0], index=False)