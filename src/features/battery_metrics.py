import pandas as pd
from datetime import datetime, timedelta, time

def truncateTime(df, segment_column, new_day_segment, datetime_column, date_column, new_time):
    df.loc[:, segment_column] = new_day_segment
    df.loc[:, datetime_column] = df[date_column].apply(lambda date: datetime.combine(date, new_time))
    return df

def computeTruncatedBatteryTimeDifferences(df):
    df["truncated_time_diff"] = df["local_end_date_time"] - df["local_start_date_time"]
    df["truncated_time_diff"] = df["truncated_time_diff"].apply(lambda time: time.total_seconds()/3600)
    df["battery_diff"] = df["battery_diff"] * (df["truncated_time_diff"] / df["time_diff"])
    del df["time_diff"]
    df.rename(columns={"truncated_time_diff": "time_diff"}, inplace=True)
    return df

def splitOvernightEpisodes(battery_data):
    overnight = battery_data[(battery_data["local_start_date"] + timedelta(days=1)) == battery_data["local_end_date"]]
    not_overnight = battery_data[battery_data["local_start_date"] == battery_data["local_end_date"]]

    if not overnight.empty:
        today = overnight[["battery_diff", "time_diff", "local_start_date_time", "local_start_date", "local_start_day_segment"]].copy()
        tomorrow = overnight[["battery_diff", "time_diff", "local_end_date_time", "local_end_date", "local_end_day_segment"]].copy()

        # truncate the end time of all overnight periods to midnight
        today = truncateTime(today, "local_end_day_segment", "evening", "local_end_date_time", "local_start_date", time(23,59,59))
        today["local_end_date"] = overnight["local_start_date"]
        
        # set the start time of all periods after midnight to midnight
        tomorrow = truncateTime(tomorrow, "local_start_day_segment", "night", "local_start_date_time", "local_end_date", time(0,0,0))
        tomorrow["local_start_date"] = overnight["local_end_date"]

        overnight = pd.concat([today, tomorrow], axis=0, sort=False)

        # calculate new battery_diff and time_diff for split overnight periods
        overnight = computeTruncatedBatteryTimeDifferences(overnight)

    return pd.concat([not_overnight, overnight], axis=0, sort=False)

def splitMultiSegmentEpisodes(battery_data, day_segment):
    # extract episodes that start and end at the same epochs
    exact_segments = battery_data.query("local_start_day_segment == local_end_day_segment and local_start_day_segment == @day_segment").copy()

    # extract episodes that start and end at different epochs
    across_segments = battery_data.query("local_start_day_segment != local_end_day_segment").copy()
    # 1) if start time is in current day_segment
    start_segment = across_segments[across_segments["local_start_day_segment"] == day_segment].copy()
    if not start_segment.empty:
        start_segment = truncateTime(start_segment, "local_end_day_segment", day_segment, "local_end_date_time", "local_end_date", time(EPOCH_TIMES[day_segment][1],59,59))
    # 2) if end time is in current day_segment
    end_segment = across_segments[across_segments["local_end_day_segment"] == day_segment].copy()
    if not end_segment.empty:
        end_segment = truncateTime(end_segment, "local_start_day_segment", day_segment, "local_start_date_time", "local_start_date", time(EPOCH_TIMES[day_segment][0],0,0))
    # 3) if current episode comtains day_segment
    across_segments.loc[:,"start_segment"] = across_segments["local_start_day_segment"].apply(lambda seg: SEGMENT[seg])
    across_segments.loc[:,"end_segment"] = across_segments["local_end_day_segment"].apply(lambda seg: SEGMENT[seg])
    day_segment_num = SEGMENT[day_segment]
    within_segments = across_segments.query("start_segment < @day_segment_num and end_segment > @day_segment_num")
    del across_segments["start_segment"], across_segments["end_segment"]
    del within_segments["start_segment"], within_segments["end_segment"]
    
    if not within_segments.empty:
        within_segments = truncateTime(within_segments, "local_start_day_segment", day_segment, "local_start_date_time", "local_start_date", time(EPOCH_TIMES[day_segment][0],0,0))

        within_segments = truncateTime(within_segments, "local_end_day_segment", day_segment, "local_end_date_time", "local_end_date", time(EPOCH_TIMES[day_segment][1],59,59))

    across_segments = pd.concat([start_segment, end_segment, within_segments], axis=0, sort=False)

    if not across_segments.empty:
        accross_segments = computeTruncatedBatteryTimeDifferences(across_segments)

    return pd.concat([exact_segments, across_segments], axis=0, sort=False)


battery_data = pd.read_csv(snakemake.input[0], parse_dates=["local_start_date_time", "local_end_date_time", "local_start_date", "local_end_date"])
day_segment = snakemake.params["day_segment"]
metrics = snakemake.params["metrics"]
SEGMENT = {"night": 0, "morning": 1, "afternoon": 2, "evening": 3}
EPOCH_TIMES = {"night": [0,5], "morning": [6,11], "afternoon": [12,17], "evening": [18,23]}

if battery_data.empty:
    battery_features = pd.DataFrame(columns=["local_date"] + ["battery_" + day_segment + "_" + x for x in metrics])
else:
    battery_data = splitOvernightEpisodes(battery_data)

    if day_segment != "daily":
        battery_data = splitMultiSegmentEpisodes(battery_data, day_segment)

    battery_data["battery_consumption_rate"] = battery_data["battery_diff"] / battery_data["time_diff"]

    # for battery_data_discharge:
    battery_data_discharge = battery_data[battery_data["battery_diff"] > 0]
    battery_discharge_features = pd.DataFrame()
    if "countdischarge" in metrics:
        battery_discharge_features["battery_"+day_segment+"_countdischarge"] = battery_data_discharge.groupby(["local_start_date"])["local_start_date"].count()
    if "sumdurationdischarge" in metrics:
        battery_discharge_features["battery_"+day_segment+"_sumdurationdischarge"] = battery_data_discharge.groupby(["local_start_date"])["time_diff"].sum()
    if "avgconsumptionrate" in metrics:
        battery_discharge_features["battery_"+day_segment+"_avgconsumptionrate"] = battery_data_discharge.groupby(["local_start_date"])["battery_consumption_rate"].mean()
    if "maxconsumptionrate" in metrics:
        battery_discharge_features["battery_"+day_segment+"_maxconsumptionrate"] = battery_data_discharge.groupby(["local_start_date"])["battery_consumption_rate"].max()

    # for battery_data_charge:
    battery_data_charge = battery_data[battery_data["battery_diff"] <= 0]
    battery_charge_features = pd.DataFrame()
    if "countcharge" in metrics:
        battery_charge_features["battery_"+day_segment+"_countcharge"] = battery_data_charge.groupby(["local_start_date"])["local_start_date"].count()
    if "sumdurationcharge" in metrics:
        battery_charge_features["battery_"+day_segment+"_sumdurationcharge"] = battery_data_charge.groupby(["local_start_date"])["time_diff"].sum()

    # combine discharge features and charge features; fill the missing values with ZERO
    battery_features = pd.concat([battery_discharge_features, battery_charge_features], axis=1, sort=True).fillna(0)
    battery_features.index.rename("local_date", inplace=True)
battery_features.to_csv(snakemake.output[0], index=True)