import pandas as pd
import itertools



def dailyFeaturesFromSummaryData(sleep_summary_data, sleep_type):
    if sleep_type == "main":
        sleep_summary_data = sleep_summary_data[sleep_summary_data["is_main_sleep"] == 1]
    elif sleep_type == "nap":
        sleep_summary_data = sleep_summary_data[sleep_summary_data["is_main_sleep"] == 0]
    elif sleep_type == "all":
        pass
    else:
        raise ValueError("sleep_type can only be one of ['main', 'nap', 'all'].")

    features_sum = sleep_summary_data[["minutes_after_wakeup", "minutes_asleep", "minutes_awake", "minutes_to_fall_asleep", "minutes_in_bed", "local_end_date"]].groupby(["local_end_date"]).sum()
    features_sum.index.rename("local_date", inplace=True)
    if "sumdurationafterwakeup" in daily_features_from_summary_data:
        sleep_daily_features["sleep_daily_sumdurationafterwakeup" + sleep_type] = features_sum["minutes_after_wakeup"]
    if "sumdurationasleep" in daily_features_from_summary_data:
        sleep_daily_features["sleep_daily_sumdurationasleep" + sleep_type] = features_sum["minutes_asleep"]
    if "sumdurationawake" in daily_features_from_summary_data:
        sleep_daily_features["sleep_daily_sumdurationawake" + sleep_type] = features_sum["minutes_awake"]
    if "sumdurationtofallasleep" in daily_features_from_summary_data:
        sleep_daily_features["sleep_daily_sumdurationtofallasleep" + sleep_type] = features_sum["minutes_to_fall_asleep"]
    if "sumdurationinbed" in daily_features_from_summary_data:
        sleep_daily_features["sleep_daily_sumdurationinbed" + sleep_type] = features_sum["minutes_in_bed"]

    features_avg = sleep_summary_data[["efficiency", "local_end_date"]].groupby(["local_end_date"]).mean()
    features_avg.index.rename("local_date", inplace=True)
    if "avgefficiency" in daily_features_from_summary_data:
        sleep_daily_features["sleep_daily_avgefficiency" + sleep_type] = features_avg["efficiency"]
    
    features_count = sleep_summary_data[["local_start_date_time", "local_end_date"]].groupby(["local_end_date"]).count()
    features_count.index.rename("local_date", inplace=True)
    if "countepisode" in daily_features_from_summary_data:
        sleep_daily_features["sleep_daily_count" + sleep_type] = features_count["local_start_date_time"]

    return sleep_daily_features



sleep_summary_data = pd.read_csv(snakemake.input["sleep_summary_data"])
sleep_types = snakemake.params["sleep_types"]
daily_features_from_summary_data = snakemake.params["daily_features_from_summary_data"]
day_segment = snakemake.params["day_segment"]

daily_features_can_be_zero = list(set(daily_features_from_summary_data) - set(["avgefficiency"]))
colnames_can_be_zero = ["sleep_daily_" + x for x in ["".join(feature) for feature in itertools.product(daily_features_can_be_zero, sleep_types)]]

colnames = ["sleep_daily_" + x for x in ["".join(feature) for feature in itertools.product(daily_features_from_summary_data, sleep_types)]]

if sleep_summary_data.empty:
    sleep_daily_features = pd.DataFrame(columns=["local_date"] + colnames)
else:
    sleep_daily_features = pd.DataFrame(columns=colnames)
    for sleep_type in sleep_types:
        sleep_daily_features = dailyFeaturesFromSummaryData(sleep_summary_data, sleep_type)

    sleep_daily_features[colnames_can_be_zero] = sleep_daily_features[colnames_can_be_zero].fillna(0)



if day_segment == "daily":
    sleep_daily_features.to_csv(snakemake.output[0])
else:
    pd.DataFrame().to_csv(snakemake.output[0])
