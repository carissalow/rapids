import pandas as pd
import itertools



def dailyFeaturesFromSummaryData(sleep_daily_features, sleep_summary_data, summary_features, sleep_type):
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
    if "sumdurationafterwakeup" in summary_features:
        sleep_daily_features = sleep_daily_features.join(features_sum[["minutes_after_wakeup"]], how="outer").rename(columns={"minutes_after_wakeup": "sleep_daily_sumdurationafterwakeup" + sleep_type})
    if "sumdurationasleep" in summary_features:
        sleep_daily_features = sleep_daily_features.join(features_sum[["minutes_asleep"]], how="outer").rename(columns={"minutes_asleep": "sleep_daily_sumdurationasleep" + sleep_type})
    if "sumdurationawake" in summary_features:
        sleep_daily_features = sleep_daily_features.join(features_sum[["minutes_awake"]], how="outer").rename(columns={"minutes_awake": "sleep_daily_sumdurationawake" + sleep_type})
    if "sumdurationtofallasleep" in summary_features:
        sleep_daily_features = sleep_daily_features.join(features_sum[["minutes_to_fall_asleep"]], how="outer").rename(columns={"minutes_to_fall_asleep": "sleep_daily_sumdurationtofallasleep" + sleep_type})
    if "sumdurationinbed" in summary_features:
        sleep_daily_features = sleep_daily_features.join(features_sum[["minutes_in_bed"]], how="outer").rename(columns={"minutes_in_bed": "sleep_daily_sumdurationinbed" + sleep_type})

    features_avg = sleep_summary_data[["efficiency", "local_end_date"]].groupby(["local_end_date"]).mean()
    features_avg.index.rename("local_date", inplace=True)
    if "avgefficiency" in summary_features:
        sleep_daily_features = sleep_daily_features.join(features_avg[["efficiency"]], how="outer").rename(columns={"efficiency": "sleep_daily_avgefficiency" + sleep_type})
    
    features_count = sleep_summary_data[["local_start_date_time", "local_end_date"]].groupby(["local_end_date"]).count()
    features_count.index.rename("local_date", inplace=True)
    if "countepisode" in summary_features:
        sleep_daily_features = sleep_daily_features.join(features_count[["local_start_date_time"]], how="outer").rename(columns={"local_start_date_time": "sleep_daily_countepisode" + sleep_type})

    return sleep_daily_features

def base_fitbit_sleep_features(sleep_summary_data, day_segment, requested_summary_features, requested_sleep_type):
    if not day_segment == "daily":
        return pd.DataFrame(columns=["local_date"])
    else:
        # name of the features this function can compute
        base_summary_features_names = ["sumdurationafterwakeup", "sumdurationasleep", "sumdurationawake", "sumdurationtofallasleep", "sumdurationinbed", "avgefficiency", "countepisode"]
        base_sleep_type = ["main", "nap", "all"]
        # the subset of requested features this function can compute
        summary_features_to_compute = list(set(requested_summary_features) & set(base_summary_features_names))
        sleep_type_to_compute = list(set(requested_sleep_type) & set(base_sleep_type))
        # full names
        features_fullnames_to_compute = ["".join(feature) for feature in itertools.product(summary_features_to_compute, sleep_type_to_compute)]

        colnames_can_be_zero = ["sleep_daily_" + x for x in [col for col in features_fullnames_to_compute if "avgefficiency" not in col]]

        if sleep_summary_data.empty:
            sleep_summary_features = pd.DataFrame(columns=["local_date"] + ["sleep_daily_" + x for x in features_fullnames_to_compute])
        else:

            sleep_summary_features = pd.DataFrame()

            for sleep_type in sleep_type_to_compute:
                sleep_summary_features = dailyFeaturesFromSummaryData(sleep_summary_features, sleep_summary_data, summary_features_to_compute, sleep_type)
            
            sleep_summary_features[colnames_can_be_zero] = sleep_summary_features[colnames_can_be_zero].fillna(0)

            sleep_summary_features = sleep_summary_features.reset_index()
        
        return sleep_summary_features

