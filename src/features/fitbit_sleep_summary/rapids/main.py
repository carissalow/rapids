import pandas as pd
import itertools

def extractSleepFeaturesFromSummaryData(sleep_summary_data, summary_features, sleep_type, sleep_summary_features):
    if sleep_type == "main":
        sleep_summary_data = sleep_summary_data[sleep_summary_data["is_main_sleep"] == 1]
    elif sleep_type == "nap":
        sleep_summary_data = sleep_summary_data[sleep_summary_data["is_main_sleep"] == 0]
    elif sleep_type == "all":
        pass
    else:
        raise ValueError("sleep_type can only be one of ['main', 'nap', 'all'].")

    features_sum = sleep_summary_data[["local_segment", "minutes_after_wakeup", "minutes_asleep", "minutes_awake", "minutes_to_fall_asleep", "minutes_in_bed"]].groupby(["local_segment"]).sum()

    if "summarysumdurationafterwakeup" in summary_features:
        sleep_summary_features = sleep_summary_features.join(features_sum[["minutes_after_wakeup"]], how="outer").rename(columns={"minutes_after_wakeup": "sleep_rapids_summarysumdurationafterwakeup" + sleep_type})
    if "summarysumdurationasleep" in summary_features:
        sleep_summary_features = sleep_summary_features.join(features_sum[["minutes_asleep"]], how="outer").rename(columns={"minutes_asleep": "sleep_rapids_summarysumdurationasleep" + sleep_type})
    if "summarysumdurationawake" in summary_features:
        sleep_summary_features = sleep_summary_features.join(features_sum[["minutes_awake"]], how="outer").rename(columns={"minutes_awake": "sleep_rapids_summarysumdurationawake" + sleep_type})
    if "summarysumdurationtofallasleep" in summary_features:
        sleep_summary_features = sleep_summary_features.join(features_sum[["minutes_to_fall_asleep"]], how="outer").rename(columns={"minutes_to_fall_asleep": "sleep_rapids_summarysumdurationtofallasleep" + sleep_type})
    if "summarysumdurationinbed" in summary_features:
        sleep_summary_features = sleep_summary_features.join(features_sum[["minutes_in_bed"]], how="outer").rename(columns={"minutes_in_bed": "sleep_rapids_summarysumdurationinbed" + sleep_type})

    features_avg = sleep_summary_data[["local_segment", "efficiency", "minutes_after_wakeup", "minutes_asleep", "minutes_awake", "minutes_to_fall_asleep", "minutes_in_bed"]].groupby(["local_segment"]).mean()

    if "summaryavgefficiency" in summary_features:
        sleep_summary_features = sleep_summary_features.join(features_avg[["efficiency"]], how="outer").rename(columns={"efficiency": "sleep_rapids_summaryavgefficiency" + sleep_type})
    if "summaryavgdurationafterwakeup" in summary_features:
        sleep_summary_features = sleep_summary_features.join(features_avg[["minutes_after_wakeup"]], how="outer").rename(columns={"minutes_after_wakeup": "sleep_rapids_summaryavgdurationafterwakeup" + sleep_type})
    if "summaryavgdurationasleep" in summary_features:
        sleep_summary_features = sleep_summary_features.join(features_avg[["minutes_asleep"]], how="outer").rename(columns={"minutes_asleep": "sleep_rapids_summaryavgdurationasleep" + sleep_type})
    if "summaryavgdurationawake" in summary_features:
        sleep_summary_features = sleep_summary_features.join(features_avg[["minutes_awake"]], how="outer").rename(columns={"minutes_awake": "sleep_rapids_summaryavgdurationawake" + sleep_type})
    if "summaryavgdurationtofallasleep" in summary_features:
        sleep_summary_features = sleep_summary_features.join(features_avg[["minutes_to_fall_asleep"]], how="outer").rename(columns={"minutes_to_fall_asleep": "sleep_rapids_summaryavgdurationtofallasleep" + sleep_type})
    if "summaryavgdurationinbed" in summary_features:
        sleep_summary_features = sleep_summary_features.join(features_avg[["minutes_in_bed"]], how="outer").rename(columns={"minutes_in_bed": "sleep_rapids_summaryavgdurationinbed" + sleep_type})
    
    features_count = sleep_summary_data[["local_segment", "timestamp"]].groupby(["local_segment"]).count()
    
    if "summarycountepisode" in summary_features:
        sleep_summary_features = sleep_summary_features.join(features_count[["timestamp"]], how="outer").rename(columns={"timestamp": "sleep_rapids_summarycountepisode" + sleep_type})

    return sleep_summary_features


def rapids_features(sensor_data_files, day_segment, provider, filter_data_by_segment, *args, **kwargs):

    sleep_summary_data = pd.read_csv(sensor_data_files["sensor_data"])

    requested_summary_features = ["summary" + x for x in provider["FEATURES"]]
    requested_sleep_types = provider["SLEEP_TYPES"]

    # name of the features this function can compute
    base_summary_features = ["summarycountepisode", "summaryavgefficiency", "summarysumdurationafterwakeup", "summarysumdurationasleep", "summarysumdurationawake", "summarysumdurationtofallasleep", "summarysumdurationinbed", "summaryavgdurationafterwakeup", "summaryavgdurationasleep", "summaryavgdurationawake", "summaryavgdurationtofallasleep", "summaryavgdurationinbed"]
    base_sleep_types = ["main", "nap", "all"]
    # the subset of requested features this function can compute
    summary_features_to_compute = list(set(requested_summary_features) & set(base_summary_features))
    sleep_types_to_compute = list(set(requested_sleep_types) & set(base_sleep_types))
    # full names
    features_fullnames_to_compute = ["".join(feature) for feature in itertools.product(summary_features_to_compute, sleep_types_to_compute)]
    
    colnames_can_be_zero = ["sleep_rapids_" + x for x in [col for col in features_fullnames_to_compute if "summaryavgefficiency" not in col]]
    
    # extract features from summary data
    sleep_summary_features = pd.DataFrame(columns=["local_segment"] + ["sleep_rapids_" + x for x in features_fullnames_to_compute])
    if not sleep_summary_data.empty:
        sleep_summary_data = filter_data_by_segment(sleep_summary_data, day_segment)

        if not sleep_summary_data.empty:
            # only keep the segments start at 00:00:00 and end at 23:59:59
            datetime_start_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2} 00:00:00"
            datetime_end_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2} 23:59:59"

            segment_regex = "{}#{},{}".format(day_segment, datetime_start_regex, datetime_end_regex)
            sleep_summary_data = sleep_summary_data[sleep_summary_data["local_segment"].str.match(segment_regex)]

            if not sleep_summary_data.empty:
                sleep_summary_features = pd.DataFrame()

                for sleep_type in sleep_types_to_compute:
                    sleep_summary_features = extractSleepFeaturesFromSummaryData(sleep_summary_data, summary_features_to_compute, sleep_type, sleep_summary_features)

                sleep_summary_features[colnames_can_be_zero] = sleep_summary_features[colnames_can_be_zero].fillna(0)

                sleep_summary_features = sleep_summary_features.reset_index()
    
    return sleep_summary_features
