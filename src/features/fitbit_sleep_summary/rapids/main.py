import pandas as pd
import numpy as np
import itertools

def extractSleepFeaturesFromSummaryData(sleep_summary_data, summary_features, sleep_type, waketime_features, sleep_summary_features):
    if sleep_type == "main":
        sleep_summary_data = sleep_summary_data[sleep_summary_data["is_main_sleep"] == 1]
    elif sleep_type == "nap":
        sleep_summary_data = sleep_summary_data[sleep_summary_data["is_main_sleep"] == 0]
    elif sleep_type == "all":
        pass
    else:
        raise ValueError("sleep_type can only be one of ['main', 'nap', 'all'].")

    features_sum = sleep_summary_data[["local_segment", "minutes_after_wakeup", "minutes_asleep", "minutes_awake", "minutes_to_fall_asleep", "minutes_in_bed"]].groupby(["local_segment"]).sum(min_count=1)

    if "sumdurationafterwakeup" in summary_features:
        sleep_summary_features = sleep_summary_features.join(features_sum[["minutes_after_wakeup"]], how="outer").rename(columns={"minutes_after_wakeup": "sumdurationafterwakeup" + sleep_type})
    if "sumdurationasleep" in summary_features:
        sleep_summary_features = sleep_summary_features.join(features_sum[["minutes_asleep"]], how="outer").rename(columns={"minutes_asleep": "sumdurationasleep" + sleep_type})
    if "sumdurationawake" in summary_features:
        sleep_summary_features = sleep_summary_features.join(features_sum[["minutes_awake"]], how="outer").rename(columns={"minutes_awake": "sumdurationawake" + sleep_type})
    if "sumdurationtofallasleep" in summary_features:
        sleep_summary_features = sleep_summary_features.join(features_sum[["minutes_to_fall_asleep"]], how="outer").rename(columns={"minutes_to_fall_asleep": "sumdurationtofallasleep" + sleep_type})
    if "sumdurationinbed" in summary_features:
        sleep_summary_features = sleep_summary_features.join(features_sum[["minutes_in_bed"]], how="outer").rename(columns={"minutes_in_bed": "sumdurationinbed" + sleep_type})

    features_avg = sleep_summary_data[["local_segment", "efficiency", "minutes_after_wakeup", "minutes_asleep", "minutes_awake", "minutes_to_fall_asleep", "minutes_in_bed"]].groupby(["local_segment"]).mean()

    if "avgefficiency" in summary_features:
        sleep_summary_features = sleep_summary_features.join(features_avg[["efficiency"]], how="outer").rename(columns={"efficiency": "avgefficiency" + sleep_type})
    if "avgdurationafterwakeup" in summary_features:
        sleep_summary_features = sleep_summary_features.join(features_avg[["minutes_after_wakeup"]], how="outer").rename(columns={"minutes_after_wakeup": "avgdurationafterwakeup" + sleep_type})
    if "avgdurationasleep" in summary_features:
        sleep_summary_features = sleep_summary_features.join(features_avg[["minutes_asleep"]], how="outer").rename(columns={"minutes_asleep": "avgdurationasleep" + sleep_type})
    if "avgdurationawake" in summary_features:
        sleep_summary_features = sleep_summary_features.join(features_avg[["minutes_awake"]], how="outer").rename(columns={"minutes_awake": "avgdurationawake" + sleep_type})
    if "avgdurationtofallasleep" in summary_features:
        sleep_summary_features = sleep_summary_features.join(features_avg[["minutes_to_fall_asleep"]], how="outer").rename(columns={"minutes_to_fall_asleep": "avgdurationtofallasleep" + sleep_type})
    if "avgdurationinbed" in summary_features:
        sleep_summary_features = sleep_summary_features.join(features_avg[["minutes_in_bed"]], how="outer").rename(columns={"minutes_in_bed": "avgdurationinbed" + sleep_type})


    if "countepisode" in summary_features:
        features_count = sleep_summary_data[["local_segment", "timestamp"]].groupby(["local_segment"]).count()
        sleep_summary_features = sleep_summary_features.join(features_count[["timestamp"]], how="outer").rename(columns={"timestamp": "countepisode" + sleep_type})
    
    
    if "firstbedtime" in summary_features:
        features_first = sleep_summary_data[["local_segment", "minutes_start_episode"]].groupby(["local_segment"]).min()
        sleep_summary_features = sleep_summary_features.join(features_first[["minutes_start_episode"]], how="outer").rename(columns={"minutes_start_episode": "firstbedtime" + sleep_type})
    if "lastbedtime" in summary_features:
        features_last = sleep_summary_data[["local_segment", "minutes_start_episode"]].groupby(["local_segment"]).max()
        sleep_summary_features = sleep_summary_features.join(features_last[["minutes_start_episode"]], how="outer").rename(columns={"minutes_start_episode": "lastbedtime" + sleep_type})

    if "firstwaketime" in summary_features:
        sleep_summary_features = sleep_summary_features.join(waketime_features[["firstwaketime" + sleep_type]], how="outer")
    if "lastwaketime" in summary_features:
        sleep_summary_features = sleep_summary_features.join(waketime_features[["lastwaketime" + sleep_type]], how="outer")


    return sleep_summary_features


def rapids_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):

    sleep_summary_data = pd.read_csv(sensor_data_files["sensor_data"])

    requested_summary_features = provider["FEATURES"]
    requested_sleep_types = provider["SLEEP_TYPES"]

    # name of the features this function can compute
    base_summary_features = ["firstwaketime", "lastwaketime", "firstbedtime", "lastbedtime", "countepisode", "avgefficiency", "sumdurationafterwakeup", "sumdurationasleep", "sumdurationawake", "sumdurationtofallasleep", "sumdurationinbed", "avgdurationafterwakeup", "avgdurationasleep", "avgdurationawake", "avgdurationtofallasleep", "avgdurationinbed"]
    base_sleep_types = ["main", "nap", "all"]
    # the subset of requested features this function can compute
    summary_features_to_compute = list(set(requested_summary_features) & set(base_summary_features))
    sleep_types_to_compute = list(set(requested_sleep_types) & set(base_sleep_types))
    # full names
    features_fullnames_to_compute = ["".join(feature) for feature in itertools.product(summary_features_to_compute, sleep_types_to_compute)]
    
    colnames_can_be_zero = ["".join(feature) for feature in itertools.product(set(summary_features_to_compute) - set(["firstwaketime", "lastwaketime", "firstbedtime", "lastbedtime", "avgefficiency"]), sleep_types_to_compute)]
    
    # extract features from summary data
    sleep_summary_features = pd.DataFrame(columns=["local_segment"] + features_fullnames_to_compute)
    if not sleep_summary_data.empty:

        # calculate number of minutes after segment's start date time
        dt_cols = ["local_start_date_time", "local_end_date_time", "local_date_time"]
        sleep_summary_data[dt_cols] = sleep_summary_data[dt_cols].apply(pd.to_datetime)
        sleep_summary_data["minutes_start_episode"] = (sleep_summary_data["local_start_date_time"] - sleep_summary_data["local_date_time"]) / pd.Timedelta(minutes=1)
        sleep_summary_data["minutes_end_episode"] = (sleep_summary_data["local_end_date_time"] - (sleep_summary_data["local_date_time"] + pd.Timedelta(days=1))) / pd.Timedelta(minutes=1)

        # get minutes_end_episode for different sleep types
        sleep_summary_data["minutes_end_episode_main"] = sleep_summary_data.apply(lambda row: row["minutes_end_episode"] if row["is_main_sleep"] == 1 else np.nan, axis=1)
        sleep_summary_data["minutes_end_episode_nap"] = sleep_summary_data.apply(lambda row: row["minutes_end_episode"] if row["is_main_sleep"] == 0 else np.nan, axis=1)

        # extract daily wake time
        sleep_summary_data = sleep_summary_data.merge(sleep_summary_data.groupby("local_date_time")[["minutes_end_episode", "minutes_end_episode_main", "minutes_end_episode_nap"]].agg(
            firstwaketimeall=("minutes_end_episode", "min"),
            lastwaketimeall=("minutes_end_episode", "max"),
            firstwaketimemain=("minutes_end_episode_main", "min"),
            lastwaketimemain=("minutes_end_episode_main", "max"),
            firstwaketimenap=("minutes_end_episode_nap", "min"),
            lastwaketimenap=("minutes_end_episode_nap", "max"),
            ).shift(), right_index=True, left_on="local_date_time")

        # filter by segment
        sleep_summary_data = filter_data_by_segment(sleep_summary_data, time_segment)

        notna_segments = sleep_summary_data[sleep_summary_data["type"].notna()]["local_segment"].unique()

        if not sleep_summary_data.empty:
            # only keep the segments start at 00:00:00 and end at 23:59:59
            datetime_start_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2} 00:00:00"
            datetime_end_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2} 23:59:59"

            segment_regex = "{}#{},{}".format(time_segment, datetime_start_regex, datetime_end_regex)
            sleep_summary_data = sleep_summary_data[sleep_summary_data["local_segment"].str.match(segment_regex)]

            if not sleep_summary_data.empty:

                waketime_features = sleep_summary_data.drop_duplicates(subset=["local_date_time"]).groupby(["local_segment"]).agg(
                    firstwaketimeall=("firstwaketimeall", "mean"),
                    lastwaketimeall=("lastwaketimeall", "mean"),
                    firstwaketimemain=("firstwaketimemain", "mean"),
                    lastwaketimemain=("lastwaketimemain", "mean"),
                    firstwaketimenap=("firstwaketimenap", "mean"),
                    lastwaketimenap=("lastwaketimenap", "mean"),
                )
                sleep_summary_data.dropna(subset=["is_main_sleep"], inplace=True)

                if not sleep_summary_data.empty:
                    sleep_summary_features = pd.DataFrame()

                    for sleep_type in sleep_types_to_compute:
                        sleep_summary_features = extractSleepFeaturesFromSummaryData(sleep_summary_data, summary_features_to_compute, sleep_type, waketime_features, sleep_summary_features)
                    
                    sleep_summary_features.dropna(how="all", inplace=True)
                    sleep_summary_features.loc[notna_segments, colnames_can_be_zero] = sleep_summary_features.loc[notna_segments, colnames_can_be_zero].fillna(0)

                    sleep_summary_features = sleep_summary_features.reset_index()
    
    return sleep_summary_features
