import pandas as pd
import numpy as np

def statsFeatures(steps_data, features_to_compute, features_type, steps_features):
    if features_type == "steps" or features_type == "sumsteps":
        col_name = "steps"
    elif features_type == "durationsedentarybout" or features_type == "durationactivebout":
        col_name = "duration"
    else:
        raise ValueError("features_type can only be one of ['steps', 'sumsteps', 'durationsedentarybout', 'durationactivebout'].")

    if ("summarycount" if features_type == "sumsteps" else "intradaycount") + features_type.replace("duration", "episode") in features_to_compute:
        steps_features["steps_rapids_" + ("summarycount" if features_type == "sumsteps" else "intradaycount") + features_type.replace("duration", "episode")] = steps_data.groupby(["local_segment"])[col_name].count()
    if ("summarysum" if features_type == "sumsteps" else "intradaysum") + features_type in features_to_compute:
        steps_features["steps_rapids_" + ("summarysum" if features_type == "sumsteps" else "intradaysum") + features_type] = steps_data.groupby(["local_segment"])[col_name].sum()
    if ("summarymax" if features_type == "sumsteps" else "intradaymax") + features_type in features_to_compute:
        steps_features["steps_rapids_" + ("summarymax" if features_type == "sumsteps" else "intradaymax") + features_type] = steps_data.groupby(["local_segment"])[col_name].max()
    if ("summarymin" if features_type == "sumsteps" else "intradaymin") + features_type in features_to_compute:
        steps_features["steps_rapids_" + ("summarymin" if features_type == "sumsteps" else "intradaymin") + features_type] = steps_data.groupby(["local_segment"])[col_name].min()
    if ("summaryavg" if features_type == "sumsteps" else "intradayavg") + features_type in features_to_compute:
        steps_features["steps_rapids_" + ("summaryavg" if features_type == "sumsteps" else "intradayavg") + features_type] = steps_data.groupby(["local_segment"])[col_name].mean()
    if ("summarymedian" if features_type == "sumsteps" else "intradaymedian") + features_type in features_to_compute:
        steps_features["steps_rapids_" + ("summarymedian" if features_type == "sumsteps" else "intradaymedian") + features_type] = steps_data.groupby(["local_segment"])[col_name].median()
    if ("summarystd" if features_type == "sumsteps" else "intradaystd") + features_type in features_to_compute:
        steps_features["steps_rapids_" + ("summarystd" if features_type == "sumsteps" else "intradaystd") + features_type] = steps_data.groupby(["local_segment"])[col_name].std()

    return steps_features

def getBouts(steps_data):

    # put consecutive rows into the same group if they have the same values of "isactivebout", "local_timezone", and "local_segment"
    steps_data["group_idx"] =  (steps_data[["isactivebout", "local_timezone", "local_segment"]].shift() != steps_data[["isactivebout", "local_timezone", "local_segment"]]).any(axis=1).cumsum()
    
    # get bouts: duration column contains the number of minutes (rows) of sedentary and active activity for each episode
    grouped = steps_data.groupby("group_idx")
    bouts = grouped["local_segment"].agg(duration="count")
    bouts[["local_segment", "isactivebout"]] = grouped[["local_segment", "isactivebout"]].first()

    return bouts

def extractStepsFeaturesFromSummaryData(steps_summary_data, summary_features_to_compute):
    steps_summary_features = pd.DataFrame()

    # statistics features of daily steps count
    steps_summary_features = statsFeatures(steps_summary_data, summary_features_to_compute, "sumsteps", steps_summary_features)

    steps_summary_features.reset_index(inplace=True)
    
    return steps_summary_features

def extractStepsFeaturesFromIntradayData(steps_intraday_data, threshold_active_bout, intraday_features_to_compute_steps, intraday_features_to_compute_sedentarybout, intraday_features_to_compute_activebout, steps_intraday_features):
    steps_intraday_features = pd.DataFrame()

    # statistics features of steps count
    steps_intraday_features = statsFeatures(steps_intraday_data, intraday_features_to_compute_steps, "steps", steps_intraday_features)

    # sedentary bout: less than THRESHOLD_ACTIVE_BOUT (default: 10) steps in a minute
    # active bout: greater or equal to THRESHOLD_ACTIVE_BOUT (default: 10) steps in a minute
    isactivebout = np.where(steps_intraday_data["steps"] < int(threshold_active_bout), 0, 1)
    steps_intraday_data = steps_intraday_data.assign(isactivebout = isactivebout)
    bouts = getBouts(steps_intraday_data)

    # statistics features of sedentary bout
    sedentary_bout = bouts[bouts["isactivebout"] == 0]
    steps_intraday_features = statsFeatures(sedentary_bout, intraday_features_to_compute_sedentarybout, "durationsedentarybout", steps_intraday_features)

    # statistics features of active bout
    active_bout = bouts[bouts["isactivebout"] == 1]
    steps_intraday_features = statsFeatures(active_bout, intraday_features_to_compute_activebout, "durationactivebout", steps_intraday_features)

    steps_intraday_features.reset_index(inplace=True)

    return steps_intraday_features


def rapids_features(sensor_data_files, day_segment, provider, filter_data_by_segment, *args, **kwargs):

    threshold_active_bout = provider["THRESHOLD_ACTIVE_BOUT"]
    include_zero_step_rows = provider["INCLUDE_ZERO_STEP_ROWS"]

    steps_summary_data = pd.read_csv(sensor_data_files["sensor_data"][0])
    steps_intraday_data = pd.read_csv(sensor_data_files["sensor_data"][1])

    requested_summary_features = ["summary" + x for x in provider["FEATURES"]["SUMMARY"]]
    requested_intraday_features = provider["FEATURES"]["INTRADAY"]

    requested_intraday_features_steps = ["intraday" + x + "steps" for x in requested_intraday_features["STEPS"]]
    requested_intraday_features_sedentarybout = ["intraday" + x + "sedentarybout" for x in requested_intraday_features["SEDENTARY_BOUT"]]
    requested_intraday_features_activebout = ["intraday" + x + "activebout" for x in requested_intraday_features["ACTIVE_BOUT"]]
    # name of the features this function can compute
    base_summary_features = ["summarymaxsumsteps", "summaryminsumsteps", "summaryavgsumsteps", "summarymediansumsteps", "summarystdsumsteps"]
    base_intraday_features_steps = ["intradaysumsteps", "intradaymaxsteps", "intradayminsteps", "intradayavgsteps", "intradaystdsteps"]
    base_intraday_features_sedentarybout = ["intradaycountepisodesedentarybout", "intradaysumdurationsedentarybout", "intradaymaxdurationsedentarybout", "intradaymindurationsedentarybout", "intradayavgdurationsedentarybout", "intradaystddurationsedentarybout"]
    base_intraday_features_activebout = ["intradaycountepisodeactivebout", "intradaysumdurationactivebout", "intradaymaxdurationactivebout", "intradaymindurationactivebout", "intradayavgdurationactivebout", "intradaystddurationactivebout"]
    # the subset of requested features this function can compute
    intraday_features_to_compute_steps = list(set(requested_intraday_features_steps) & set(base_intraday_features_steps))
    intraday_features_to_compute_sedentarybout = list(set(requested_intraday_features_sedentarybout) & set(base_intraday_features_sedentarybout))
    intraday_features_to_compute_activebout = list(set(requested_intraday_features_activebout) & set(base_intraday_features_activebout))

    summary_features_to_compute = list(set(requested_summary_features) & set(base_summary_features))
    intraday_features_to_compute = intraday_features_to_compute_steps + intraday_features_to_compute_sedentarybout + intraday_features_to_compute_activebout

    # extract features from summary data
    steps_summary_features = pd.DataFrame(columns=["local_segment"] + ["steps_rapids_" + x for x in summary_features_to_compute])
    if not steps_summary_data.empty:
        steps_summary_data = filter_data_by_segment(steps_summary_data, day_segment)

        if not steps_summary_data.empty:
            # only keep the segments start at 00:00:00 and end at 23:59:59
            datetime_start_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2} 00:00:00"
            datetime_end_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2} 23:59:59"

            segment_regex = "{}#{},{}".format(day_segment, datetime_start_regex, datetime_end_regex)
            steps_summary_data = steps_summary_data[steps_summary_data["local_segment"].str.match(segment_regex)]

            if not steps_summary_data.empty:
                steps_summary_features = extractStepsFeaturesFromSummaryData(steps_summary_data, summary_features_to_compute)

    # extract features from intraday features
    steps_intraday_features = pd.DataFrame(columns=["local_segment"] + ["steps_rapids_" + x for x in intraday_features_to_compute])
    if not steps_intraday_data.empty:
        steps_intraday_data = filter_data_by_segment(steps_intraday_data, day_segment)

        if not steps_intraday_data.empty:
            steps_intraday_features = extractStepsFeaturesFromIntradayData(steps_intraday_data, threshold_active_bout, intraday_features_to_compute_steps, intraday_features_to_compute_sedentarybout, intraday_features_to_compute_activebout, steps_intraday_features)

    # merge summary features and intraday features
    steps_features = steps_intraday_features.merge(steps_summary_features, on=["local_segment"], how="outer")


    # exclude rows when the total step count is ZERO during the whole day
    if not include_zero_step_rows:
        steps_features.index = steps_features["local_segment"].apply(lambda segment: segment.split("#")[1][:10])

        steps_features["dailycountstep"] = steps_intraday_data.groupby(["local_date"])["steps"].sum()
        steps_features = steps_features.query("dailycountstep != 0")

        del steps_features["dailycountstep"]
        steps_features.reset_index(drop=True, inplace=True)
    
    return steps_features
