import pandas as pd
import numpy as np

def statsFeatures(steps_data, features_to_compute, features_type, steps_features, *args, **kwargs):
    if features_type == "steps" or features_type == "sumsteps":
        col_name = "steps"
        reference_hour = kwargs["reference_hour"]
    elif features_type == "durationsedentarybout" or features_type == "durationactivebout":
        col_name = "duration"
    else:
        raise ValueError("features_type can only be one of ['steps', 'sumsteps', 'durationsedentarybout', 'durationactivebout'].")

    if "count" + features_type.replace("duration", "episode") in features_to_compute:
        steps_features["count" + features_type.replace("duration", "episode")] = steps_data.groupby(["local_segment"])[col_name].count()
    if "sum" + features_type in features_to_compute:
        steps_features["sum" + features_type] = steps_data.groupby(["local_segment"])[col_name].sum()
    if "max" + features_type in features_to_compute:
        steps_features["max" + features_type] = steps_data.groupby(["local_segment"])[col_name].max()
    if "min" + features_type in features_to_compute:
        steps_features["min" + features_type] = steps_data.groupby(["local_segment"])[col_name].min()
    if "avg" + features_type in features_to_compute:
        steps_features["avg" + features_type] = steps_data.groupby(["local_segment"])[col_name].mean()
    if "median" + features_type in features_to_compute:
        steps_features["median" + features_type] = steps_data.groupby(["local_segment"])[col_name].median()
    if "std" + features_type in features_to_compute:
        steps_features["std" + features_type] = steps_data.groupby(["local_segment"])[col_name].std()
    if (col_name == "steps") and ("firststeptime" in features_to_compute):
        steps_features["firststeptime"] = steps_data[steps_data["steps"].ne(0)].groupby(["local_segment"])["local_time"].first().apply(lambda x: (int(x.split(":")[0]) - reference_hour) * 60 + int(x.split(":")[1]) + (int(x.split(":")[2]) / 60))
    if (col_name == "steps") and ("laststeptime" in features_to_compute):
        steps_features["laststeptime"] = steps_data[steps_data["steps"].ne(0)].groupby(["local_segment"])["local_time"].last().apply(lambda x: (int(x.split(":")[0]) - reference_hour) * 60 + int(x.split(":")[1]) + (int(x.split(":")[2]) / 60))

    return steps_features

def getBouts(steps_data):

    # put consecutive rows into the same group if they have the same values of "isactivebout", "local_timezone", and "local_segment"
    steps_data["group_idx"] =  (steps_data[["isactivebout", "local_timezone", "local_segment"]].shift() != steps_data[["isactivebout", "local_timezone", "local_segment"]]).any(axis=1).cumsum()
    
    # get bouts: duration column contains the number of minutes (rows) of sedentary and active activity for each episode
    grouped = steps_data.groupby("group_idx")
    bouts = grouped["local_segment"].agg(duration="count")
    bouts[["local_segment", "isactivebout"]] = grouped[["local_segment", "isactivebout"]].first()

    return bouts

def extractStepsFeaturesFromIntradayData(steps_intraday_data, reference_hour, threshold_active_bout, intraday_features_to_compute_steps, intraday_features_to_compute_sedentarybout, intraday_features_to_compute_activebout, steps_intraday_features):
    steps_intraday_features = pd.DataFrame()

    # statistics features of steps count
    steps_intraday_features = statsFeatures(steps_intraday_data, intraday_features_to_compute_steps, "steps", steps_intraday_features, reference_hour=reference_hour)

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



def rapids_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):

    reference_hour = provider["REFERENCE_HOUR"]
    threshold_active_bout = provider["THRESHOLD_ACTIVE_BOUT"]
    include_zero_step_rows = provider["INCLUDE_ZERO_STEP_ROWS"]

    steps_intraday_data = pd.read_csv(sensor_data_files["sensor_data"])

    requested_intraday_features = provider["FEATURES"]

    requested_intraday_features_steps = [x + "steps" if x not in ["firststeptime", "laststeptime"] else x for x in requested_intraday_features["STEPS"]]
    requested_intraday_features_sedentarybout = [x + "sedentarybout" for x in requested_intraday_features["SEDENTARY_BOUT"]]
    requested_intraday_features_activebout = [x + "activebout" for x in requested_intraday_features["ACTIVE_BOUT"]]
    # name of the features this function can compute
    base_intraday_features_steps = ["sumsteps", "maxsteps", "minsteps", "avgsteps", "stdsteps", "firststeptime", "laststeptime"]
    base_intraday_features_sedentarybout = ["countepisodesedentarybout", "sumdurationsedentarybout", "maxdurationsedentarybout", "mindurationsedentarybout", "avgdurationsedentarybout", "stddurationsedentarybout"]
    base_intraday_features_activebout = ["countepisodeactivebout", "sumdurationactivebout", "maxdurationactivebout", "mindurationactivebout", "avgdurationactivebout", "stddurationactivebout"]
    # the subset of requested features this function can compute
    intraday_features_to_compute_steps = list(set(requested_intraday_features_steps) & set(base_intraday_features_steps))
    intraday_features_to_compute_sedentarybout = list(set(requested_intraday_features_sedentarybout) & set(base_intraday_features_sedentarybout))
    intraday_features_to_compute_activebout = list(set(requested_intraday_features_activebout) & set(base_intraday_features_activebout))

    intraday_features_to_compute = intraday_features_to_compute_steps + intraday_features_to_compute_sedentarybout + intraday_features_to_compute_activebout

    # exclude rows when the total step count is ZERO during the whole day
    if (not steps_intraday_data.empty) and (not include_zero_step_rows):
        dailycountstep = steps_intraday_data.groupby(["local_date"])[["steps"]].sum()
        zerocountdates = dailycountstep[dailycountstep["steps"] == 0].index.tolist()
        steps_intraday_data = steps_intraday_data[~steps_intraday_data["local_date"].isin(zerocountdates)]

    # extract features from intraday features
    steps_intraday_features = pd.DataFrame(columns=["local_segment"] + intraday_features_to_compute)
    if not steps_intraday_data.empty:
        steps_intraday_data = filter_data_by_segment(steps_intraday_data, time_segment)

        if not steps_intraday_data.empty:
            steps_intraday_features = extractStepsFeaturesFromIntradayData(steps_intraday_data, reference_hour, threshold_active_bout, intraday_features_to_compute_steps, intraday_features_to_compute_sedentarybout, intraday_features_to_compute_activebout, steps_intraday_features)
    
    return steps_intraday_features
