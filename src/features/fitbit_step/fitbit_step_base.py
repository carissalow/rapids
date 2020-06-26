import pandas as pd
import numpy as np

def getBouts(step_data, time_interval):
    # resample the data into time_interval minute bins, set "isactivebout" column to be NA if it is missing
    resampled_step_minute = pd.DataFrame(step_data.resample(str(time_interval) + "T", on="local_date_time")["isactivebout"].sum(min_count=1))

    # group rows by consecutive values of "isactivebout" column
    group = pd.DataFrame(resampled_step_minute["isactivebout"] != resampled_step_minute["isactivebout"].shift()).cumsum().rename(columns={"isactivebout": "group_idx"})

    # combine resampled_acc_minute and group column
    resampled_step_minute = pd.concat([resampled_step_minute, group], axis=1)
    
    # drop rows where "isactivebout" column is missing and reset the index
    resampled_step_minute.dropna(subset=["isactivebout"], inplace=True)
    resampled_step_minute.reset_index(inplace=True)
    resampled_step_minute.loc[:, "local_date"] = resampled_step_minute["local_date_time"].dt.date

    # duration column contains the number of minutes (rows) of active and sedentary bout
    bouts = resampled_step_minute.groupby(["isactivebout", "group_idx", "local_date"]).count().rename(columns={"local_date_time": "duration"}).reset_index()
    bouts["duration"] = bouts["duration"] * time_interval

    return bouts

def statsFeatures(step_data, day_segment, features_to_compute, features_type, step_features):
    if features_type == "allsteps":
        col_name = "steps"
    elif features_type == "durationsedentarybout" or features_type == "durationactivebout":
        col_name = "duration"
    else:
        raise ValueError("features_type can only be one of ['allsteps', 'durationsedentarybout', 'durationactivebout'].")

    if "count" + features_type.replace("duration", "episode") in features_to_compute:
        step_features["step_" + day_segment + "_count" + features_type.replace("duration", "episode")] = step_data.groupby(["local_date"])[col_name].count()
    if "sum" + features_type in features_to_compute:
        step_features["step_" + day_segment + "_sum" + features_type] = step_data.groupby(["local_date"])[col_name].sum()
    if "max" + features_type in features_to_compute:
        step_features["step_" + day_segment + "_max" + features_type] = step_data.groupby(["local_date"])[col_name].max()
    if "min" + features_type in features_to_compute:
        step_features["step_" + day_segment + "_min" + features_type] = step_data.groupby(["local_date"])[col_name].min()
    if "avg" + features_type in features_to_compute:
        step_features["step_" + day_segment + "_avg" + features_type] = step_data.groupby(["local_date"])[col_name].mean()
    if "median" + features_type in features_to_compute:
        step_features["step_" + day_segment + "_median" + features_type] = step_data.groupby(["local_date"])[col_name].median()
    if "std" + features_type in features_to_compute:
        step_features["step_" + day_segment + "_std" + features_type] = step_data.groupby(["local_date"])[col_name].std()

    return step_features

def base_fitbit_step_features(step_data, day_segment, requested_features, threshold_active_bout, include_zero_step_rows):
    requested_features_allsteps = requested_features["features_all_steps"]
    requested_features_sedentarybout = requested_features["features_sedentary_bout"]
    requested_features_activebout = requested_features["features_active_bout"]

    # name of the features this function can compute
    base_features_allsteps = ["sumallsteps", "maxallsteps", "minallsteps", "avgallsteps", "stdallsteps"]
    base_features_sedentarybout = ["countepisodesedentarybout", "sumdurationsedentarybout", "maxdurationsedentarybout", "mindurationsedentarybout", "avgdurationsedentarybout", "stddurationsedentarybout"]
    base_features_activebout = ["countepisodeactivebout", "sumdurationactivebout", "maxdurationactivebout", "mindurationactivebout", "avgdurationactivebout", "stddurationactivebout"]
    # the subset of requested features this function can compute
    features_to_compute_allsteps = list(set(requested_features_allsteps) & set(base_features_allsteps))
    features_to_compute_sedentarybout = list(set(requested_features_sedentarybout) & set(base_features_sedentarybout))
    features_to_compute_activebout = list(set(requested_features_activebout) & set(base_features_activebout))

    features_to_compute = features_to_compute_allsteps + features_to_compute_sedentarybout + features_to_compute_activebout

    step_features = pd.DataFrame(columns=["local_date"] + ["step_" + day_segment + "_" + x for x in features_to_compute])
    if not step_data.empty:
        if day_segment != "daily":
            step_data =step_data[step_data["local_day_segment"] == day_segment]
        
        if not step_data.empty:
            step_features = pd.DataFrame()

            # statistics features of step count
            step_features = statsFeatures(step_data, day_segment, features_to_compute_allsteps, "allsteps", step_features)

            # calculate time interval between two records in minutes
            time_interval = step_data["local_date_time"].diff().min().total_seconds() / 60

            # sedentary bout: less than THRESHOLD_ACTIVE_BOUT (default: 10) steps in a minute
            # active bout: greater or equal to THRESHOLD_ACTIVE_BOUT (default: 10) steps in a minute
            isactivebout = np.where(step_data["steps"] < int(threshold_active_bout) * time_interval, 0, 1)
            step_data = step_data.assign(isactivebout = isactivebout)

            bouts = getBouts(step_data, time_interval)

            # statistics features of sedentary bout
            sedentary_bout = bouts[bouts["isactivebout"] == 0]
            step_features = statsFeatures(sedentary_bout, day_segment, features_to_compute_sedentarybout, "durationsedentarybout", step_features)

            # statistics features of active bout
            active_bout = bouts[bouts["isactivebout"] == 1]
            step_features = statsFeatures(active_bout, day_segment, features_to_compute_activebout, "durationactivebout", step_features)

            # exclude data when the total step count is ZERO during the whole epoch
            if not include_zero_step_rows:
                step_features["sumallsteps_aux"] = step_data.groupby(["local_date"])["steps"].sum()
                step_features = step_features.query("sumallsteps_aux != 0")
                del step_features["sumallsteps_aux"]

            step_features = step_features.reset_index()

    return step_features
