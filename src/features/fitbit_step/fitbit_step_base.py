import pandas as pd
import numpy as np
import datetime as dt
from features_utils import splitOvernightEpisodes, splitMultiSegmentEpisodes

def base_fitbit_step_features(step_data, day_segment, requested_features, threshold_active_bout, include_zero_step_rows):
    requested_features_allsteps = requested_features["features_all_steps"]
    requested_features_sedentarybout = requested_features["features_sedentary_bout"]
    requested_features_activebout = requested_features["features_active_bout"]

    # name of the features this function can compute
    base_features_allsteps = ["sumallsteps", "maxallsteps", "minallsteps", "avgallsteps", "stdallsteps"]
    base_features_sedentarybout = ["countsedentarybout", "maxdurationsedentarybout", "mindurationsedentarybout", "avgdurationsedentarybout", "stddurationsedentarybout", "sumdurationsedentarybout"]
    base_features_activebout = ["countactivebout", "maxdurationactivebout", "mindurationactivebout", "avgdurationactivebout", "stddurationactivebout"]
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

            resampled_data = step_data.set_index(step_data.local_date_time)
            resampled_data.index.names = ["datetime"]

            # Replace the first element of time_diff_minutes with its second element
            resampled_data["time_diff_minutes"] = resampled_data["local_date_time"].diff().fillna(resampled_data["local_date_time"].diff()[1]).dt.total_seconds().div(60).astype(int)

            # Sedentary Bout when you have less than 10 steps in a minute
            # Active Bout when you have greater or equal to 10 steps in a minute
            resampled_data["active_sedentary"] = np.where(resampled_data["steps"] < int(threshold_active_bout) * resampled_data["time_diff_minutes"],"sedentary","active")

            # Time Calculations of sedentary/active bouts:
            resampled_data["active_sedentary_groups"] = (resampled_data.active_sedentary != resampled_data.active_sedentary.shift()).cumsum().values

            # Get the total minutes for each episode
            minutes_per_episode = resampled_data.groupby(["local_date","active_sedentary","active_sedentary_groups"])["time_diff_minutes"].sum()
            
            # Get Stats for all episodes in terms of minutes
            stats_per_episode = minutes_per_episode.groupby(["local_date", "active_sedentary"]).agg([max, min, np.mean, np.std, np.sum])
            mux = pd.MultiIndex.from_product([stats_per_episode.index.levels[0], stats_per_episode.index.levels[1]], names=["local_date", "active_sedentary"])
            stats_per_episode = stats_per_episode.reindex(mux, fill_value=None).reset_index()
            stats_per_episode.set_index("local_date", inplace = True)
            
            # Descriptive Statistics Features:
            if "sumallsteps" in features_to_compute_allsteps:
                step_features["step_" + str(day_segment) + "_sumallsteps"] = resampled_data["steps"].resample("D").sum()
            if "maxallsteps" in features_to_compute_allsteps:
                step_features["step_" + str(day_segment) + "_maxallsteps"] = resampled_data["steps"].resample("D").max()
            if "minallsteps" in features_to_compute_allsteps:
                step_features["step_" + str(day_segment) + "_minallsteps"] = resampled_data["steps"].resample("D").min()
            if "avgallsteps" in features_to_compute_allsteps:
                step_features["step_" + str(day_segment) + "_avgallsteps"] = resampled_data["steps"].resample("D").mean()
            if "stdallsteps" in features_to_compute_allsteps:
                step_features["step_" + str(day_segment) + "_stdallsteps"] = resampled_data["steps"].resample("D").std()

            if "countsedentarybout" in features_to_compute_sedentarybout:
                step_features["step_" + str(day_segment) + "_countsedentarybout"] = resampled_data[resampled_data["active_sedentary"] == "sedentary"]["active_sedentary_groups"].resample("D").nunique()
            if "countactivebout" in features_to_compute_activebout:
                step_features["step_" + str(day_segment) + "_countactivebout"] = resampled_data[resampled_data["active_sedentary"] == "active"]["active_sedentary_groups"].resample("D").nunique()
            if "maxdurationsedentarybout" in features_to_compute_sedentarybout:
                step_features["step_" + str(day_segment) + "_maxdurationsedentarybout"] = stats_per_episode[stats_per_episode["active_sedentary"]=="sedentary"]["max"]
            if "mindurationsedentarybout" in features_to_compute_sedentarybout:
                step_features["step_" + str(day_segment) + "_mindurationsedentarybout"] = stats_per_episode[stats_per_episode["active_sedentary"]=="sedentary"]["min"]
            if "avgdurationsedentarybout" in features_to_compute_sedentarybout:
                step_features["step_" + str(day_segment) + "_avgdurationsedentarybout"] = stats_per_episode[stats_per_episode["active_sedentary"]=="sedentary"]["mean"]
            if "stddurationsedentarybout" in features_to_compute_sedentarybout:
                step_features["step_" + str(day_segment) + "_stddurationsedentarybout"] = stats_per_episode[stats_per_episode["active_sedentary"]=="sedentary"]["std"]
            if "sumdurationsedentarybout" in features_to_compute_sedentarybout:
                step_features["step_" + str(day_segment) + "_sumdurationsedentarybout"] = stats_per_episode[stats_per_episode["active_sedentary"]=="sedentary"]["sum"]
            if "maxdurationactivebout" in features_to_compute_activebout:
                step_features["step_" + str(day_segment) + "_maxdurationactivebout"] = stats_per_episode[stats_per_episode["active_sedentary"]== "active"]["max"]
            if "mindurationactivebout" in features_to_compute_activebout:
                step_features["step_" + str(day_segment) + "_mindurationactivebout"] = stats_per_episode[stats_per_episode["active_sedentary"]== "active"]["min"]
            if "avgdurationactivebout" in features_to_compute_activebout:
                step_features["step_" + str(day_segment) + "_avgdurationactivebout"] = stats_per_episode[stats_per_episode["active_sedentary"]== "active"]["mean"]
            if "stddurationactivebout" in features_to_compute_activebout:
                step_features["step_" + str(day_segment) + "_stddurationactivebout"] = stats_per_episode[stats_per_episode["active_sedentary"]== "active"]["std"]
            
            #Exclude data when the total step count is ZERO during the whole epoch
            if not include_zero_step_rows:
                step_features["sumallsteps_aux"] = resampled_data["steps"].resample("D").sum()
                step_features = step_features.query("sumallsteps_aux != 0")
                del step_features["sumallsteps_aux"]

            step_features.index.names = ["local_date"]
            step_features = step_features.reset_index()

    return step_features
