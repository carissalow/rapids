import pandas as pd
import numpy as np
import time
from fitbit_step.fitbit_step_base import base_fitbit_step_features

def isInvalidTime(str_time):
    try:
        time.strptime(str_time, '%H:%M')
        return False
    except ValueError:
        return True

def isInMainSleep(local_date_time, sleep):
    # sleep_period_container = sleep.query("local_start_date_time <= @local_date_time <= local_end_date_time")
    sleep_period_container = sleep[(sleep["local_start_date_time"] <= local_date_time) & (local_date_time <= sleep["local_end_date_time"])]
    if sleep_period_container.shape[0] >= 1:
        return True
    else:
        return False

def getStepsOutsideFitbitMainSleep(sleep, steps):
    steps['inMainSleep'] = steps.apply(lambda row : isInMainSleep(row['local_date_time'], sleep), axis = 1) 
    return steps[steps['inMainSleep'] == False]


def getStepsOutsideFixedMainSleep(sleepStart, sleepEnd, steps):
    steps = steps.set_index('local_date_time')
    steps['inMainSleep'] = False
    steps.loc[steps.between_time(sleepStart, sleepEnd).index, 'inMainSleep'] = True
    steps.reset_index(level=0, inplace=True)
    return steps[steps['inMainSleep'] == False]

step_data = pd.read_csv(snakemake.input["step_data"], parse_dates=["local_date_time", "local_date"])
day_segment = snakemake.params["day_segment"]
threshold_active_bout = snakemake.params["threshold_active_bout"]
include_zero_step_rows = snakemake.params["include_zero_step_rows"]
exclude_sleep = snakemake.params["exclude_sleep"]
exclude_sleep_type = snakemake.params["exclude_sleep_type"]
exclude_sleep_fixed_start = snakemake.params["exclude_sleep_fixed_start"]
exclude_sleep_fixed_end = snakemake.params["exclude_sleep_fixed_end"]

step_features = pd.DataFrame(columns=["local_date"])
requested_features = {}
requested_features["features_all_steps"] = snakemake.params["features_all_steps"]
requested_features["features_sedentary_bout"] = [feature + "sedentarybout" for feature in snakemake.params["features_sedentary_bout"]]
requested_features["features_active_bout"] = [feature + "activebout" for feature in snakemake.params["features_active_bout"]]

if exclude_sleep == True:
    if exclude_sleep_type == "FIXED":
        if isInvalidTime(exclude_sleep_fixed_start):
            raise ValueError("Your fixed start time has an invalid format in your config.yml file")
        if isInvalidTime(exclude_sleep_fixed_end):
            raise ValueError("Your fixed end time has an invalid format in your config.yml file")
        step_data = getStepsOutsideFixedMainSleep(exclude_sleep_fixed_start, exclude_sleep_fixed_end, step_data)
    elif exclude_sleep_type == "FITBIT_BASED":
        sleep_data = pd.read_csv(snakemake.input["sleep_data"], parse_dates=["local_start_date_time", "local_end_date_time"])
        step_data = getStepsOutsideFitbitMainSleep(sleep_data, step_data)
    else:
        raise ValueError("We only support FIXED or FITBIT_BASED to filter step data based on sleep data. You typed " + exclude_sleep_type + ", Check your config.yaml file for typos")

step_features = step_features.merge(base_fitbit_step_features(step_data, day_segment, requested_features, threshold_active_bout, include_zero_step_rows), on="local_date", how="outer")


assert np.sum([len(x) for x in requested_features.values()]) + 1 == step_features.shape[1], "The number of features in the output dataframe (=" + str(step_features.shape[1]) + ") does not match the expected value (=" + str(np.sum([len(x) for x in requested_features.values()])) + " + 1). Verify your fitbit step feature extraction functions"

step_features.to_csv(snakemake.output[0], index=False)