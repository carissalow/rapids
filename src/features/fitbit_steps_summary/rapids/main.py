import pandas as pd
import numpy as np

def statsFeatures(steps_data, features_to_compute, features_type, steps_features):
    if features_type == "steps" or features_type == "sumsteps":
        col_name = "steps"
    elif features_type == "durationsedentarybout" or features_type == "durationactivebout":
        col_name = "duration"
    elif features_type == "volatilitysteps":
        col_name = "volatilitysteps"
    else:
        raise ValueError("features_type can only be one of ['steps', 'sumsteps', 'durationsedentarybout', 'durationactivebout', 'volatilitysteps'].")

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
    if "annualized" + features_type in features_to_compute:
        steps_features["annualized" + features_type] = steps_data.groupby(["local_segment"])[col_name].var(ddof = 1)
    return steps_features

def extractStepsVolatility(steps_data):
    #Set index to be datetime
    steps_data.index = pd.DatetimeIndex(steps_data["local_date"])

    #Add in the missing dates between date range
    date_range = pd.date_range(start = steps_data.index.min(), end = steps_data.index.max(), freq="D")
    steps_data = steps_data.reindex(date_range, fill_value = np.nan)

    #Create the denominator for volatility function
    steps_data.loc[steps_data["steps"] == 0, "steps"] = np.nan
    steps_data["last_steps"] = steps_data["steps"].shift(periods = 1, fill_value = np.nan)

    steps_data['volatilitysteps'] = np.log(steps_data["steps"] / steps_data["last_steps"])

    steps_data.dropna(subset = ["local_segment"], inplace = True)
    steps_data.reset_index(drop = True, inplace = True)

    return steps_data

def extractStepsFeaturesFromSummaryData(steps_summary_data, summary_features_to_compute):
    steps_summary_features = pd.DataFrame()

    # statistics features of daily steps count
    steps_summary_features = statsFeatures(steps_summary_data, summary_features_to_compute, "sumsteps", steps_summary_features)
    steps_summary_features = statsFeatures(steps_summary_data, summary_features_to_compute, "volatilitysteps", steps_summary_features)

    steps_summary_features.reset_index(inplace=True)
    
    return steps_summary_features



def rapids_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):

    steps_summary_data = pd.read_csv(sensor_data_files["sensor_data"])
    requested_summary_features = provider["FEATURES"]

    # name of the features this function can compute
    base_summary_features = ["maxsumsteps", "minsumsteps", "avgsumsteps", "mediansumsteps", "stdsumsteps", 
                             "maxvolatilitysteps", "minvolatilitysteps", "avgvolatilitysteps", "medianvolatilitysteps", "stdvolatilitysteps", "annualizedvolatilitysteps"]

    # the subset of requested features this function can compute
    summary_features_to_compute = list(set(requested_summary_features) & set(base_summary_features))

    # extract features from summary data
    steps_summary_features = pd.DataFrame(columns=["local_segment"] + summary_features_to_compute)

    if not steps_summary_data.empty:
        steps_summary_data = filter_data_by_segment(steps_summary_data, time_segment)

        if not steps_summary_data.empty:

            # only keep the segments start at 00:00:00 and end at 23:59:59
            datetime_start_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2} 00:00:00"
            datetime_end_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2} 23:59:59"

            segment_regex = "{}#{},{}".format(time_segment, datetime_start_regex, datetime_end_regex)
            steps_summary_data = steps_summary_data[steps_summary_data["local_segment"].str.match(segment_regex)]

            if not steps_summary_data.empty:

                steps_summary_data = extractStepsVolatility(steps_summary_data)
                steps_summary_features = extractStepsFeaturesFromSummaryData(steps_summary_data, summary_features_to_compute)
    
    return steps_summary_features

