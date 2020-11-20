import pandas as pd
from scipy.stats import entropy

def statsFeatures(heartrate_data, features, features_type, heartrate_features):

    if features_type == "hr":
        col_name = "heartrate"
    elif features_type == "restinghr":
        col_name = "heartrate_daily_restinghr"
    elif features_type == "caloriesoutofrange":
        col_name = "heartrate_daily_caloriesoutofrange"
    elif features_type == "caloriesfatburn":
        col_name = "heartrate_daily_caloriesfatburn"
    elif features_type == "caloriescardio":
        col_name = "heartrate_daily_caloriescardio"
    elif features_type == "caloriespeak":
        col_name = "heartrate_daily_caloriespeak"
    else:
        raise ValueError("features_type can only be one of ['hr', 'restinghr', 'caloriesoutofrange', 'caloriesfatburn', 'caloriescardio', 'caloriespeak'].")

    if "summarysum" + features_type in features:
        heartrate_features["heartrate_rapids_summarysum" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].sum()
    if "summarymax" + features_type in features:
        heartrate_features["heartrate_rapids_summarymax" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].max()
    if "summarymin" + features_type in features:
        heartrate_features["heartrate_rapids_summarymin" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].min()
    if "summaryavg" + features_type in features:
        heartrate_features["heartrate_rapids_summaryavg" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].mean()
    if "summarymedian" + features_type in features:
        heartrate_features["heartrate_rapids_summarymedian" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].median()
    if "summarymode" + features_type in features:
        heartrate_features["heartrate_rapids_summarymode" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].agg(lambda x: pd.Series.mode(x)[0])
    if "summarystd" + features_type in features:
        heartrate_features["heartrate_rapids_summarystd" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].std()
    if "summarydiffmaxmode" + features_type in features:
        heartrate_features["heartrate_rapids_summarydiffmaxmode" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].max() - heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].agg(lambda x: pd.Series.mode(x)[0])
    if "summarydiffminmode" + features_type in features:
        heartrate_features["heartrate_rapids_summarydiffminmode" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].agg(lambda x: pd.Series.mode(x)[0]) - heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].min()
    if "summaryentropy" + features_type in features:
        heartrate_features["heartrate_rapids_summaryentropy" + features_type] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[col_name].agg(entropy)
    
    return heartrate_features

def extractHRFeaturesFromSummaryData(heartrate_summary_data, summary_features):
    heartrate_summary_features = pd.DataFrame()

    # get stats of resting heartrate
    heartrate_summary_features = statsFeatures(heartrate_summary_data, summary_features, "restinghr", heartrate_summary_features)

    # get stats of calories features
    # calories features might be inaccurate: they depend on users' fitbit profile (weight, height, etc.)
    heartrate_summary_features = statsFeatures(heartrate_summary_data, summary_features, "caloriesoutofrange", heartrate_summary_features)
    heartrate_summary_features = statsFeatures(heartrate_summary_data, summary_features, "caloriesfatburn", heartrate_summary_features)
    heartrate_summary_features = statsFeatures(heartrate_summary_data, summary_features, "caloriescardio", heartrate_summary_features)
    heartrate_summary_features = statsFeatures(heartrate_summary_data, summary_features, "caloriespeak", heartrate_summary_features)

    heartrate_summary_features.reset_index(inplace=True)
    
    return heartrate_summary_features


def rapids_features(sensor_data_files, day_segment, provider, filter_data_by_segment, *args, **kwargs):

    heartrate_summary_data = pd.read_csv(sensor_data_files["sensor_data"])

    requested_summary_features = ["summary" + x for x in provider["FEATURES"]]
    # name of the features this function can compute
    base_summary_features_names = ["summarymaxrestinghr", "summaryminrestinghr", "summaryavgrestinghr", "summarymedianrestinghr", "summarymoderestinghr", "summarystdrestinghr", "summarydiffmaxmoderestinghr", "summarydiffminmoderestinghr", "summaryentropyrestinghr", "summarysumcaloriesoutofrange", "summarymaxcaloriesoutofrange", "summarymincaloriesoutofrange", "summaryavgcaloriesoutofrange", "summarymediancaloriesoutofrange", "summarystdcaloriesoutofrange", "summaryentropycaloriesoutofrange", "summarysumcaloriesfatburn", "summarymaxcaloriesfatburn", "summarymincaloriesfatburn", "summaryavgcaloriesfatburn", "summarymediancaloriesfatburn", "summarystdcaloriesfatburn", "summaryentropycaloriesfatburn", "summarysumcaloriescardio", "summarymaxcaloriescardio", "summarymincaloriescardio", "summaryavgcaloriescardio", "summarymediancaloriescardio", "summarystdcaloriescardio", "summaryentropycaloriescardio", "summarysumcaloriespeak", "summarymaxcaloriespeak", "summarymincaloriespeak", "summaryavgcaloriespeak", "summarymediancaloriespeak", "summarystdcaloriespeak", "summaryentropycaloriespeak"]
    # the subset of requested features this function can compute
    summary_features_to_compute = list(set(requested_summary_features) & set(base_summary_features_names))

    # extract features from summary data
    heartrate_summary_features = pd.DataFrame(columns=["local_segment"] + ["heartrate_rapids_" + x for x in summary_features_to_compute])
    if not heartrate_summary_data.empty:
        heartrate_summary_data = filter_data_by_segment(heartrate_summary_data, day_segment)

        if not heartrate_summary_data.empty:
            # only keep the segments start at 00:00:00 and end at 23:59:59
            datetime_start_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2} 00:00:00"
            datetime_end_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2} 23:59:59"

            segment_regex = "{}#{},{}".format(day_segment, datetime_start_regex, datetime_end_regex)
            heartrate_summary_data = heartrate_summary_data[heartrate_summary_data["local_segment"].str.match(segment_regex)]

            if not heartrate_summary_data.empty:
                heartrate_summary_features = extractHRFeaturesFromSummaryData(heartrate_summary_data, summary_features_to_compute)
    
    return heartrate_summary_features
