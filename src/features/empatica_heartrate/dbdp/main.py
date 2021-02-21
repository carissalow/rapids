import pandas as pd
from scipy.stats import entropy


def statsFeatures(heartrate_data, features, heartrate_features):
    col_name = "heartrate"
    if "sumhr" in features:
        heartrate_features["sumhr"] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].sum()
    if "maxhr" in features:
        heartrate_features["maxhr"] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].max()
    if "minhr" in features:
        heartrate_features["minhr"] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].min()
    if "avghr" in features:
        heartrate_features["avghr"] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].mean()
    if "medianhr" in features:
        heartrate_features["medianhr"] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].median()
    if "modehr" in features:
        heartrate_features["modehr"] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].agg(lambda x: pd.Series.mode(x)[0])
    if "stdhr" in features:
        heartrate_features["stdhr"] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].std()
    if "diffmaxmodehr" in features:
        heartrate_features["diffmaxmodehr"] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[
                                                  col_name].max() - \
                                              heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[
                                                  col_name].agg(lambda x: pd.Series.mode(x)[0])
    if "diffminmodehr" in features:
        heartrate_features["diffminmodehr"] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[
                                                  col_name].agg(lambda x: pd.Series.mode(x)[0]) - \
                                              heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[
                                                  col_name].min()
    if "entropyhr" in features:
        heartrate_features["entropyhr"] = heartrate_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].agg(entropy)

    return heartrate_features


def extractHRFeaturesFromIntradayData(heartrate_intraday_data, features, time_segment, filter_data_by_segment):
    heartrate_intraday_features = pd.DataFrame(columns=["local_segment"] + features)
    if not heartrate_intraday_data.empty:
        heartrate_intraday_data = filter_data_by_segment(heartrate_intraday_data, time_segment)

        if not heartrate_intraday_data.empty:
            heartrate_intraday_features = pd.DataFrame()

            # get stats of heartrate
            heartrate_intraday_features = statsFeatures(heartrate_intraday_data, features, heartrate_intraday_features)

            heartrate_intraday_features.reset_index(inplace=True)

    return heartrate_intraday_features


def dbdp_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):
    heartrate_intraday_data = pd.read_csv(sensor_data_files["sensor_data"])

    requested_intraday_features = provider["FEATURES"]
    # name of the features this function can compute
    base_intraday_features_names = ["maxhr", "minhr", "avghr", "medianhr", "modehr", "stdhr", "diffmaxmodehr",
                                    "diffminmodehr", "entropyhr"]
    # the subset of requested features this function can compute
    intraday_features_to_compute = list(set(requested_intraday_features) & set(base_intraday_features_names))

    # extract features from intraday data
    heartrate_intraday_features = extractHRFeaturesFromIntradayData(heartrate_intraday_data,
                                                                    intraday_features_to_compute, time_segment,
                                                                    filter_data_by_segment)

    return heartrate_intraday_features
