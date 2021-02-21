import pandas as pd
from scipy.stats import entropy


def statsFeatures(eda_data, features, eda_features):
    col_name = "electrodermal_activity"
    if "sumeda" in features:
        eda_features["sumeda"] = eda_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].sum()
    if "maxeda" in features:
        eda_features["maxeda"] = eda_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].max()
    if "mineda" in features:
        eda_features["mineda"] = eda_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].min()
    if "avgeda" in features:
        eda_features["avgeda"] = eda_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].mean()
    if "medianeda" in features:
        eda_features["medianeda"] = eda_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].median()
    if "modeeda" in features:
        eda_features["modeeda"] = eda_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].agg(lambda x: pd.Series.mode(x)[0])
    if "stdeda" in features:
        eda_features["stdeda"] = eda_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].std()
    if "diffmaxmodeeda" in features:
        eda_features["diffmaxmodeeda"] = eda_data[["local_segment", col_name]].groupby(["local_segment"])[
                                                  col_name].max() - \
                                              eda_data[["local_segment", col_name]].groupby(["local_segment"])[
                                                  col_name].agg(lambda x: pd.Series.mode(x)[0])
    if "diffminmodeeda" in features:
        eda_features["diffminmodeeda"] = eda_data[["local_segment", col_name]].groupby(["local_segment"])[
                                                  col_name].agg(lambda x: pd.Series.mode(x)[0]) - \
                                              eda_data[["local_segment", col_name]].groupby(["local_segment"])[
                                                  col_name].min()
    if "entropyeda" in features:
        eda_features["entropyeda"] = eda_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].agg(entropy)

    return eda_features


def extractEDAFeaturesFromIntradayData(eda_intraday_data, features, time_segment, filter_data_by_segment):
    eda_intraday_features = pd.DataFrame(columns=["local_segment"] + features)
    if not eda_intraday_data.empty:
        eda_intraday_data = filter_data_by_segment(eda_intraday_data, time_segment)

        if not eda_intraday_data.empty:
            eda_intraday_features = pd.DataFrame()

            # get stats of eda
            eda_intraday_features = statsFeatures(eda_intraday_data, features, eda_intraday_features)

            eda_intraday_features.reset_index(inplace=True)

    return eda_intraday_features


def dbdp_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):
    eda_intraday_data = pd.read_csv(sensor_data_files["sensor_data"])

    requested_intraday_features = provider["FEATURES"]
    # name of the features this function can compute
    base_intraday_features_names = ["maxeda", "mineda", "avgeda", "medianeda", "modeeda", "stdeda", "diffmaxmodeeda",
                                    "diffminmodeeda", "entropyeda"]
    # the subset of requested features this function can compute
    intraday_features_to_compute = list(set(requested_intraday_features) & set(base_intraday_features_names))

    # extract features from intraday data
    eda_intraday_features = extractEDAFeaturesFromIntradayData(eda_intraday_data,
                                                                    intraday_features_to_compute, time_segment,
                                                                    filter_data_by_segment)

    return eda_intraday_features