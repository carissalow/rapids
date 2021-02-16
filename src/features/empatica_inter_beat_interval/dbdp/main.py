import pandas as pd
from scipy.stats import entropy


def statsFeatures(ibi_data, features, ibi_features):
    col_name = "inter_beat_interval"
    if "sumibi" in features:
        ibi_features["sumibi"] = ibi_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].sum()
    if "maxibi" in features:
        ibi_features["maxibi"] = ibi_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].max()
    if "minibi" in features:
        ibi_features["minibi"] = ibi_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].min()
    if "avgibi" in features:
        ibi_features["avgibi"] = ibi_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].mean()
    if "medianibi" in features:
        ibi_features["medianibi"] = ibi_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].median()
    if "modeibi" in features:
        ibi_features["modeibi"] = ibi_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].agg(lambda x: pd.Series.mode(x)[0])
    if "stdibi" in features:
        ibi_features["stdibi"] = ibi_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].std()
    if "diffmaxmodeibi" in features:
        ibi_features["diffmaxmodeibi"] = ibi_data[["local_segment", col_name]].groupby(["local_segment"])[
                                                  col_name].max() - \
                                              ibi_data[["local_segment", col_name]].groupby(["local_segment"])[
                                                  col_name].agg(lambda x: pd.Series.mode(x)[0])
    if "diffminmodeibi" in features:
        ibi_features["diffminmodeibi"] = ibi_data[["local_segment", col_name]].groupby(["local_segment"])[
                                                  col_name].agg(lambda x: pd.Series.mode(x)[0]) - \
                                              ibi_data[["local_segment", col_name]].groupby(["local_segment"])[
                                                  col_name].min()
    if "entropyibi" in features:
        ibi_features["entropyibi"] = ibi_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].agg(entropy)

    return ibi_features


def extractIBIFeaturesFromIntradayData(ibi_intraday_data, features, time_segment, filter_data_by_segment):
    ibi_intraday_features = pd.DataFrame(columns=["local_segment"] + features)
    if not ibi_intraday_data.empty:
        ibi_intraday_data = filter_data_by_segment(ibi_intraday_data, time_segment)

        if not ibi_intraday_data.empty:
            ibi_intraday_features = pd.DataFrame()

            # get stats of ibi
            ibi_intraday_features = statsFeatures(ibi_intraday_data, features, ibi_intraday_features)

            ibi_intraday_features.reset_index(inplace=True)

    return ibi_intraday_features


def dbdp_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):
    ibi_intraday_data = pd.read_csv(sensor_data_files["sensor_data"])

    requested_intraday_features = provider["FEATURES"]
    # name of the features this function can compute
    base_intraday_features_names = ["maxibi", "minibi", "avgibi", "medianibi", "modeibi", "stdibi", "diffmaxmodeibi",
                                    "diffminmodeibi", "entropyibi"]
    # the subset of requested features this function can compute
    intraday_features_to_compute = list(set(requested_intraday_features) & set(base_intraday_features_names))

    # extract features from intraday data
    ibi_intraday_features = extractIBIFeaturesFromIntradayData(ibi_intraday_data,
                                                                    intraday_features_to_compute, time_segment,
                                                                    filter_data_by_segment)

    return ibi_intraday_features