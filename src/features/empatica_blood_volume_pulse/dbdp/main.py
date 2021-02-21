import pandas as pd
from scipy.stats import entropy


def statsFeatures(bvp_data, features, bvp_features):
    col_name = "blood_volume_pulse"
    if "sumbvp" in features:
        bvp_features["sumbvp"] = bvp_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].sum()
    if "maxbvp" in features:
        bvp_features["maxbvp"] = bvp_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].max()
    if "minbvp" in features:
        bvp_features["minbvp"] = bvp_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].min()
    if "avgbvp" in features:
        bvp_features["avgbvp"] = bvp_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].mean()
    if "medianbvp" in features:
        bvp_features["medianbvp"] = bvp_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].median()
    if "modebvp" in features:
        bvp_features["modebvp"] = bvp_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].agg(lambda x: pd.Series.mode(x)[0])
    if "stdbvp" in features:
        bvp_features["stdbvp"] = bvp_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].std()
    if "diffmaxmodebvp" in features:
        bvp_features["diffmaxmodebvp"] = bvp_data[["local_segment", col_name]].groupby(["local_segment"])[
                                                  col_name].max() - \
                                              bvp_data[["local_segment", col_name]].groupby(["local_segment"])[
                                                  col_name].agg(lambda x: pd.Series.mode(x)[0])
    if "diffminmodebvp" in features:
        bvp_features["diffminmodebvp"] = bvp_data[["local_segment", col_name]].groupby(["local_segment"])[
                                                  col_name].agg(lambda x: pd.Series.mode(x)[0]) - \
                                              bvp_data[["local_segment", col_name]].groupby(["local_segment"])[
                                                  col_name].min()
    if "entropybvp" in features:
        bvp_features["entropybvp"] = bvp_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].agg(entropy)

    return bvp_features


def extractBVPFeaturesFromIntradayData(bvp_intraday_data, features, time_segment, filter_data_by_segment):
    bvp_intraday_features = pd.DataFrame(columns=["local_segment"] + features)
    if not bvp_intraday_data.empty:
        bvp_intraday_data = filter_data_by_segment(bvp_intraday_data, time_segment)

        if not bvp_intraday_data.empty:
            bvp_intraday_features = pd.DataFrame()

            # get stats of bvp
            bvp_intraday_features = statsFeatures(bvp_intraday_data, features, bvp_intraday_features)

            bvp_intraday_features.reset_index(inplace=True)

    return bvp_intraday_features


def dbdp_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):
    bvp_intraday_data = pd.read_csv(sensor_data_files["sensor_data"])

    requested_intraday_features = provider["FEATURES"]
    # name of the features this function can compute
    base_intraday_features_names = ["maxbvp", "minbvp", "avgbvp", "medianbvp", "modebvp", "stdbvp", "diffmaxmodebvp",
                                    "diffminmodebvp", "entropybvp"]
    # the subset of requested features this function can compute
    intraday_features_to_compute = list(set(requested_intraday_features) & set(base_intraday_features_names))

    # extract features from intraday data
    bvp_intraday_features = extractBVPFeaturesFromIntradayData(bvp_intraday_data,
                                                                    intraday_features_to_compute, time_segment,
                                                                    filter_data_by_segment)

    return bvp_intraday_features