import pandas as pd
from scipy.stats import entropy


def statsFeatures(temperature_data, features, temperature_features):
    col_name = "temperature"
    if "sumtemp" in features:
        temperature_features["sumtemp"] = temperature_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].sum()
    if "maxtemp" in features:
        temperature_features["maxtemp"] = temperature_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].max()
    if "mintemp" in features:
        temperature_features["mintemp"] = temperature_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].min()
    if "avgtemp" in features:
        temperature_features["avgtemp"] = temperature_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].mean()
    if "mediantemp" in features:
        temperature_features["mediantemp"] = temperature_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].median()
    if "modetemp" in features:
        temperature_features["modetemp"] = temperature_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].agg(lambda x: pd.Series.mode(x)[0])
    if "stdtemp" in features:
        temperature_features["stdtemp"] = temperature_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].std()
    if "diffmaxmodetemp" in features:
        temperature_features["diffmaxmodetemp"] = temperature_data[["local_segment", col_name]].groupby(["local_segment"])[
                                                  col_name].max() - \
                                              temperature_data[["local_segment", col_name]].groupby(["local_segment"])[
                                                  col_name].agg(lambda x: pd.Series.mode(x)[0])
    if "diffminmodetemp" in features:
        temperature_features["diffminmodetemp"] = temperature_data[["local_segment", col_name]].groupby(["local_segment"])[
                                                  col_name].agg(lambda x: pd.Series.mode(x)[0]) - \
                                              temperature_data[["local_segment", col_name]].groupby(["local_segment"])[
                                                  col_name].min()
    if "entropytemp" in features:
        temperature_features["entropytemp"] = temperature_data[["local_segment", col_name]].groupby(["local_segment"])[
            col_name].agg(entropy)

    return temperature_features


def extractTempFeaturesFromIntradayData(temperature_intraday_data, features, time_segment, filter_data_by_segment):
    temperature_intraday_features = pd.DataFrame(columns=["local_segment"] + features)
    if not temperature_intraday_data.empty:
        temperature_intraday_data = filter_data_by_segment(temperature_intraday_data, time_segment)

        if not temperature_intraday_data.empty:
            temperature_intraday_features = pd.DataFrame()

            # get stats of temperature
            temperature_intraday_features = statsFeatures(temperature_intraday_data, features, temperature_intraday_features)

            temperature_intraday_features.reset_index(inplace=True)

    return temperature_intraday_features


def dbdp_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):
    temperature_intraday_data = pd.read_csv(sensor_data_files["sensor_data"])

    requested_intraday_features = provider["FEATURES"]
    # name of the features this function can compute
    base_intraday_features_names = ["maxtemp", "mintemp", "avgtemp", "mediantemp", "modetemp", "stdtemp", "diffmaxmodetemp",
                                    "diffminmodetemp", "entropytemp"]
    # the subset of requested features this function can compute
    intraday_features_to_compute = list(set(requested_intraday_features) & set(base_intraday_features_names))

    # extract features from intraday data
    temperature_intraday_features = extractTempFeaturesFromIntradayData(temperature_intraday_data,
                                                                    intraday_features_to_compute, time_segment,
                                                                    filter_data_by_segment)

    return temperature_intraday_features