from barnett_library import *
from statistics import mode
import warnings

def barnett_daily_features(snakemake):
    location_data = pd.read_csv(snakemake.input["sensor_data"])
    accuracy_limit = snakemake.params["provider"]["ACCURACY_LIMIT"]

    datetime_start_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2} 00:00:00"
    datetime_end_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2} 23:59:59"
    segment_regex = ".*#{},{}".format(datetime_start_regex, datetime_end_regex)
    location_data = location_data[location_data["assigned_segments"].str.match(segment_regex)]
    loc_daily_data_len = len(location_data)

    location_data.query("accuracy < @accuracy_limit", inplace=True)

    features_to_compute = ["local_date", "hometime", "disttravelled", "rog", "maxdiam", "maxhomedist", "siglocsvisited", "avgflightlen", "stdflightlen", "avgflightdur", "stdflightdur", "probpause", "siglocentropy", "minsmissing", "circdnrtn", "wkenddayrtn", "minutes_data_used"]

    nrows = len(location_data)
    if loc_daily_data_len == 0:
        warnings.warn("Barnett's location features cannot be computed for data or time segments that do not span one or more entire days (00:00:00 to 23:59:59). This participant does not have location data or it spans less than 24 hours.")
        location_features = pd.DataFrame(columns=features_to_compute)
    elif nrows == 0:
        warnings.warn("Barnett's location features cannot be computed because there are no rows with an accuracy value lower than ACCURACY_LIMIT: {}".format(accuracy_limit))
        location_features = pd.DataFrame(columns=features_to_compute)
    else:
        location_minutes_used = location_data.groupby(["local_date", "local_hour"])[["local_minute"]].nunique().reset_index().groupby("local_date").sum()[["local_minute"]].rename(columns={"local_minute": "minutes_data_used"})
        timezone = mode(location_data["local_timezone"].values)
        location_df = location_data[["timestamp", "double_latitude", "double_longitude", "double_altitude", "accuracy"]]
        location_df.rename(columns={"double_latitude": "latitude", "double_longitude": "longitude", "double_altitude": "altitude"})
        output_mobility = run_barnett_features_for_rapids(location_df, accuracy_limit=accuracy_limit, timezone=timezone) #make local_date as the index for the output_mobility dataframe
        location_features = output_mobility.merge(location_minutes_used, on="local_date", how="left")

    location_features.reset_index(inplace=True)
    location_features.to_csv(snakemake.output[0], index=False)

barnett_daily_features(snakemake)