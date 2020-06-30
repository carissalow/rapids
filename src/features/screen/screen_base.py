import pandas as pd
import itertools
from features_utils import splitOvernightEpisodes, splitMultiSegmentEpisodes

EPOCH2HOUR = {"night": ("0_", "1_", "2_", "3_", "4_", "5_"),
              "morning": ("6_", "7_", "8_", "9_", "10_", "11_"),
              "afternoon": ("12_", "13_", "14_", "15_", "16_", "17_"),
              "evening": ("18_", "19_", "20_", "21_", "22_", "23_")}

def getEpisodeDurationFeatures(screen_data, day_segment, episode, features, phone_sensed_bins, bin_size, reference_hour_first_use):
    screen_data_episode = screen_data[screen_data["episode"] == episode]
    duration_helper = pd.DataFrame()
    if "countepisode" in features:
        duration_helper = pd.concat([duration_helper, screen_data_episode[["time_diff"]].groupby(["local_start_date"]).count().rename(columns = {"time_diff": "screen_" + day_segment + "_countepisode" + episode})], axis = 1)
    if "episodepersensedminutes" in features:
        for date, row in screen_data_episode[["time_diff"]].groupby(["local_start_date"]).count().iterrows():

            try:
                if day_segment == "daily":
                    sensed_minutes = phone_sensed_bins.loc[date, :].sum() * bin_size
                else:
                    sensed_minutes = phone_sensed_bins.loc[date, phone_sensed_bins.columns.str.startswith(EPOCH2HOUR[day_segment])].sum() * bin_size
            except:
                raise ValueError("You need to include the screen sensor in the list for phone_sensed_bins.")
                
            episode_per_sensedminutes = row["time_diff"] / (1 if sensed_minutes == 0 else sensed_minutes)
            duration_helper.loc[date, "screen_" + day_segment + "_episodepersensedminutes" + episode] = episode_per_sensedminutes
    if "sumduration" in features:
        duration_helper = pd.concat([duration_helper, screen_data_episode[["time_diff"]].groupby(["local_start_date"]).sum().rename(columns = {"time_diff": "screen_" + day_segment + "_sumduration" + episode})], axis = 1)
    if "maxduration" in features:
        duration_helper = pd.concat([duration_helper, screen_data_episode[["time_diff"]].groupby(["local_start_date"]).max().rename(columns = {"time_diff": "screen_" + day_segment + "_maxduration" + episode})], axis = 1)        
    if "minduration" in features:
        duration_helper = pd.concat([duration_helper, screen_data_episode[["time_diff"]].groupby(["local_start_date"]).min().rename(columns = {"time_diff": "screen_" + day_segment + "_minduration" + episode})], axis = 1)
    if "avgduration" in features:
        duration_helper = pd.concat([duration_helper, screen_data_episode[["time_diff"]].groupby(["local_start_date"]).mean().rename(columns = {"time_diff":"screen_" + day_segment + "_avgduration" + episode})], axis = 1)
    if "stdduration" in features:
        duration_helper = pd.concat([duration_helper, screen_data_episode[["time_diff"]].groupby(["local_start_date"]).std().rename(columns = {"time_diff":"screen_" + day_segment + "_stdduration" + episode})], axis = 1)
    if "firstuseafter" + "{0:0=2d}".format(reference_hour_first_use) in features:
        duration_helper = pd.concat([duration_helper, pd.DataFrame(screen_data_episode.groupby(["local_start_date"]).first()[["local_start_date_time"]].local_start_date_time.apply(lambda x: (x.to_pydatetime().hour - reference_hour_first_use) * 60 + x.to_pydatetime().minute + (x.to_pydatetime().second / 60))).rename(columns = {"local_start_date_time":"screen_" + day_segment + "_firstuseafter" + "{0:0=2d}".format(reference_hour_first_use) + episode})], axis = 1)
    return duration_helper


def base_screen_features(screen_data, phone_sensed_bins, day_segment, params):

    reference_hour_first_use = params["reference_hour_first_use"]
    bin_size = params["bin_size"]
    requested_features_deltas = params["requested_features_deltas"]
    requested_episode_types = params["requested_episode_types"]

    # name of the features this function can compute
    base_features_deltas = ["countepisode", "episodepersensedminutes", "sumduration", "maxduration", "minduration", "avgduration", "stdduration", "firstuseafter"]
    base_episode_type = ["unlock"]
    # the subset of requested features this function can compute
    features_deltas_to_compute = list(set(requested_features_deltas) & set(base_features_deltas))
    episode_type_to_compute = list(set(requested_episode_types) & set(base_episode_type))
    
    features_deltas_to_compute = ["firstuseafter" + "{0:0=2d}".format(reference_hour_first_use) if feature_name == "firstuseafter" else feature_name for feature_name in features_deltas_to_compute]
    features_to_compute = ["".join(feature) for feature in itertools.product(features_deltas_to_compute, episode_type_to_compute)]

    screen_features = pd.DataFrame(columns=["local_date"]+["screen_" + day_segment + "_" + x for x in features_to_compute])
    if not screen_data.empty:
        # preprocess day_segment and episodes
        screen_data = splitOvernightEpisodes(screen_data, [], ["episode"])
        if (not screen_data.empty) and (day_segment != "daily"):
            screen_data = splitMultiSegmentEpisodes(screen_data, day_segment, [])
        screen_data.set_index(["local_start_date"],inplace=True)

        if not screen_data.empty:
            screen_features = pd.DataFrame()
            for episode in episode_type_to_compute:
                screen_features = pd.concat([screen_features, getEpisodeDurationFeatures(screen_data, day_segment, episode, features_deltas_to_compute, phone_sensed_bins, bin_size, reference_hour_first_use)], axis=1)

        if not screen_features.empty:
            screen_features = screen_features.rename_axis("local_date").reset_index()

    return screen_features
