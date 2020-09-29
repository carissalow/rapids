import pandas as pd
import itertools

def getEpisodeDurationFeatures(screen_data, day_segment, episode, features, reference_hour_first_use):
    screen_data_episode = screen_data[screen_data["episode"] == episode]
    duration_helper = pd.DataFrame()
    if "countepisode" in features:
        duration_helper = pd.concat([duration_helper, screen_data_episode.groupby(["local_segment"])[["time_diff"]].count().rename(columns = {"time_diff": "screen_rapids_countepisode" + episode})], axis = 1)
    if "sumduration" in features:
        duration_helper = pd.concat([duration_helper, screen_data_episode.groupby(["local_segment"])[["time_diff"]].sum().rename(columns = {"time_diff": "screen_rapids_sumduration" + episode})], axis = 1)
    if "maxduration" in features:
        duration_helper = pd.concat([duration_helper, screen_data_episode.groupby(["local_segment"])[["time_diff"]].max().rename(columns = {"time_diff": "screen_rapids_maxduration" + episode})], axis = 1)        
    if "minduration" in features:
        duration_helper = pd.concat([duration_helper, screen_data_episode.groupby(["local_segment"])[["time_diff"]].min().rename(columns = {"time_diff": "screen_rapids_minduration" + episode})], axis = 1)
    if "avgduration" in features:
        duration_helper = pd.concat([duration_helper, screen_data_episode.groupby(["local_segment"])[["time_diff"]].mean().rename(columns = {"time_diff":"screen_rapids_avgduration" + episode})], axis = 1)
    if "stdduration" in features:
        duration_helper = pd.concat([duration_helper, screen_data_episode.groupby(["local_segment"])[["time_diff"]].std().rename(columns = {"time_diff":"screen_rapids_stdduration" + episode})], axis = 1)
    if "firstuseafter" + "{0:0=2d}".format(reference_hour_first_use) in features:
        screen_data_episode_after_hour = screen_data_episode.copy()
        screen_data_episode_after_hour["hour"] = pd.to_datetime(screen_data_episode["local_start_date_time"]).dt.hour
        screen_data_episode_after_hour = screen_data_episode_after_hour[screen_data_episode_after_hour["hour"] >= reference_hour_first_use]

        duration_helper = pd.concat([duration_helper, pd.DataFrame(screen_data_episode_after_hour.groupby(["local_segment"])[["local_start_date_time"]].min().local_start_date_time.apply(lambda x: (x.to_pydatetime().hour - reference_hour_first_use) * 60 + x.to_pydatetime().minute + (x.to_pydatetime().second / 60))).rename(columns = {"local_start_date_time":"screen_rapids_firstuseafter" + "{0:0=2d}".format(reference_hour_first_use) + episode})], axis = 1)
    return duration_helper


def rapids_features(screen_data, day_segment, provider, filter_data_by_segment, *args, **kwargs):

    reference_hour_first_use = provider["REFERENCE_HOUR_FIRST_USE"]
    requested_features_episodes = provider["FEATURES"]
    requested_episode_types = provider["EPISODE_TYPES"]
    ignore_episodes_shorter_than = provider["IGNORE_EPISODES_SHORTER_THAN"]
    ignore_episodes_longer_than = provider["IGNORE_EPISODES_LONGER_THAN"]

    # name of the features this function can compute
    base_features_episodes = ["countepisode", "episodepersensedminutes", "sumduration", "maxduration", "minduration", "avgduration", "stdduration", "firstuseafter"]
    base_episode_type = ["unlock"]
    # the subset of requested features this function can compute
    features_episodes_to_compute = list(set(requested_features_episodes) & set(base_features_episodes))
    episode_type_to_compute = list(set(requested_episode_types) & set(base_episode_type))
    
    features_episodes_to_compute = ["firstuseafter" + "{0:0=2d}".format(reference_hour_first_use) if feature_name == "firstuseafter" else feature_name for feature_name in features_episodes_to_compute]
    features_to_compute = ["".join(feature) for feature in itertools.product(features_episodes_to_compute, episode_type_to_compute)]

    screen_features = pd.DataFrame(columns=["local_segment"]+["screen_rapids_" + x for x in features_to_compute])
    if not screen_data.empty:

        screen_data = filter_data_by_segment(screen_data, day_segment)
        if not screen_data.empty:
            # chunk_episodes
            screen_data = kwargs["chunk_episodes"](screen_data)

            if ignore_episodes_shorter_than > 0:
                screen_data = screen_data.query('@ignore_episodes_shorter_than <= time_diff')
            if ignore_episodes_longer_than > 0:
                screen_data = screen_data.query('time_diff <= @ignore_episodes_longer_than')

        if not screen_data.empty:
            screen_features = pd.DataFrame()
            for episode in episode_type_to_compute:
                screen_features = pd.concat([screen_features, getEpisodeDurationFeatures(screen_data, day_segment, episode, features_episodes_to_compute, reference_hour_first_use)], axis=1)

        if not screen_features.empty:
            screen_features = screen_features.reset_index()

    return screen_features
