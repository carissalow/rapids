import pandas as pd
import itertools

def getEpisodeDurationFeatures(screen_data, time_segment, episode, features, reference_hour_first_use):
    screen_data_episode = screen_data[screen_data["episode"] == episode]
    duration_helper = pd.DataFrame()
    if "countepisode" in features:
        duration_helper = pd.concat([duration_helper, screen_data_episode.groupby(["local_segment"])[["duration"]].count().rename(columns = {"duration": "countepisode" + episode})], axis = 1)
    if "sumduration" in features:
        duration_helper = pd.concat([duration_helper, screen_data_episode.groupby(["local_segment"])[["duration"]].sum().rename(columns = {"duration": "sumduration" + episode})], axis = 1)
    if "maxduration" in features:
        duration_helper = pd.concat([duration_helper, screen_data_episode.groupby(["local_segment"])[["duration"]].max().rename(columns = {"duration": "maxduration" + episode})], axis = 1)        
    if "minduration" in features:
        duration_helper = pd.concat([duration_helper, screen_data_episode.groupby(["local_segment"])[["duration"]].min().rename(columns = {"duration": "minduration" + episode})], axis = 1)
    if "avgduration" in features:
        duration_helper = pd.concat([duration_helper, screen_data_episode.groupby(["local_segment"])[["duration"]].mean().rename(columns = {"duration":"avgduration" + episode})], axis = 1)
    if "stdduration" in features:
        duration_helper = pd.concat([duration_helper, screen_data_episode.groupby(["local_segment"])[["duration"]].std().rename(columns = {"duration":"stdduration" + episode})], axis = 1)
    if "firstuseafter" + "{0:0=2d}".format(reference_hour_first_use) in features:
        screen_data_episode_after_hour = screen_data_episode.copy()
        screen_data_episode_after_hour["hour"] = pd.to_datetime(screen_data_episode["local_start_date_time"]).dt.hour
        screen_data_episode_after_hour = screen_data_episode_after_hour[screen_data_episode_after_hour["hour"] >= reference_hour_first_use]

        duration_helper = pd.concat([duration_helper, pd.DataFrame(screen_data_episode_after_hour.groupby(["local_segment"])[["local_start_date_time"]].min().local_start_date_time.apply(lambda x: (x.to_pydatetime().hour - reference_hour_first_use) * 60 + x.to_pydatetime().minute + (x.to_pydatetime().second / 60))).rename(columns = {"local_start_date_time":"firstuseafter" + "{0:0=2d}".format(reference_hour_first_use) + episode})], axis = 1)
    return duration_helper


def rapids_features(sensor_data_files, time_segment, provider, filter_data_by_segment, *args, **kwargs):

    screen_data = pd.read_csv(sensor_data_files["sensor_episodes"])

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

    screen_features = pd.DataFrame(columns=["local_segment"] + features_to_compute)
    if not screen_data.empty:

        screen_data = filter_data_by_segment(screen_data, time_segment)
        
        if not screen_data.empty:
            if ignore_episodes_shorter_than > 0:
                screen_data = screen_data.query('@ignore_episodes_shorter_than <= duration')
            if ignore_episodes_longer_than > 0:
                screen_data = screen_data.query('duration <= @ignore_episodes_longer_than')

        if not screen_data.empty:
            screen_features = pd.DataFrame()
            for episode in episode_type_to_compute:
                screen_features = pd.concat([screen_features, getEpisodeDurationFeatures(screen_data, time_segment, episode, features_episodes_to_compute, reference_hour_first_use)], axis=1)

        if not screen_features.empty:
            screen_features.fillna(value={feature_name: 0 for feature_name in screen_features.columns if not feature_name.startswith(("stdduration", "firstuseafter"))}, inplace=True)
            screen_features = screen_features.reset_index()

    return screen_features

