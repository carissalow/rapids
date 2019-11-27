import pandas as pd
import numpy as np
import datetime
from datetime import datetime, timedelta, time
from features_utils import splitOvernightEpisodes, splitMultiSegmentEpisodes

def getEpisodeDurationFeatures(screen_deltas, episode, metrics):
    screen_deltas_episode = screen_deltas[screen_deltas["episode"] == episode]
    duration_helper = pd.DataFrame()
    if "sumduration" in metrics:
        duration_helper = pd.concat([duration_helper, screen_deltas_episode.groupby(["local_start_date"]).sum()[["time_diff"]].rename(columns = {"time_diff": "screen_" + day_segment + "_sumduration" + episode})], axis = 1)
    if "maxduration" in metrics:
        duration_helper = pd.concat([duration_helper, screen_deltas_episode.groupby(["local_start_date"]).max()[["time_diff"]].rename(columns = {"time_diff": "screen_" + day_segment + "_maxduration" + episode})], axis = 1)        
    if "minduration" in metrics:
        duration_helper = pd.concat([duration_helper, screen_deltas_episode.groupby(["local_start_date"]).min()[["time_diff"]].rename(columns = {"time_diff": "screen_" + day_segment + "_minduration" + episode})], axis = 1)
    if "avgduration" in metrics:
        duration_helper = pd.concat([duration_helper, screen_deltas_episode.groupby(["local_start_date"]).mean()[["time_diff"]].rename(columns = {"time_diff":"screen_" + day_segment + "_avgduration" + episode})], axis = 1)
    if "stdduration" in metrics:
        duration_helper = pd.concat([duration_helper, screen_deltas_episode.groupby(["local_start_date"]).std()[["time_diff"]].rename(columns = {"time_diff":"screen_" + day_segment + "_stdduration" + episode})], axis = 1)
    
    duration_helper = duration_helper.fillna(0)
    return duration_helper

def getEventFeatures(screen_data, metrics_event):
    # get count_helper
    screen_status = screen_data.groupby(["local_date", "screen_status"]).count()[["timestamp"]].reset_index()
    count_on = screen_status[screen_status["screen_status"] == 0].set_index("local_date")[["timestamp"]].rename(columns = {"timestamp": "count_on"})
    count_off = screen_status[screen_status["screen_status"] == 1].set_index("local_date")[["timestamp"]].rename(columns = {"timestamp": "count_off"})
    count_lock = screen_status[screen_status["screen_status"] == 2].set_index("local_date")[["timestamp"]].rename(columns = {"timestamp": "count_lock"})
    count_unlock = screen_status[screen_status["screen_status"] == 3].set_index("local_date")[["timestamp"]].rename(columns = {"timestamp": "count_unlock"})
    
    count_helper = pd.concat([count_on, count_off, count_lock, count_unlock], axis = 1)
    count_helper = count_helper.fillna(0).astype(np.int64)

    # count on-off; unlock-lock
    count_helper["diff_count_on_off"] = count_helper["count_on"] - count_helper["count_off"]
    count_helper["diff_count_unlock_lock"] = count_helper["count_unlock"] - count_helper["count_lock"]
    
    event_features = pd.DataFrame()
    if "counton" in metrics_event:
        event_features["screen_" + day_segment + "_counton"] = count_helper[["count_on", "count_off"]].max(axis=1)
    if "countunlock" in metrics_event:
        event_features["screen_" + day_segment + "_countunlock"] = count_helper[["count_lock", "count_unlock"]].max(axis=1)

    ############################################################################################
    # check missing values
    event_features["screen_" + day_segment + "_diffcountonoff"] = count_helper["diff_count_on_off"]
    event_features["screen_" + day_segment + "_diffcountunlocklock"] = count_helper["diff_count_unlock_lock"]
    ############################################################################################

    return event_features

screen_data = pd.read_csv(snakemake.input["screen_events"], parse_dates=["local_date_time", "local_date"])
screen_deltas = pd.read_csv(snakemake.input["screen_deltas"], parse_dates=["local_start_date_time", "local_end_date_time", "local_start_date", "local_end_date"])
day_segment = snakemake.params["day_segment"]
metrics_event = snakemake.params["metrics_event"]
metrics_episode = snakemake.params["metrics_episode"]
episodes = snakemake.params["episodes"]

if screen_data.empty:
    metrics_episode_name = ["".join(metric) for metric in itertools.product(metrics_episode,episodes)]
    screen_features = pd.DataFrame(columns=["local_date"]+["screen_" + day_segment + "_" + x for x in metrics_event + metrics_episode_name])
else:    
    # drop consecutive duplicates of screen_status keeping the last one
    screen_data = screen_data.loc[(screen_data[["screen_status"]].shift(-1) != screen_data[["screen_status"]]).any(axis=1)].reset_index(drop=True)

    # preprocess day_segment and episodes
    screen_deltas = splitOvernightEpisodes(screen_deltas, [], ["episode"])
    if day_segment != "daily":
        screen_data = screen_data[screen_data["local_day_segment"] == day_segment]
        screen_deltas = splitMultiSegmentEpisodes(screen_deltas, day_segment, [])
    screen_deltas.set_index(["local_start_date"],inplace=True)

    # extract features for events and episodes
    event_features = getEventFeatures(screen_data, metrics_event)
    duration_features = pd.DataFrame()
    for episode in episodes:
        duration_features = pd.concat([duration_features, getEpisodeDurationFeatures(screen_deltas, episode, metrics_episode)], axis=1)

    screen_features = pd.concat([event_features, duration_features], axis = 1).fillna(0)
    screen_features.reset_index(inplace=True)

screen_features.to_csv(snakemake.output[0], index=False)