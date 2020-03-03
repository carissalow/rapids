import pandas as pd
import numpy as np
import datetime
import itertools
from datetime import datetime, timedelta, time
from features_utils import splitOvernightEpisodes, splitMultiSegmentEpisodes

def getEpisodeDurationFeatures(screen_deltas, episode, metrics, phone_sensed_bins, bin_size):
    screen_deltas_episode = screen_deltas[screen_deltas["episode"] == episode]
    duration_helper = pd.DataFrame()
    if "countepisode" in metrics:
        duration_helper = pd.concat([duration_helper, screen_deltas_episode.groupby(["local_start_date"]).count()[["time_diff"]].rename(columns = {"time_diff": "screen_" + day_segment + "_countepisode" + episode})], axis = 1)
    if "episodepersensedminutes" in metrics:
        for date, row in screen_deltas_episode.groupby(["local_start_date"]).count()[["time_diff"]].iterrows():
            sensed_minutes = phone_sensed_bins.loc[date, :].sum() * bin_size
            episode_per_sensedminutes = row["time_diff"] / (1 if sensed_minutes == 0 else sensed_minutes)
            duration_helper.loc[date, "screen_" + day_segment + "_episodepersensedminutes" + episode] = episode_per_sensedminutes
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


screen_deltas = pd.read_csv(snakemake.input["screen_deltas"], parse_dates=["local_start_date_time", "local_end_date_time", "local_start_date", "local_end_date"])
phone_sensed_bins = pd.read_csv(snakemake.input["phone_sensed_bins"], parse_dates=["local_date"], index_col="local_date")
phone_sensed_bins[phone_sensed_bins > 0] = 1

day_segment = snakemake.params["day_segment"]
metrics_deltas = snakemake.params["metrics_deltas"]
episode_types = snakemake.params["episode_types"]
bin_size = snakemake.params["bin_size"]

metrics_deltas_name = ["".join(metric) for metric in itertools.product(metrics_deltas, episode_types)]

screen_features = pd.DataFrame(columns=["local_date"]+["screen_" + day_segment + "_" + x for x in metrics_deltas_name])
if not screen_deltas.empty:
    # preprocess day_segment and episodes
    screen_deltas = splitOvernightEpisodes(screen_deltas, [], ["episode"])
    if (not screen_deltas.empty) and (day_segment != "daily"):
        screen_deltas = splitMultiSegmentEpisodes(screen_deltas, day_segment, [])
    screen_deltas.set_index(["local_start_date"],inplace=True)

    if not screen_deltas.empty:
        screen_features = pd.DataFrame()
        for episode in episode_types:
            screen_features = pd.concat([screen_features, getEpisodeDurationFeatures(screen_deltas, episode, metrics_deltas, phone_sensed_bins, bin_size)], axis=1)

    screen_features = screen_features.rename_axis("local_date").reset_index()

screen_features.to_csv(snakemake.output[0], index=False)