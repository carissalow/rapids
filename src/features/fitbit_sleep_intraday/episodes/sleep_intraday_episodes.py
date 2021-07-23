import pandas as pd
import numpy as np


def mergeSleepEpisodes(sleep_data, cols_for_groupby):
    sleep_episodes = pd.DataFrame(columns=["device_id", "type_episode_id", "level_episode_id", "level", "unified_level", "is_main_sleep", "type", "timestamp", "duration"])
    if not sleep_data.empty:
        sleep_data = sleep_data.groupby(by=cols_for_groupby, sort=False)
        sleep_episodes = sleep_data[["timestamp"]].first()
        sleep_episodes["duration"] = sleep_data["duration"].sum()
    
    return sleep_episodes


sleep_intraday = pd.read_csv(snakemake.input["sleep_intraday"])

# discard useless columns
for col in ["local_timezone", "local_date_time", "local_date", "local_time", "local_hour", "local_minute", "assigned_segments"]:
    del sleep_intraday[col]

# Extract "unified_level" based on "level" field
# For "classic" type, "unified_level" is one of {0, 1} where 0: awake {"awake" + "restless"}, 1: asleep {"asleep"}
# For "stages" type, "unified_level" is one of {0, 1} where 0: awake {"wake"}, 1: asleep {"deep" + "light" + "rem"}
sleep_intraday["unified_level"] = np.where(sleep_intraday["level"].isin(["awake", "restless", "wake"]), 0, 1)

# Put consecutive rows with the same "level" field together and merge episodes
sleep_intraday.insert(2, "level_episode_id", (sleep_intraday[["type_episode_id", "level"]] != sleep_intraday[["type_episode_id", "level"]].shift()).any(axis=1).cumsum())
sleep_intraday_episodes = mergeSleepEpisodes(sleep_intraday, ["device_id", "type_episode_id", "level_episode_id", "level", "unified_level", "is_main_sleep", "type"])


# Generate "start_timestamp" and "end_timestamp"
sleep_intraday_episodes["end_timestamp"] = sleep_intraday_episodes["timestamp"] + ((sleep_intraday_episodes["duration"] - 1) * 1000) + 999
sleep_intraday_episodes.rename(columns={"timestamp": "start_timestamp"}, inplace=True)

del sleep_intraday_episodes["duration"]

sleep_intraday_episodes.to_csv(snakemake.output[0], index=True)
