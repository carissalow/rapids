import pandas as pd


def resample_screen_deltas(screen_deltas):

    column_names = ("episode_id", "episode", "screen_sequence", "timestamp", "duration")
    records_resampled = []

    for _, row in screen_deltas.iterrows():
        episode_id, episode, screen_sequence = row["episode_id"], row["episode"], row["screen_sequence"]
        start_timestamp, end_timestamp = row["start_timestamp"], row["end_timestamp"]

        for timestamp in range(start_timestamp, end_timestamp, 1000 * 60):
            records_resampled.append((episode_id, episode, screen_sequence, timestamp, min(1, (end_timestamp - timestamp) / (1000 * 60))))

    return records_resampled, column_names


def resample_battery_deltas(battery_deltas):
    column_names = ("battery_diff", "timestamp")
    records_resampled = []

    for _, row in battery_deltas.iterrows():
        start_timestamp, end_timestamp = row["start_timestamp"], row["end_timestamp"]
        battery_diff = row["battery_diff"] / row["time_diff"]

        for timestamp in range(start_timestamp, end_timestamp, 1000 * 60):
            records_resampled.append((battery_diff, timestamp))

    return records_resampled, column_names


deltas = pd.read_csv(snakemake.input[0])
sensor = snakemake.params["sensor"]

if sensor == "battery":
    records_resampled, column_names = resample_battery_deltas(deltas)

if sensor == "screen":
    records_resampled, column_names = resample_screen_deltas(deltas)


deltas_resampled = pd.DataFrame(data=records_resampled, columns=column_names)
deltas_resampled.to_csv(snakemake.output[0], index=False)
