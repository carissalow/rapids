import pandas as pd

calls = pd.read_csv(snakemake.input["calls"]).rename(columns={"timestamp": "start_timestamp"})
calls["end_timestamp"] = calls["start_timestamp"] + calls["call_duration"] * 1000
calls["episode_id"] = calls.index

calls[["episode_id", "device_id", "call_type", "trace", "start_timestamp", "end_timestamp"]].to_csv(snakemake.output[0], index=False)
