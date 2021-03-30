import pandas as pd
import numpy as np
from importlib import import_module, util
from pathlib import Path


# import filter_data_by_segment from src/features/utils/utils.py
spec = util.spec_from_file_location("util", str(Path(snakemake.scriptdir).parent.parent / "features" / "utils" / "utils.py"))
mod = util.module_from_spec(spec)
spec.loader.exec_module(mod)
filter_data_by_segment = getattr(mod,  "filter_data_by_segment")

targets = pd.read_csv(snakemake.input["targets"])
time_segments_labels = pd.read_csv(snakemake.input["time_segments_labels"], header=0)

all_targets = pd.DataFrame(columns=["local_segment"])
for time_segment in time_segments_labels["label"]:
    filtered_targets = filter_data_by_segment(targets, time_segment)
    all_targets = all_targets.merge(filtered_targets, how="outer")

segment_colums = pd.DataFrame()
all_targets["local_segment"] = all_targets["local_segment"].str.replace(r'_RR\d+SS', '')
split_segemnt_columns = all_targets["local_segment"].str.split(pat="(.*)#(.*),(.*)", expand=True)
new_segment_columns = split_segemnt_columns.iloc[:,1:4] if split_segemnt_columns.shape[1] == 5 else pd.DataFrame(columns=["local_segment_label", "local_segment_start_datetime","local_segment_end_datetime"])
segment_colums[["local_segment_label", "local_segment_start_datetime", "local_segment_end_datetime"]] = new_segment_columns
for i in range(segment_colums.shape[1]):
    all_targets.insert(1 + i, segment_colums.columns[i], segment_colums[segment_colums.columns[i]])

all_targets.to_csv(snakemake.output[0], index=False)
