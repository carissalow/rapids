import pandas as pd
import numpy as np

pid = snakemake.params["pid"]
summarised = snakemake.params["summarised"]
targets_ratio_threshold = snakemake.params["targets_ratio_threshold"]
targets_value_threshold = snakemake.params["targets_value_threshold"]

if summarised == "summarised":
    targets = pd.DataFrame(columns=["pid", "target"])
    participant_info = pd.read_csv(snakemake.input["participant_info"])

    if not participant_info.empty:
        cesds = participant_info.loc[0, ["preop_cesd_total", "inpatient_cesd_total", "postop_cesd_total", "3month_cesd_total"]]
        # targets: 1 => 50% (ceiling) or more of available CESD scores were 16 or higher; 0 => otherwise
        num_threshold = int((cesds.count() + 1) * targets_ratio_threshold)
        target = 1 if cesds.apply(lambda x : 1 if x >= targets_value_threshold else 0).sum() >= num_threshold else 0
        targets.loc[0, :] = [pid, target]

targets.to_csv(snakemake.output[0], index=False)
