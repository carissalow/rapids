import pandas as pd
import numpy as np

pid = snakemake.params["pid"]
summarised = snakemake.params["summarised"]
participant_info = pd.read_csv(snakemake.input["participant_info"])

if summarised == "summarised":
    raise ValueError("Do not support summarised features for example dataset.")

elif summarised == "notsummarised":
    targets = participant_info[["local_date", "target"]]
    targets.insert(0, "pid", pid)

else:
    raise ValueError("SUMMARISED parameter in config.yaml can only be 'summarised' or 'notsummarised'")

targets.to_csv(snakemake.output[0], index=False)
