import pandas as pd

participant_info = pd.read_csv(snakemake.input["participant_info"])
summarised = snakemake.params["summarised"]
pid = snakemake.input["participant_info"].split("/")[2]

targets = pd.DataFrame({"pid": [pid], "target": [None]})
if summarised == "summarised":
    if not participant_info.empty:
        cesds = participant_info.loc[0, ["preop_cesd_total", "inpatient_cesd_total", "postop_cesd_total", "3month_cesd_total"]]
        # targets: 1 => 50% (ceiling) or more of available CESD scores were 16 or higher; 0 => otherwise
        threshold_num = (cesds.count() + 1) // 2
        threshold_cesd = 16
        target = 1 if cesds.apply(lambda x : 1 if x >= threshold_cesd else 0).sum() >= threshold_num else 0
        targets.loc[0, "target"] = target
targets.to_csv(snakemake.output[0], index=False)
