import pandas as pd

pid = snakemake.params["pid"]
requested_features = snakemake.params["features"]
demographic_features = pd.DataFrame(columns=requested_features)

participant_info = pd.read_csv(snakemake.input["participant_info"], parse_dates=["surgery_date", "discharge_date"])
if not participant_info.empty:
    if "age" in requested_features:
        demographic_features.loc[0, "age"] = participant_info.loc[0, "age"]
    if "gender" in requested_features:
        demographic_features.loc[0, "gender"] = participant_info.loc[0, "gender"]
    if "inpatientdays" in requested_features:
        demographic_features.loc[0, "inpatientdays"] = (participant_info.loc[0, "discharge_date"] - participant_info.loc[0, "surgery_date"]).days

demographic_features.to_csv(snakemake.output[0], index=False)
