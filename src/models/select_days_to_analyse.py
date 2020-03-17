import pandas as pd
from datetime import timedelta

def appendDaysInRange(days_to_analyse, start_date, end_date):
    num_of_days = (end_date - start_date).days
    for day in range(num_of_days + 1):
            days_to_analyse = days_to_analyse.append({"days_to_analyse": start_date + timedelta(days = day)}, ignore_index=True)
    return days_to_analyse

days_before_surgery = snakemake.params["days_before_surgery"]
days_in_hospital = snakemake.params["days_in_hospital"]
days_after_discharge = snakemake.params["days_after_discharge"]

participant_info = pd.read_csv(snakemake.input["participant_info"], parse_dates=["surgery_date", "discharge_date"])
with open(snakemake.input["pid_file"], encoding="ISO-8859-1") as external_file:
    pid_file_content = external_file.readlines()
device_ids = pid_file_content[0].strip().split(",")

days_to_analyse = pd.DataFrame(columns = ["days_to_analyse"])
participant_info = participant_info[participant_info["device_id"].isin(device_ids)]
try:
    surgery_date, discharge_date = participant_info["surgery_date"].iloc[0].date(), participant_info["discharge_date"].iloc[0].date()
except:
    pass
else:
    start_date = surgery_date - timedelta(days = days_before_surgery)
    end_date = discharge_date + timedelta(days = days_after_discharge)

    days_to_analyse = appendDaysInRange(days_to_analyse, start_date, surgery_date - timedelta(days = 1))
    if days_in_hospital == "T":
        days_to_analyse = appendDaysInRange(days_to_analyse, surgery_date, discharge_date)
    days_to_analyse = appendDaysInRange(days_to_analyse, discharge_date + timedelta(days = 1), end_date)
    
days_to_analyse.to_csv(snakemake.output[0], index=False)
