rule days_to_analyse:
    input:
        participant_info = "data/external/participant_info.csv",
        pid_file = "data/external/{pid}"
    params:
        days_before_surgery = config["METRICS_FOR_ANALYSIS"]["DAYS_BEFORE_SURGERY"],
        days_after_discharge = config["METRICS_FOR_ANALYSIS"]["DAYS_AFTER_DISCHARGE"],
        days_in_hospital= config["METRICS_FOR_ANALYSIS"]["DAYS_IN_HOSPITAL"]
    output:
        "data/interim/{pid}/days_to_analyse_{days_before_surgery}_{days_in_hospital}_{days_after_discharge}.csv"
    script:
        "../src/models/select_days_to_analyse.py"
