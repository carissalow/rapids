rule days_to_analyse:
    input:
        participant_info = "data/raw/{pid}/" + config["METRICS_FOR_ANALYSIS"]["GROUNDTRUTH_TABLE"] + "_raw.csv"
    params:
        days_before_surgery = "{days_before_surgery}",
        days_in_hospital = "{days_in_hospital}",
        days_after_discharge= "{days_after_discharge}"
    output:
        "data/interim/{pid}/days_to_analyse_{days_before_surgery}_{days_in_hospital}_{days_after_discharge}.csv"
    script:
        "../src/models/select_days_to_analyse.py"
