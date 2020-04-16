rule days_to_analyse:
    input:
        participant_info = "data/raw/{pid}/" + config["PARAMS_FOR_ANALYSIS"]["GROUNDTRUTH_TABLE"] + "_raw.csv"
    params:
        days_before_surgery = "{days_before_surgery}",
        days_in_hospital = "{days_in_hospital}",
        days_after_discharge= "{days_after_discharge}"
    output:
        "data/interim/{pid}/days_to_analyse_{days_before_surgery}_{days_in_hospital}_{days_after_discharge}.csv"
    script:
        "../src/models/select_days_to_analyse.py"

rule targets:
    input:
        participant_info = "data/raw/{pid}/" + config["PARAMS_FOR_ANALYSIS"]["GROUNDTRUTH_TABLE"] + "_raw.csv"
    params:
        pid = "{pid}",
        summarised = "{summarised}",
        targets_ratio_threshold = config["PARAMS_FOR_ANALYSIS"]["TARGETS_RATIO_THRESHOLD"],
        targets_value_threshold = config["PARAMS_FOR_ANALYSIS"]["TARGETS_VALUE_THRESHOLD"]
    output:
        "data/processed/{pid}/targets_{summarised}.csv"
    script:
        "../src/models/targets.py"

rule demographic_features:
    input:
        participant_info = "data/raw/{pid}/" + config["PARAMS_FOR_ANALYSIS"]["GROUNDTRUTH_TABLE"] + "_raw.csv"
    params:
        pid = "{pid}",
        features = config["PARAMS_FOR_ANALYSIS"]["DEMOGRAPHIC_FEATURES"]
    output:
        "data/processed/{pid}/demographic_features.csv"
    script:
        "../src/features/demographic_features.py"
