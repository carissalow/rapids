rule download_demographic_data:
    input:
        participant_file = "data/external/participant_files/{pid}.yaml",
        data = config["PARAMS_FOR_ANALYSIS"]["DEMOGRAPHIC"]["FOLDER"] + "/" + config["PARAMS_FOR_ANALYSIS"]["DEMOGRAPHIC"]["CONTAINER"]
    output:
        "data/raw/{pid}/participant_info_raw.csv"
    script:
        "../src/data/workflow_example/download_demographic_data.R"

rule demographic_features:
    input:
        participant_info = "data/raw/{pid}/participant_info_raw.csv"
    params:
        pid = "{pid}",
        features = config["PARAMS_FOR_ANALYSIS"]["DEMOGRAPHIC"]["FEATURES"]
    output:
        "data/processed/features/{pid}/demographic_features.csv"
    script:
        "../src/features/workflow_example/demographic_features.py"

rule download_target_data:
    input:
        participant_file = "data/external/participant_files/{pid}.yaml",
        data = config["PARAMS_FOR_ANALYSIS"]["TARGET"]["FOLDER"] + "/" + config["PARAMS_FOR_ANALYSIS"]["TARGET"]["CONTAINER"]
    output:
        "data/raw/{pid}/participant_target_raw.csv"
    script:
        "../src/data/workflow_example/download_target_data.R"

rule target_readable_datetime:
    input:
        sensor_input = "data/raw/{pid}/participant_target_raw.csv",
        time_segments = "data/interim/time_segments/{pid}_time_segments.csv",
        pid_file = "data/external/participant_files/{pid}.yaml",
        tzcodes_file = input_tzcodes_file,
    params:
        device_type = "fitbit",
        timezone_parameters = config["TIMEZONE"],
        pid = "{pid}",
        time_segments_type = config["TIME_SEGMENTS"]["TYPE"],
        include_past_periodic_segments = config["TIME_SEGMENTS"]["INCLUDE_PAST_PERIODIC_SEGMENTS"]
    output:
        target_with_datetime = "data/raw/{pid}/participant_target_with_datetime.csv",
        flag_file = touch("data/raw/{pid}/participant_target_with_datetime.done")
    script:
        "../src/data/datetime/readable_datetime.R"

rule parse_targets:
    input:
        targets = "data/raw/{pid}/participant_target_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    output:
        "data/processed/targets/{pid}/parsed_targets.csv"
    script:
        "../src/models/workflow_example/parse_targets.py"

rule merge_features_and_targets_for_individual_model:
    input:
        cleaned_sensor_features = "data/processed/features/{pid}/all_sensor_features_cleaned_rapids.csv",
        targets = "data/processed/targets/{pid}/parsed_targets.csv",
    output:
        "data/processed/models/individual_model/{pid}/input.csv"
    script:
        "../src/models/workflow_example/merge_features_and_targets_for_individual_model.py"

rule merge_features_and_targets_for_population_model:
    input:
        cleaned_sensor_features = "data/processed/features/all_participants/all_sensor_features_cleaned_rapids.csv",
        demographic_features = expand("data/processed/features/{pid}/demographic_features.csv", pid=config["PIDS"]),
        targets = expand("data/processed/targets/{pid}/parsed_targets.csv", pid=config["PIDS"]),
    output:
        "data/processed/models/population_model/input.csv"
    script:
        "../src/models/workflow_example/merge_features_and_targets_for_population_model.py"

rule baselines_for_individual_model:
    input:
        "data/processed/models/individual_model/{pid}/input.csv"
    params:
        cv_method = "{cv_method}",
        colnames_demographic_features = config["PARAMS_FOR_ANALYSIS"]["DEMOGRAPHIC"]["FEATURES"],
    output:
        "data/processed/models/individual_model/{pid}/output_{cv_method}/baselines.csv"
    log:
        "data/processed/models/individual_model/{pid}/output_{cv_method}/baselines_notes.log"
    script:
        "../src/models/workflow_example/baselines.py"

rule baselines_for_population_model:
    input:
        "data/processed/models/population_model/input.csv"
    params:
        cv_method = "{cv_method}",
        colnames_demographic_features = config["PARAMS_FOR_ANALYSIS"]["DEMOGRAPHIC"]["FEATURES"],
    output:
        "data/processed/models/population_model/output_{cv_method}/baselines.csv"
    log:
        "data/processed/models/population_model/output_{cv_method}/baselines_notes.log"
    script:
        "../src/models/workflow_example/baselines.py"

rule modelling_for_individual_participants:
    input:
        data = "data/processed/models/individual_model/{pid}/input.csv"
    params:
        model = "{model}",
        cv_method = "{cv_method}",
        scaler = "{scaler}",
        categorical_operators = config["PARAMS_FOR_ANALYSIS"]["CATEGORICAL_OPERATORS"],
        categorical_demographic_features = config["PARAMS_FOR_ANALYSIS"]["DEMOGRAPHIC"]["CATEGORICAL_FEATURES"],
        model_hyperparams = config["PARAMS_FOR_ANALYSIS"]["MODEL_HYPERPARAMS"],
    output:
        fold_predictions = "data/processed/models/individual_model/{pid}/output_{cv_method}/{model}/{scaler}/fold_predictions.csv",
        fold_metrics = "data/processed/models/individual_model/{pid}/output_{cv_method}/{model}/{scaler}/fold_metrics.csv",
        overall_results = "data/processed/models/individual_model/{pid}/output_{cv_method}/{model}/{scaler}/overall_results.csv",
        fold_feature_importances = "data/processed/models/individual_model/{pid}/output_{cv_method}/{model}/{scaler}/fold_feature_importances.csv"
    log:
        "data/processed/models/individual_model/{pid}/output_{cv_method}/{model}/{scaler}/notes.log"
    script:
        "../src/models/workflow_example/modelling.py"

rule modelling_for_all_participants:
    input:
        data = "data/processed/models/population_model/input.csv"
    params:
        model = "{model}",
        cv_method = "{cv_method}",
        scaler = "{scaler}",
        categorical_operators = config["PARAMS_FOR_ANALYSIS"]["CATEGORICAL_OPERATORS"],
        categorical_demographic_features = config["PARAMS_FOR_ANALYSIS"]["DEMOGRAPHIC"]["CATEGORICAL_FEATURES"],
        model_hyperparams = config["PARAMS_FOR_ANALYSIS"]["MODEL_HYPERPARAMS"],
    output:
        fold_predictions = "data/processed/models/population_model/output_{cv_method}/{model}/{scaler}/fold_predictions.csv",
        fold_metrics = "data/processed/models/population_model/output_{cv_method}/{model}/{scaler}/fold_metrics.csv",
        overall_results = "data/processed/models/population_model/output_{cv_method}/{model}/{scaler}/overall_results.csv",
        fold_feature_importances = "data/processed/models/population_model/output_{cv_method}/{model}/{scaler}/fold_feature_importances.csv"
    log:
        "data/processed/models/population_model/output_{cv_method}/{model}/{scaler}/notes.log"
    script:
        "../src/models/workflow_example/modelling.py"
