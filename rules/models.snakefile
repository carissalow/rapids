def input_merge_features_of_single_participant(wildcards):
    if wildcards.source == "phone_fitbit_features":
        return expand("data/processed/{pid}/{features}_{day_segment}.csv", pid=wildcards.pid, features=config["PARAMS_FOR_ANALYSIS"]["PHONE_FEATURES"] + config["PARAMS_FOR_ANALYSIS"]["FITBIT_FEATURES"], day_segment=wildcards.day_segment)
    else:
        return expand("data/processed/{pid}/{features}_{day_segment}.csv", pid=wildcards.pid, features=config["PARAMS_FOR_ANALYSIS"][wildcards.source.upper()], day_segment=wildcards.day_segment)

def optional_input_days_to_include(wildcards):
    if config["PARAMS_FOR_ANALYSIS"]["DAYS_TO_ANALYSE"]["ENABLED"]:
        # This input automatically trigers the rule days_to_analyse in mystudy.snakefile
        return ["data/interim/{pid}/days_to_analyse" + \
                    "_" + str(config["PARAMS_FOR_ANALYSIS"]["DAYS_TO_ANALYSE"]["DAYS_BEFORE_SURGERY"]) + \
                    "_" + str(config["PARAMS_FOR_ANALYSIS"]["DAYS_TO_ANALYSE"]["DAYS_IN_HOSPITAL"]) + \
                    "_" + str(config["PARAMS_FOR_ANALYSIS"]["DAYS_TO_ANALYSE"]["DAYS_AFTER_DISCHARGE"]) + ".csv"]
    else:
        return []

def optional_input_valid_sensed_days(wildcards):
    if config["PARAMS_FOR_ANALYSIS"]["DROP_VALID_SENSED_DAYS"]["ENABLED"]:
        # This input automatically trigers the rule phone_valid_sensed_days in preprocessing.snakefile
        return ["data/interim/{pid}/phone_valid_sensed_days.csv"]
    else:
        return []

rule merge_features_for_individual_model:
    input:
        feature_files = input_merge_features_of_single_participant,
        phone_valid_sensed_days = optional_input_valid_sensed_days,
        days_to_include = optional_input_days_to_include
    params:
        source = "{source}"
    output:
        "data/processed/{pid}/data_for_individual_model/{source}_{day_segment}_original.csv"
    script:
        "../src/models/merge_features_for_individual_model.R"

rule merge_features_for_population_model:
    input:
        feature_files = expand("data/processed/{pid}/data_for_individual_model/{{source}}_{{day_segment}}_original.csv", pid=config["PIDS"])
    output:
        "data/processed/data_for_population_model/{source}_{day_segment}_original.csv"
    script:
        "../src/models/merge_features_for_population_model.R"

rule merge_demographicfeatures_for_population_model:
    input:
        data_files = expand("data/processed/{pid}/demographic_features.csv", pid=config["PIDS"])
    output:
        "data/processed/data_for_population_model/demographic_features.csv"
    script:
        "../src/models/merge_data_for_population_model.py"

rule merge_targets_for_population_model:
    input:
        data_files = expand("data/processed/{pid}/targets_{{summarised}}.csv", pid=config["PIDS"])
    output:
        "data/processed/data_for_population_model/targets_{summarised}.csv"
    script:
        "../src/models/merge_data_for_population_model.py"

rule clean_features_for_individual_model:
    input:
        rules.merge_features_for_individual_model.output
    params:
        cols_nan_threshold = "{cols_nan_threshold}",
        cols_var_threshold = "{cols_var_threshold}",
        days_before_threshold = "{days_before_threshold}",
        days_after_threshold = "{days_after_threshold}",
        rows_nan_threshold = "{rows_nan_threshold}",
    output:
        "data/processed/{pid}/data_for_individual_model/{rows_nan_threshold}|{cols_nan_threshold}_{days_before_threshold}|{days_after_threshold}_{cols_var_threshold}/{source}_{day_segment}_clean.csv"
    script:
        "../src/models/clean_features_for_model.R"

rule clean_features_for_population_model:
    input:
        rules.merge_features_for_population_model.output
    params:
        cols_nan_threshold = "{cols_nan_threshold}",
        cols_var_threshold = "{cols_var_threshold}",
        days_before_threshold = "{days_before_threshold}",
        days_after_threshold = "{days_after_threshold}",
        rows_nan_threshold = "{rows_nan_threshold}",
    output:
        "data/processed/data_for_population_model/{rows_nan_threshold}|{cols_nan_threshold}_{days_before_threshold}|{days_after_threshold}_{cols_var_threshold}/{source}_{day_segment}_clean.csv"
    script:
        "../src/models/clean_features_for_model.R"

rule modeling:
    input:
        cleaned_features = "data/processed/data_for_population_model/{source}_{day_segment}_clean.csv",
        demographic_features = "data/processed/data_for_population_model/demographic_features.csv",
        targets = "data/processed/data_for_population_model/targets_{summarised}.csv",
    params:
        model = "{model}",
        cv_method = "{cv_method}",
        source = "{source}",
        day_segment = "{day_segment}",
        summarised = "{summarised}",
        scaler = "{scaler}",
        cols_var_threshold = config["PARAMS_FOR_ANALYSIS"]["COLS_VAR_THRESHOLD"],
        numerical_operators = config["PARAMS_FOR_ANALYSIS"]["NUMERICAL_OPERATORS"],
        categorical_operators = config["PARAMS_FOR_ANALYSIS"]["CATEGORICAL_OPERATORS"],
        categorical_demographic_features = config["PARAMS_FOR_ANALYSIS"]["CATEGORICAL_DEMOGRAPHIC_FEATURES"],
        model_hyperparams = config["PARAMS_FOR_ANALYSIS"]["MODEL_HYPERPARAMS"],
        rowsnan_colsnan_days_colsvar_threshold = "{rows_nan_threshold}_{cols_nan_threshold}_{days_before_threshold}|{days_after_threshold}_{cols_var_threshold}"
    output:
        fold_predictions = "data/processed/output_population_model/{rows_nan_threshold}_{cols_nan_threshold}_{days_before_threshold}|{days_after_threshold}_{cols_var_threshold}/{model}/{cv_method}/{source}_{day_segment}_{summarised}_{scaler}/fold_predictions.csv",
        fold_metrics = "data/processed/output_population_model/{rows_nan_threshold}_{cols_nan_threshold}_{days_before_threshold}|{days_after_threshold}_{cols_var_threshold}/{model}/{cv_method}/{source}_{day_segment}_{summarised}_{scaler}/fold_metrics.csv",
        overall_results = "data/processed/output_population_model/{rows_nan_threshold}_{cols_nan_threshold}_{days_before_threshold}|{days_after_threshold}_{cols_var_threshold}/{model}/{cv_method}/{source}_{day_segment}_{summarised}_{scaler}/overall_results.csv",
        fold_feature_importances = "data/processed/output_population_model/{rows_nan_threshold}_{cols_nan_threshold}_{days_before_threshold}|{days_after_threshold}_{cols_var_threshold}/{model}/{cv_method}/{source}_{day_segment}_{summarised}_{scaler}/fold_feature_importances.csv"
    script:
        "../src/models/modeling.py"