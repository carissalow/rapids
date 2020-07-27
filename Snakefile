configfile: "config.yaml"
include: "rules/renv.snakefile"
include: "rules/preprocessing.snakefile"
include: "rules/features.snakefile"
include: "rules/models.snakefile"
include: "rules/reports.snakefile"
include: "rules/mystudy.snakefile" # You can add snakfiles with rules tailored to your project

import itertools

files_to_compute = []

if len(config["PIDS"]) == 0:
    raise ValueError("Add participants IDs to PIDS in config.yaml. Remember to create their participant files in data/external")

if config["PHONE_VALID_SENSED_BINS"]["COMPUTE"]:
    if len(config["PHONE_VALID_SENSED_BINS"]["TABLES"]) == 0:
            raise ValueError("If you want to compute PHONE_VALID_SENSED_BINS, you need to add at least one table to [PHONE_VALID_SENSED_BINS][TABLES] in config.yaml")
    files_to_compute.extend(expand("data/interim/{pid}/phone_sensed_bins.csv", pid=config["PIDS"]))

if config["PHONE_VALID_SENSED_DAYS"]["COMPUTE"]:
    if len(config["PHONE_VALID_SENSED_BINS"]["TABLES"]) == 0:
            raise ValueError("If you want to compute PHONE_VALID_SENSED_DAYS, you need to add at least one table to [PHONE_VALID_SENSED_BINS][TABLES] in config.yaml")
    files_to_compute.extend(expand("data/interim/{pid}/phone_sensed_bins.csv", pid=config["PIDS"]))
    files_to_compute.extend(expand("data/interim/{pid}/phone_valid_sensed_days_{min_valid_hours_per_day}hours_{min_valid_bins_per_hour}bins.csv",
                                pid=config["PIDS"],
                                min_valid_hours_per_day=config["PHONE_VALID_SENSED_DAYS"]["MIN_VALID_HOURS_PER_DAY"],
                                min_valid_bins_per_hour=config["PHONE_VALID_SENSED_DAYS"]["MIN_VALID_BINS_PER_HOUR"]))

if config["MESSAGES"]["COMPUTE"]:
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["MESSAGES"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["MESSAGES"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/processed/{pid}/messages_{messages_type}_{day_segment}.csv", pid=config["PIDS"], messages_type = config["MESSAGES"]["TYPES"], day_segment = config["MESSAGES"]["DAY_SEGMENTS"]))

if config["CALLS"]["COMPUTE"]:
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["CALLS"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["CALLS"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime_unified.csv", pid=config["PIDS"], sensor=config["CALLS"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/processed/{pid}/calls_{call_type}_{day_segment}.csv", pid=config["PIDS"], call_type=config["CALLS"]["TYPES"], day_segment = config["CALLS"]["DAY_SEGMENTS"]))

if config["BARNETT_LOCATION"]["COMPUTE"]:
    # TODO add files_to_compute.extend(optional_location_input(None))
    if config["BARNETT_LOCATION"]["LOCATIONS_TO_USE"] == "RESAMPLE_FUSED":
        if config["BARNETT_LOCATION"]["DB_TABLE"] in config["PHONE_VALID_SENSED_BINS"]["TABLES"]:
            files_to_compute.extend(expand("data/interim/{pid}/phone_sensed_bins.csv", pid=config["PIDS"]))
            files_to_compute.extend(expand("data/raw/{pid}/{sensor}_resampled.csv", pid=config["PIDS"], sensor=config["BARNETT_LOCATION"]["DB_TABLE"]))
        else:
            raise ValueError("Error: Add your locations table (and as many sensor tables as you have) to [PHONE_VALID_SENSED_BINS][TABLES] in config.yaml. This is necessary to compute phone_sensed_bins (bins of time when the smartphone was sensing data) which is used to resample fused location data (RESAMPLED_FUSED)")            
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["BARNETT_LOCATION"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["BARNETT_LOCATION"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/processed/{pid}/location_barnett_{day_segment}.csv", pid=config["PIDS"], day_segment = config["BARNETT_LOCATION"]["DAY_SEGMENTS"]))

if config["BLUETOOTH"]["COMPUTE"]:
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["BLUETOOTH"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["BLUETOOTH"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/processed/{pid}/bluetooth_{day_segment}.csv", pid=config["PIDS"], day_segment = config["BLUETOOTH"]["DAY_SEGMENTS"]))

if config["ACTIVITY_RECOGNITION"]["COMPUTE"]:
    # TODO add files_to_compute.extend(optional_ar_input(None)), the Android or iOS table gets processed depending on each participant
    files_to_compute.extend(expand("data/processed/{pid}/activity_recognition_{day_segment}.csv",pid=config["PIDS"], day_segment = config["ACTIVITY_RECOGNITION"]["DAY_SEGMENTS"]))

if config["BATTERY"]["COMPUTE"]:
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["BATTERY"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["BATTERY"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime_unified.csv", pid=config["PIDS"], sensor=config["BATTERY"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/processed/{pid}/battery_deltas.csv", pid=config["PIDS"]))
    files_to_compute.extend(expand("data/processed/{pid}/battery_{day_segment}.csv", pid = config["PIDS"], day_segment = config["BATTERY"]["DAY_SEGMENTS"]))

if config["SCREEN"]["COMPUTE"]:
    if config["SCREEN"]["DB_TABLE"] in config["PHONE_VALID_SENSED_BINS"]["TABLES"]:
        files_to_compute.extend(expand("data/interim/{pid}/phone_sensed_bins.csv", pid=config["PIDS"]))
    else:
        raise ValueError("Error: Add your screen table (and as many sensor tables as you have) to [PHONE_VALID_SENSED_BINS][TABLES] in config.yaml. This is necessary to compute phone_sensed_bins (bins of time when the smartphone was sensing data)")
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["SCREEN"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["SCREEN"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime_unified.csv", pid=config["PIDS"], sensor=config["SCREEN"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/processed/{pid}/screen_deltas.csv", pid=config["PIDS"]))
    files_to_compute.extend(expand("data/processed/{pid}/screen_{day_segment}.csv", pid = config["PIDS"], day_segment = config["SCREEN"]["DAY_SEGMENTS"]))

if config["LIGHT"]["COMPUTE"]:
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["LIGHT"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["LIGHT"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/processed/{pid}/light_{day_segment}.csv", pid = config["PIDS"], day_segment = config["LIGHT"]["DAY_SEGMENTS"]))

if config["ACCELEROMETER"]["COMPUTE"]:
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["ACCELEROMETER"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["ACCELEROMETER"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/processed/{pid}/accelerometer_{day_segment}.csv", pid = config["PIDS"], day_segment = config["ACCELEROMETER"]["DAY_SEGMENTS"]))

if config["APPLICATIONS_FOREGROUND"]["COMPUTE"]:
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["APPLICATIONS_FOREGROUND"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["APPLICATIONS_FOREGROUND"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/interim/{pid}/{sensor}_with_datetime_with_genre.csv", pid=config["PIDS"], sensor=config["APPLICATIONS_FOREGROUND"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/processed/{pid}/applications_foreground_{day_segment}.csv", pid = config["PIDS"], day_segment = config["APPLICATIONS_FOREGROUND"]["DAY_SEGMENTS"]))

if config["WIFI"]["COMPUTE"]:
    if len(config["WIFI"]["DB_TABLE"]["VISIBLE_ACCESS_POINTS"]) > 0:
        files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["WIFI"]["DB_TABLE"]["VISIBLE_ACCESS_POINTS"]))
        files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["WIFI"]["DB_TABLE"]["VISIBLE_ACCESS_POINTS"]))
        files_to_compute.extend(expand("data/processed/{pid}/wifi_{day_segment}.csv", pid = config["PIDS"], day_segment = config["WIFI"]["DAY_SEGMENTS"]))

    if len(config["WIFI"]["DB_TABLE"]["CONNECTED_ACCESS_POINTS"]) > 0:
        files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["WIFI"]["DB_TABLE"]["CONNECTED_ACCESS_POINTS"]))
        files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["WIFI"]["DB_TABLE"]["CONNECTED_ACCESS_POINTS"]))
        files_to_compute.extend(expand("data/processed/{pid}/wifi_{day_segment}.csv", pid = config["PIDS"], day_segment = config["WIFI"]["DAY_SEGMENTS"]))

if config["HEARTRATE"]["COMPUTE"]:
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["HEARTRATE"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/fitbit_heartrate_{fitbit_data_type}_with_datetime.csv", pid=config["PIDS"], fitbit_data_type=["summary", "intraday"]))
    files_to_compute.extend(expand("data/processed/{pid}/fitbit_heartrate_{day_segment}.csv", pid = config["PIDS"], day_segment = config["HEARTRATE"]["DAY_SEGMENTS"]))

if config["STEP"]["COMPUTE"]:
    if config["STEP"]["EXCLUDE_SLEEP"]["EXCLUDE"] == True and config["STEP"]["EXCLUDE_SLEEP"]["TYPE"] == "FITBIT_BASED":
        files_to_compute.extend(expand("data/raw/{pid}/fitbit_sleep_{fitbit_data_type}_with_datetime.csv", pid=config["PIDS"], fitbit_data_type=["summary"]))
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["STEP"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/fitbit_step_{fitbit_data_type}_with_datetime.csv", pid=config["PIDS"], fitbit_data_type=["intraday"]))
    files_to_compute.extend(expand("data/processed/{pid}/fitbit_step_{day_segment}.csv", pid = config["PIDS"], day_segment = config["STEP"]["DAY_SEGMENTS"]))

if config["SLEEP"]["COMPUTE"]:
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["SLEEP"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/fitbit_sleep_{fitbit_data_type}_with_datetime.csv", pid=config["PIDS"], fitbit_data_type=["intraday", "summary"]))
    files_to_compute.extend(expand("data/processed/{pid}/fitbit_sleep_{day_segment}.csv", pid = config["PIDS"], day_segment = config["SLEEP"]["DAY_SEGMENTS"]))

if config["CONVERSATION"]["COMPUTE"]:
    # TODO add files_to_compute.extend(optional_conversation_input(None)), the Android or iOS table gets processed depending on each participant
    files_to_compute.extend(expand("data/processed/{pid}/conversation_{day_segment}.csv",pid=config["PIDS"], day_segment = config["CONVERSATION"]["DAY_SEGMENTS"]))

if config["DORYAB_LOCATION"]["COMPUTE"]:
    if config["DORYAB_LOCATION"]["LOCATIONS_TO_USE"] == "RESAMPLE_FUSED":
        if config["DORYAB_LOCATION"]["DB_TABLE"] in config["PHONE_VALID_SENSED_BINS"]["TABLES"]:
            files_to_compute.extend(expand("data/interim/{pid}/phone_sensed_bins.csv", pid=config["PIDS"]))
            files_to_compute.extend(expand("data/raw/{pid}/{sensor}_resampled.csv", pid=config["PIDS"], sensor=config["DORYAB_LOCATION"]["DB_TABLE"]))
        else:
            raise ValueError("Error: Add your locations table (and as many sensor tables as you have) to [PHONE_VALID_SENSED_BINS][TABLES] in config.yaml. This is necessary to compute phone_sensed_bins (bins of time when the smartphone was sensing data) which is used to resample fused location data (RESAMPLED_FUSED)")      
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["DORYAB_LOCATION"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["DORYAB_LOCATION"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/processed/{pid}/location_doryab_{segment}.csv", pid=config["PIDS"], segment = config["DORYAB_LOCATION"]["DAY_SEGMENTS"]))

# visualization for data exploration
if config["HEATMAP_FEATURES_CORRELATIONS"]["PLOT"]:
    files_to_compute.extend(expand("reports/data_exploration/{min_valid_hours_per_day}hours_{min_valid_bins_per_hour}bins/heatmap_features_correlations.html", min_valid_hours_per_day=config["HEATMAP_FEATURES_CORRELATIONS"]["MIN_VALID_HOURS_PER_DAY"], min_valid_bins_per_hour=config["PHONE_VALID_SENSED_DAYS"]["MIN_VALID_BINS_PER_HOUR"]))
    
if config["HISTOGRAM_VALID_SENSED_HOURS"]["PLOT"]:
    files_to_compute.extend(expand("reports/data_exploration/{min_valid_hours_per_day}hours_{min_valid_bins_per_hour}bins/histogram_valid_sensed_hours.html", min_valid_hours_per_day=config["HISTOGRAM_VALID_SENSED_HOURS"]["MIN_VALID_HOURS_PER_DAY"], min_valid_bins_per_hour=config["PHONE_VALID_SENSED_DAYS"]["MIN_VALID_BINS_PER_HOUR"]))

if config["HEATMAP_DAYS_BY_SENSORS"]["PLOT"]:
    files_to_compute.extend(expand("reports/interim/{min_valid_hours_per_day}hours_{min_valid_bins_per_hour}bins/{pid}/heatmap_days_by_sensors.html", pid=config["PIDS"], min_valid_hours_per_day=config["HEATMAP_DAYS_BY_SENSORS"]["MIN_VALID_HOURS_PER_DAY"], min_valid_bins_per_hour=config["PHONE_VALID_SENSED_DAYS"]["MIN_VALID_BINS_PER_HOUR"]))
    files_to_compute.extend(expand("reports/data_exploration/{min_valid_hours_per_day}hours_{min_valid_bins_per_hour}bins/heatmap_days_by_sensors_all_participants.html", min_valid_hours_per_day=config["HEATMAP_DAYS_BY_SENSORS"]["MIN_VALID_HOURS_PER_DAY"], min_valid_bins_per_hour=config["PHONE_VALID_SENSED_DAYS"]["MIN_VALID_BINS_PER_HOUR"]))

if config["HEATMAP_SENSED_BINS"]["PLOT"]:
    files_to_compute.extend(expand("reports/interim/heatmap_sensed_bins/{pid}/heatmap_sensed_bins.html", pid=config["PIDS"]))
    files_to_compute.extend(["reports/data_exploration/heatmap_sensed_bins_all_participants.html"])

if config["OVERALL_COMPLIANCE_HEATMAP"]["PLOT"]:
    files_to_compute.extend(expand("reports/data_exploration/{min_valid_hours_per_day}hours_{min_valid_bins_per_hour}bins/overall_compliance_heatmap.html", min_valid_hours_per_day=config["OVERALL_COMPLIANCE_HEATMAP"]["MIN_VALID_HOURS_PER_DAY"], min_valid_bins_per_hour=config["PHONE_VALID_SENSED_DAYS"]["MIN_VALID_BINS_PER_HOUR"]))

# analysis example
if config["PARAMS_FOR_ANALYSIS"]["COMPUTE"]:
    rows_nan_threshold = config["PARAMS_FOR_ANALYSIS"]["ROWS_NAN_THRESHOLD"]
    cols_nan_threshold = config["PARAMS_FOR_ANALYSIS"]["COLS_NAN_THRESHOLD"]
    models, scalers, rows_nan_thresholds, cols_nan_thresholds = [], [], [], []
    for model_name in config["PARAMS_FOR_ANALYSIS"]["MODEL_NAMES"]:
        models = models + [model_name] * len(config["PARAMS_FOR_ANALYSIS"]["MODEL_SCALER"][model_name]) * len(rows_nan_threshold)
        scalers = scalers + config["PARAMS_FOR_ANALYSIS"]["MODEL_SCALER"][model_name] * len(rows_nan_threshold)
        rows_nan_thresholds = rows_nan_thresholds + list(itertools.chain.from_iterable([threshold] * len(config["PARAMS_FOR_ANALYSIS"]["MODEL_SCALER"][model_name]) for threshold in rows_nan_threshold))
        cols_nan_thresholds = cols_nan_thresholds + list(itertools.chain.from_iterable([threshold] * len(config["PARAMS_FOR_ANALYSIS"]["MODEL_SCALER"][model_name]) for threshold in cols_nan_threshold))
    results = config["PARAMS_FOR_ANALYSIS"]["RESULT_COMPONENTS"] + ["merged_population_model_results"]

    files_to_compute.extend(expand("data/processed/{pid}/data_for_individual_model/{min_valid_hours_per_day}hours_{min_valid_bins_per_hour}bins/{source}_{day_segment}_original.csv",
                                pid = config["PIDS"],
                                min_valid_hours_per_day=config["OVERALL_COMPLIANCE_HEATMAP"]["MIN_VALID_HOURS_PER_DAY"],
                                min_valid_bins_per_hour=config["PHONE_VALID_SENSED_DAYS"]["MIN_VALID_BINS_PER_HOUR"],
                                source = config["PARAMS_FOR_ANALYSIS"]["SOURCES"],
                                day_segment = config["PARAMS_FOR_ANALYSIS"]["DAY_SEGMENTS"]))
    files_to_compute.extend(expand("data/processed/data_for_population_model/{min_valid_hours_per_day}hours_{min_valid_bins_per_hour}bins/{source}_{day_segment}_original.csv",
                                min_valid_hours_per_day=config["OVERALL_COMPLIANCE_HEATMAP"]["MIN_VALID_HOURS_PER_DAY"],
                                min_valid_bins_per_hour=config["PHONE_VALID_SENSED_DAYS"]["MIN_VALID_BINS_PER_HOUR"],
                                source = config["PARAMS_FOR_ANALYSIS"]["SOURCES"],
                                day_segment = config["PARAMS_FOR_ANALYSIS"]["DAY_SEGMENTS"]))
    files_to_compute.extend(expand(
                                expand("data/processed/{pid}/data_for_individual_model/{min_valid_hours_per_day}hours_{min_valid_bins_per_hour}bins/{{rows_nan_threshold}}|{{cols_nan_threshold}}_{days_before_threshold}|{days_after_threshold}_{cols_var_threshold}/{source}_{day_segment}_clean.csv",
                                    pid = config["PIDS"],
                                    min_valid_hours_per_day=config["OVERALL_COMPLIANCE_HEATMAP"]["MIN_VALID_HOURS_PER_DAY"],
                                    min_valid_bins_per_hour=config["PHONE_VALID_SENSED_DAYS"]["MIN_VALID_BINS_PER_HOUR"],
                                    days_before_threshold = config["PARAMS_FOR_ANALYSIS"]["PARTICIPANT_DAYS_BEFORE_THRESHOLD"],
                                    days_after_threshold = config["PARAMS_FOR_ANALYSIS"]["PARTICIPANT_DAYS_AFTER_THRESHOLD"],
                                    cols_var_threshold = config["PARAMS_FOR_ANALYSIS"]["COLS_VAR_THRESHOLD"],
                                    source = config["PARAMS_FOR_ANALYSIS"]["SOURCES"],
                                    day_segment = config["PARAMS_FOR_ANALYSIS"]["DAY_SEGMENTS"]),
                                zip,
                                rows_nan_threshold = config["PARAMS_FOR_ANALYSIS"]["ROWS_NAN_THRESHOLD"],
                                cols_nan_threshold = config["PARAMS_FOR_ANALYSIS"]["COLS_NAN_THRESHOLD"]))
    files_to_compute.extend(expand(
                                expand("data/processed/data_for_population_model/{min_valid_hours_per_day}hours_{min_valid_bins_per_hour}bins/{{rows_nan_threshold}}|{{cols_nan_threshold}}_{days_before_threshold}|{days_after_threshold}_{cols_var_threshold}/{source}_{day_segment}_clean.csv",
                                    min_valid_hours_per_day=config["OVERALL_COMPLIANCE_HEATMAP"]["MIN_VALID_HOURS_PER_DAY"],
                                    min_valid_bins_per_hour=config["PHONE_VALID_SENSED_DAYS"]["MIN_VALID_BINS_PER_HOUR"],
                                    days_before_threshold = config["PARAMS_FOR_ANALYSIS"]["PARTICIPANT_DAYS_BEFORE_THRESHOLD"],
                                    days_after_threshold = config["PARAMS_FOR_ANALYSIS"]["PARTICIPANT_DAYS_AFTER_THRESHOLD"],
                                    cols_var_threshold = config["PARAMS_FOR_ANALYSIS"]["COLS_VAR_THRESHOLD"],
                                    source = config["PARAMS_FOR_ANALYSIS"]["SOURCES"],
                                    day_segment = config["PARAMS_FOR_ANALYSIS"]["DAY_SEGMENTS"]),
                                zip,
                                rows_nan_threshold = config["PARAMS_FOR_ANALYSIS"]["ROWS_NAN_THRESHOLD"],
                                cols_nan_threshold = config["PARAMS_FOR_ANALYSIS"]["COLS_NAN_THRESHOLD"]))
    files_to_compute.extend(expand("data/processed/data_for_population_model/demographic_features.csv"))
    files_to_compute.extend(expand("data/processed/data_for_population_model/targets_{summarised}.csv",
                                summarised = config["PARAMS_FOR_ANALYSIS"]["SUMMARISED"]))
    files_to_compute.extend(expand(
                                expand("data/processed/data_for_population_model/{{rows_nan_threshold}}|{{cols_nan_threshold}}_{days_before_threshold}|{days_after_threshold}_{cols_var_threshold}/{source}_{day_segment}_nancellsratio.csv",
                                    days_before_threshold = config["PARAMS_FOR_ANALYSIS"]["PARTICIPANT_DAYS_BEFORE_THRESHOLD"],
                                    days_after_threshold = config["PARAMS_FOR_ANALYSIS"]["PARTICIPANT_DAYS_AFTER_THRESHOLD"],
                                    cols_var_threshold = config["PARAMS_FOR_ANALYSIS"]["COLS_VAR_THRESHOLD"],
                                    source = config["PARAMS_FOR_ANALYSIS"]["SOURCES"],
                                    day_segment = config["PARAMS_FOR_ANALYSIS"]["DAY_SEGMENTS"]),
                                zip,
                                rows_nan_threshold = config["PARAMS_FOR_ANALYSIS"]["ROWS_NAN_THRESHOLD"],
                                cols_nan_threshold = config["PARAMS_FOR_ANALYSIS"]["COLS_NAN_THRESHOLD"]))
    files_to_compute.extend(expand(
                                expand("data/processed/data_for_population_model/{{rows_nan_threshold}}|{{cols_nan_threshold}}_{days_before_threshold}|{days_after_threshold}_{cols_var_threshold}/{source}_{day_segment}_{summarised}.csv",
                                    days_before_threshold = config["PARAMS_FOR_ANALYSIS"]["PARTICIPANT_DAYS_BEFORE_THRESHOLD"],
                                    days_after_threshold = config["PARAMS_FOR_ANALYSIS"]["PARTICIPANT_DAYS_AFTER_THRESHOLD"],
                                    cols_var_threshold = config["PARAMS_FOR_ANALYSIS"]["COLS_VAR_THRESHOLD"],
                                    source = config["PARAMS_FOR_ANALYSIS"]["SOURCES"],
                                    day_segment = config["PARAMS_FOR_ANALYSIS"]["DAY_SEGMENTS"],
                                    summarised = config["PARAMS_FOR_ANALYSIS"]["SUMMARISED"]),
                                zip,
                                rows_nan_threshold = config["PARAMS_FOR_ANALYSIS"]["ROWS_NAN_THRESHOLD"],
                                cols_nan_threshold = config["PARAMS_FOR_ANALYSIS"]["COLS_NAN_THRESHOLD"]))
    files_to_compute.extend(expand(
                                expand("data/processed/output_population_model/{{rows_nan_threshold}}|{{cols_nan_threshold}}_{days_before_threshold}|{days_after_threshold}_{cols_var_threshold}/{source}_{day_segment}_{summarised}_{cv_method}_baseline.csv",
                                    days_before_threshold = config["PARAMS_FOR_ANALYSIS"]["PARTICIPANT_DAYS_BEFORE_THRESHOLD"],
                                    days_after_threshold = config["PARAMS_FOR_ANALYSIS"]["PARTICIPANT_DAYS_AFTER_THRESHOLD"],
                                    cols_var_threshold = config["PARAMS_FOR_ANALYSIS"]["COLS_VAR_THRESHOLD"],
                                    cv_method = config["PARAMS_FOR_ANALYSIS"]["CV_METHODS"],
                                    source = config["PARAMS_FOR_ANALYSIS"]["SOURCES"],
                                    day_segment = config["PARAMS_FOR_ANALYSIS"]["DAY_SEGMENTS"],
                                    summarised = config["PARAMS_FOR_ANALYSIS"]["SUMMARISED"]),
                                zip,
                                rows_nan_threshold = config["PARAMS_FOR_ANALYSIS"]["ROWS_NAN_THRESHOLD"],
                                cols_nan_threshold = config["PARAMS_FOR_ANALYSIS"]["COLS_NAN_THRESHOLD"]))
    files_to_compute.extend(expand(
                                expand("data/processed/output_population_model/{{rows_nan_threshold}}|{{cols_nan_threshold}}_{days_before_threshold}|{days_after_threshold}_{cols_var_threshold}/{{model}}/{cv_method}/{source}_{day_segment}_{summarised}_{{scaler}}/{result}.csv",
                                    days_before_threshold = config["PARAMS_FOR_ANALYSIS"]["PARTICIPANT_DAYS_BEFORE_THRESHOLD"],
                                    days_after_threshold = config["PARAMS_FOR_ANALYSIS"]["PARTICIPANT_DAYS_AFTER_THRESHOLD"],
                                    cols_var_threshold = config["PARAMS_FOR_ANALYSIS"]["COLS_VAR_THRESHOLD"],
                                    cv_method = config["PARAMS_FOR_ANALYSIS"]["CV_METHODS"],
                                    source = config["PARAMS_FOR_ANALYSIS"]["SOURCES"],
                                    day_segment = config["PARAMS_FOR_ANALYSIS"]["DAY_SEGMENTS"],
                                    summarised = config["PARAMS_FOR_ANALYSIS"]["SUMMARISED"],
                                    result = results), 
                                zip,
                                rows_nan_threshold = rows_nan_thresholds,
                                cols_nan_threshold = cols_nan_thresholds,
                                model = models,
                                scaler = scalers))

rule all:
    input:
        files_to_compute

rule clean:
    shell:
        "rm -rf data/raw/* && rm -rf data/interim/* && rm -rf data/processed/* && rm -rf reports/figures/* && rm -rf reports/*.zip && rm -rf reports/compliance/*"