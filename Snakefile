configfile: "config.yaml"
include: "rules/renv.snakefile"
include: "rules/preprocessing.snakefile"
include: "rules/features.snakefile"
include: "rules/models.snakefile"
include: "rules/reports.snakefile"
include: "rules/mystudy.snakefile" # You can add snakfiles with rules tailored to your project

files_to_compute = []

if len(config["PIDS"]) == 0:
    raise ValueError("Add participants IDs to PIDS in config.yaml. Remember to create their participant files in data/external")

if config["MESSAGES"]["COMPUTE"]:
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["MESSAGES"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["MESSAGES"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/processed/{pid}/messages_{messages_type}_{day_segment}.csv", pid=config["PIDS"], messages_type = config["MESSAGES"]["TYPES"], day_segment = config["MESSAGES"]["DAY_SEGMENTS"]))

if config["CALLS"]["COMPUTE"]:
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["CALLS"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["CALLS"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime_unified.csv", pid=config["PIDS"], sensor=config["CALLS"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/processed/{pid}/calls_{call_type}_{segment}.csv", pid=config["PIDS"], call_type=config["CALLS"]["TYPES"], segment = config["CALLS"]["DAY_SEGMENTS"]))

if config["BARNETT_LOCATION"]["COMPUTE"]:
    # TODO add files_to_compute.extend(optional_location_input(None))
    if config["BARNETT_LOCATION"]["LOCATIONS_TO_USE"] == "RESAMPLE_FUSED" and config["BARNETT_LOCATION"]["DB_TABLE"] not in config["TABLES_FOR_SENSED_BINS"]:
        raise ValueError("Error: Add your locations table (and as many sensor tables as you have) to TABLES_FOR_SENSED_BINS in config.yaml. This is necessary to compute phone_sensed_bins (bins of time when the smartphone was sensing data) which is used to resample fused location data (RESAMPLED_FUSED)")
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["BARNETT_LOCATION"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["BARNETT_LOCATION"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/processed/{pid}/location_barnett_{segment}.csv", pid=config["PIDS"], segment = config["BARNETT_LOCATION"]["DAY_SEGMENTS"]))

if config["BLUETOOTH"]["COMPUTE"]:
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["BLUETOOTH"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["BLUETOOTH"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/processed/{pid}/bluetooth_{segment}.csv", pid=config["PIDS"], segment = config["BLUETOOTH"]["DAY_SEGMENTS"]))

if config["ACTIVITY_RECOGNITION"]["COMPUTE"]:
    # TODO add files_to_compute.extend(optional_ar_input(None)), the Android or iOS table gets processed depending on each participant
    files_to_compute.extend(expand("data/processed/{pid}/activity_recognition_{segment}.csv",pid=config["PIDS"], segment = config["ACTIVITY_RECOGNITION"]["DAY_SEGMENTS"]))

if config["BATTERY"]["COMPUTE"]:
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["BATTERY"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["BATTERY"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime_unified.csv", pid=config["PIDS"], sensor=config["BATTERY"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/processed/{pid}/battery_deltas.csv", pid=config["PIDS"]))
    files_to_compute.extend(expand("data/processed/{pid}/battery_{day_segment}.csv", pid = config["PIDS"], day_segment = config["BATTERY"]["DAY_SEGMENTS"]))

if config["SCREEN"]["COMPUTE"]:
    if config["SCREEN"]["DB_TABLE"] not in config["TABLES_FOR_SENSED_BINS"]:
        raise ValueError("Error: Add your screen table (and as many sensor tables as you have) to TABLES_FOR_SENSED_BINS in config.yaml. This is necessary to compute phone_sensed_bins (bins of time when the smartphone was sensing data)")
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["SCREEN"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["SCREEN"]["DB_TABLE"]))
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
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["WIFI"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["WIFI"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/processed/{pid}/wifi_{day_segment}.csv", pid = config["PIDS"], day_segment = config["WIFI"]["DAY_SEGMENTS"]))

if config["HEARTRATE"]["COMPUTE"]:
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["HEARTRATE"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/fitbit_heartrate_{fitbit_data_type}_with_datetime.csv", pid=config["PIDS"], fitbit_data_type=["summary", "intraday"]))
    files_to_compute.extend(expand("data/processed/{pid}/fitbit_heartrate_{day_segment}.csv", pid = config["PIDS"], day_segment = config["HEARTRATE"]["DAY_SEGMENTS"]))

if config["STEP"]["COMPUTE"]:
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["STEP"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/fitbit_step_{fitbit_data_type}_with_datetime.csv", pid=config["PIDS"], fitbit_data_type=["intraday"]))
    files_to_compute.extend(expand("data/processed/{pid}/fitbit_step_{day_segment}.csv", pid = config["PIDS"], day_segment = config["STEP"]["DAY_SEGMENTS"]))

if config["SLEEP"]["COMPUTE"]:
    files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["SLEEP"]["DB_TABLE"]))
    files_to_compute.extend(expand("data/raw/{pid}/fitbit_sleep_{fitbit_data_type}_with_datetime.csv", pid=config["PIDS"], fitbit_data_type=["intraday"]))
    files_to_compute.extend(expand("data/processed/{pid}/fitbit_sleep_{day_segment}.csv", pid = config["PIDS"], day_segment = config["SLEEP"]["DAY_SEGMENTS"]))

if config["CONVERSATION"]["COMPUTE"]:
    # TODO add files_to_compute.extend(optional_conversation_input(None)), the Android or iOS table gets processed depending on each participant
    files_to_compute.extend(expand("data/processed/{pid}/conversation_{segment}.csv",pid=config["PIDS"], segment = config["CONVERSATION"]["DAY_SEGMENTS"]))

rule all:
    input:
        files_to_compute

rule clean:
    shell:
        "rm -rf data/raw/* && rm -rf data/interim/* && rm -rf data/processed/* && rm -rf reports/figures/* && rm -rf reports/*.zip && rm -rf reports/compliance/*"