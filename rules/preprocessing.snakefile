rule download_participants:
    params:
        group = config["DOWNLOAD_PARTICIPANTS"]["GROUP"],
        ignored_device_ids = config["DOWNLOAD_PARTICIPANTS"]["IGNORED_DEVICE_IDS"],
        timezone = config["TIMEZONE"]
    priority: 1
    script:
        "../src/data/download_participants.R"

rule download_dataset:
    input:
        "data/external/{pid}"
    params:
        group = config["DOWNLOAD_DATASET"]["GROUP"],
        table = "{sensor}",
        timezone = config["TIMEZONE"]
    output:
        "data/raw/{pid}/{sensor}_raw.csv"
    script:
        "../src/data/download_dataset.R"

PHONE_SENSORS = []
PHONE_SENSORS.extend([config["MESSAGES"]["DB_TABLE"], config["CALLS"]["DB_TABLE"], config["BARNETT_LOCATION"]["DB_TABLE"], config["BLUETOOTH"]["DB_TABLE"], config["BATTERY"]["DB_TABLE"], config["SCREEN"]["DB_TABLE"], config["LIGHT"]["DB_TABLE"], config["ACCELEROMETER"]["DB_TABLE"], config["APPLICATIONS_FOREGROUND"]["DB_TABLE"],config["WIFI"]["DB_TABLE"], config["CONVERSATION"]["DB_TABLE"]["ANDROID"], config["CONVERSATION"]["DB_TABLE"]["IOS"], config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["ANDROID"], config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["IOS"]])

rule readable_datetime:
    input:
        sensor_input = rules.download_dataset.output
    params:
        timezones = None,
        fixed_timezone = config["READABLE_DATETIME"]["FIXED_TIMEZONE"]
    wildcard_constraints:
        sensor = '.*(' + '|'.join([re.escape(x) for x in PHONE_SENSORS]) + ').*' # only process smartphone sensors, not fitbit
    output:
        "data/raw/{pid}/{sensor}_with_datetime.csv"
    script:
        "../src/data/readable_datetime.R"

rule phone_valid_sensed_days:
    input:
        all_sensors =  expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["TABLES_FOR_SENSED_BINS"])
    params:
        bin_size = config["PHONE_VALID_SENSED_DAYS"]["BIN_SIZE"],
        min_valid_hours = config["PHONE_VALID_SENSED_DAYS"]["MIN_VALID_HOURS"],
        min_bins_per_hour = config["PHONE_VALID_SENSED_DAYS"]["MIN_BINS_PER_HOUR"]
    output:
        "data/interim/{pid}/phone_valid_sensed_days.csv"
    script:
        "../src/data/phone_valid_sensed_days.R"

rule phone_sensed_bins:
    input:
        all_sensors =  expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["TABLES_FOR_SENSED_BINS"])
    params:
        bin_size = config["PHONE_VALID_SENSED_DAYS"]["BIN_SIZE"]
    output:
        "data/interim/{pid}/phone_sensed_bins.csv"
    script:
        "../src/data/phone_sensed_bins.R"

rule unify_ios_android:
    input:
        sensor_data = "data/raw/{pid}/{sensor}_with_datetime.csv",
        participant_info = "data/external/{pid}"
    params:
        sensor = "{sensor}"
    output:
        "data/raw/{pid}/{sensor}_with_datetime_unified.csv"
    script:
        "../src/data/unify_ios_android.R"

rule resample_fused_location:
    input:
        locations = "data/raw/{pid}/{sensor}_raw.csv",
        phone_sensed_bins = rules.phone_sensed_bins.output
    params:
        bin_size = config["PHONE_VALID_SENSED_DAYS"]["BIN_SIZE"],
        timezone = config["RESAMPLE_FUSED_LOCATION"]["TIMEZONE"],
        consecutive_threshold = config["RESAMPLE_FUSED_LOCATION"]["CONSECUTIVE_THRESHOLD"],
        time_since_valid_location = config["RESAMPLE_FUSED_LOCATION"]["TIME_SINCE_VALID_LOCATION"]
    output:
        "data/raw/{pid}/{sensor}_resampled.csv"
    script:
        "../src/data/resample_fused_location.R"

rule application_genres:
    input:
        "data/raw/{pid}/{sensor}_with_datetime.csv"
    params:
        catalogue_source = config["APPLICATION_GENRES"]["CATALOGUE_SOURCE"],
        catalogue_file = config["APPLICATION_GENRES"]["CATALOGUE_FILE"],
        update_catalogue_file = config["APPLICATION_GENRES"]["UPDATE_CATALOGUE_FILE"],
        scrape_missing_genres = config["APPLICATION_GENRES"]["SCRAPE_MISSING_GENRES"]
    output:
        "data/interim/{pid}/{sensor}_with_datetime_with_genre.csv"
    script:
        "../src/data/application_genres.R"

rule fitbit_heartrate_with_datetime:
    input:
        expand("data/raw/{{pid}}/{fitbit_table}_raw.csv", fitbit_table=config["HEARTRATE"]["DB_TABLE"])
    params:
        local_timezone = config["READABLE_DATETIME"]["FIXED_TIMEZONE"],
        fitbit_sensor = "heartrate"
    output:
        summary_data = "data/raw/{pid}/fitbit_heartrate_summary_with_datetime.csv",
        intraday_data = "data/raw/{pid}/fitbit_heartrate_intraday_with_datetime.csv"
    script:
        "../src/data/fitbit_readable_datetime.py"

rule fitbit_step_with_datetime:
    input:
        expand("data/raw/{{pid}}/{fitbit_table}_raw.csv", fitbit_table=config["STEP"]["DB_TABLE"])
    params:
        local_timezone = config["READABLE_DATETIME"]["FIXED_TIMEZONE"],
        fitbit_sensor = "steps"
    output:
        intraday_data = "data/raw/{pid}/fitbit_step_intraday_with_datetime.csv"
    script:
        "../src/data/fitbit_readable_datetime.py"

rule fitbit_sleep_with_datetime:
    input:
        expand("data/raw/{{pid}}/{fitbit_table}_raw.csv", fitbit_table=config["SLEEP"]["DB_TABLE"])
    params:
        local_timezone = config["READABLE_DATETIME"]["FIXED_TIMEZONE"],
        fitbit_sensor = "sleep"
    output:
        summary_data = "data/raw/{pid}/fitbit_sleep_summary_with_datetime.csv",
        intraday_data = "data/raw/{pid}/fitbit_sleep_intraday_with_datetime.csv"
    script:
        "../src/data/fitbit_readable_datetime.py"
