rule restore_sql_file:
    input:
        sql_file = "data/external/rapids_example.sql",
        db_credentials = ".env"
    params:
        group = config["DOWNLOAD_PARTICIPANTS"]["GROUP"]
    output:
        touch("data/interim/restore_sql_file.done")
    script:
        "../src/data/restore_sql_file.py"

rule create_example_participant_files:
    output:
        expand("data/external/{pid}", pid = ["example01", "example02"])
    shell:
        "echo 'a748ee1a-1d0b-4ae9-9074-279a2b6ba524\nandroid\ntest01\n2020/04/23,2020/05/04\n' >> ./data/external/example01 && echo '13dbc8a3-dae3-4834-823a-4bc96a7d459d\nios\ntest02\n2020/04/23,2020/05/04\n' >> ./data/external/example02"

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
        timezone = config["TIMEZONE"],
        aware_multiplatform_tables = config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["ANDROID"] + "," + config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["IOS"] + "," + config["CONVERSATION"]["DB_TABLE"]["ANDROID"] + "," + config["CONVERSATION"]["DB_TABLE"]["IOS"],
        unifiable_sensors = {"calls": config["CALLS"]["DB_TABLE"], "battery": config["BATTERY"]["DB_TABLE"], "screen": config["SCREEN"]["DB_TABLE"], "ios_activity_recognition": config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["IOS"], "ios_conversation": config["CONVERSATION"]["DB_TABLE"]["IOS"]}
    output:
        "data/raw/{pid}/{sensor}_raw.csv"
    script:
        "../src/data/download_dataset.R"

rule compute_day_segments:
    input: 
        config["DAY_SEGMENTS"]["FILE"]
    params:
        day_segments_type = config["DAY_SEGMENTS"]["TYPE"],
        pid = "{pid}"
    output:
        segments_file = "data/interim/day_segments/{pid}_day_segments.csv",
        segments_labels_file = "data/interim/day_segments/{pid}_day_segments_labels.csv",
    script:
        "../src/data/compute_day_segments.py"

PHONE_SENSORS = []
PHONE_SENSORS.extend([config["MESSAGES"]["DB_TABLE"], config["CALLS"]["DB_TABLE"], config["LOCATIONS"]["DB_TABLE"], config["BLUETOOTH"]["DB_TABLE"], config["BATTERY"]["DB_TABLE"], config["SCREEN"]["DB_TABLE"], config["LIGHT"]["DB_TABLE"], config["ACCELEROMETER"]["DB_TABLE"], config["APPLICATIONS_FOREGROUND"]["DB_TABLE"], config["CONVERSATION"]["DB_TABLE"]["ANDROID"], config["CONVERSATION"]["DB_TABLE"]["IOS"], config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["ANDROID"], config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["IOS"]])
PHONE_SENSORS.extend(config["PHONE_VALID_SENSED_BINS"]["DB_TABLES"])

if len(config["WIFI"]["DB_TABLE"]["VISIBLE_ACCESS_POINTS"]) > 0:
    PHONE_SENSORS.append(config["WIFI"]["DB_TABLE"]["VISIBLE_ACCESS_POINTS"])
if len(config["WIFI"]["DB_TABLE"]["CONNECTED_ACCESS_POINTS"]) > 0:
    PHONE_SENSORS.append(config["WIFI"]["DB_TABLE"]["CONNECTED_ACCESS_POINTS"])


rule readable_datetime:
    input:
        sensor_input = "data/raw/{pid}/{sensor}_raw.csv",
        day_segments = "data/interim/day_segments/{pid}_day_segments.csv"
    params:
        timezones = None,
        fixed_timezone = config["READABLE_DATETIME"]["FIXED_TIMEZONE"],
        day_segments_type = config["DAY_SEGMENTS"]["TYPE"]
    wildcard_constraints:
        sensor = '.*(' + '|'.join([re.escape(x) for x in PHONE_SENSORS]) + ').*' # only process smartphone sensors, not fitbit
    output:
        "data/raw/{pid}/{sensor}_with_datetime.csv"
    script:
        "../src/data/readable_datetime.R"

rule phone_sensed_bins:
    input:
        all_sensors =  optional_phone_sensed_bins_input
    params:
        bin_size = config["PHONE_VALID_SENSED_BINS"]["BIN_SIZE"]
    output:
        "data/interim/{pid}/phone_sensed_bins.csv"
    script:
        "../src/data/phone_sensed_bins.R"

rule phone_valid_sensed_days:
    input:
        phone_sensed_bins =  "data/interim/{pid}/phone_sensed_bins.csv"
    params:
        min_valid_hours_per_day = "{min_valid_hours_per_day}",
        min_valid_bins_per_hour = "{min_valid_bins_per_hour}"
    output:
        "data/interim/{pid}/phone_valid_sensed_days_{min_valid_hours_per_day}hours_{min_valid_bins_per_hour}bins.csv"
    script:
        "../src/data/phone_valid_sensed_days.R"


rule unify_ios_android:
    input:
        sensor_data = "data/raw/{pid}/{sensor}_with_datetime.csv",
        participant_info = "data/external/{pid}"
    params:
        sensor = "{sensor}",
        unifiable_sensors = {"calls": config["CALLS"]["DB_TABLE"], "battery": config["BATTERY"]["DB_TABLE"], "screen": config["SCREEN"]["DB_TABLE"], "ios_activity_recognition": config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["IOS"], "ios_conversation": config["CONVERSATION"]["DB_TABLE"]["IOS"]}
    output:
        "data/raw/{pid}/{sensor}_with_datetime_unified.csv"
    script:
        "../src/data/unify_ios_android.R"

rule process_location_types:
    input:
        locations = "data/raw/{pid}/{sensor}_with_datetime.csv",
        phone_sensed_bins = rules.phone_sensed_bins.output,
        day_segments = "data/interim/day_segments/{pid}_day_segments.csv"
    params:
        bin_size = config["PHONE_VALID_SENSED_BINS"]["BIN_SIZE"],
        timezone = config["LOCATIONS"]["TIMEZONE"],
        consecutive_threshold = config["LOCATIONS"]["FUSED_RESAMPLED_CONSECUTIVE_THRESHOLD"],
        time_since_valid_location = config["LOCATIONS"]["FUSED_RESAMPLED_TIME_SINCE_VALID_LOCATION"],
        day_segments_type = config["DAY_SEGMENTS"]["TYPE"],
        locations_to_use = "{locations_to_used}"
    output:
        "data/raw/{pid}/{sensor}_processed_{locations_to_used}.csv"
    script:
        "../src/data/process_location_types.R"

rule application_genres:
    input:
        "data/raw/{pid}/{sensor}_with_datetime.csv"
    params:
        catalogue_source = config["APPLICATION_GENRES"]["CATALOGUE_SOURCE"],
        catalogue_file = config["APPLICATION_GENRES"]["CATALOGUE_FILE"],
        update_catalogue_file = config["APPLICATION_GENRES"]["UPDATE_CATALOGUE_FILE"],
        scrape_missing_genres = config["APPLICATION_GENRES"]["SCRAPE_MISSING_GENRES"]
    output:
        "data/raw/{pid}/{sensor}_with_datetime_with_genre.csv"
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

rule join_wifi_tables:
    input: 
        unpack(optional_wifi_input)
    output:
        "data/raw/{pid}/wifi_with_datetime_visibleandconnected.csv"
    script:
        "../src/data/join_visible_and_connected_wifi.R"