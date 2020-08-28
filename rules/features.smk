rule join_features_from_providers:
    input:
        location_features = find_features_files
    output:
        "data/processed/features/{pid}/{sensor_key}.csv"
    script:
        "../src/features/join_features_from_providers.R"

rule messages_features:
    input: 
        expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["MESSAGES"]["DB_TABLE"]),
        day_segments_labels = expand("data/interim/{sensor}_day_segments_labels.csv", sensor=config["MESSAGES"]["DB_TABLE"])
    params:
        messages_type = "{messages_type}",
        features = lambda wildcards: config["MESSAGES"]["FEATURES"][wildcards.messages_type]
    output:
        "data/processed/{pid}/messages_{messages_type}.csv"
    script:
        "../src/features/messages_features.R"

rule call_features:
    input: 
        expand("data/raw/{{pid}}/{sensor}_with_datetime_unified.csv", sensor=config["CALLS"]["DB_TABLE"]),
        day_segments_labels = expand("data/interim/{sensor}_day_segments_labels.csv", sensor=config["CALLS"]["DB_TABLE"])
    params:
        call_type = "{call_type}",
        features = lambda wildcards: config["CALLS"]["FEATURES"][wildcards.call_type]
    output:
        "data/processed/{pid}/calls_{call_type}.csv"
    script:
        "../src/features/call_features.R"

rule battery_deltas:
    input:
        expand("data/raw/{{pid}}/{sensor}_with_datetime_unified.csv", sensor=config["BATTERY"]["DB_TABLE"])
    output:
        "data/processed/{pid}/battery_deltas.csv"
    script:
        "../src/features/battery_deltas.R"

rule screen_deltas:
    input:
        screen = expand("data/raw/{{pid}}/{sensor}_with_datetime_unified.csv", sensor=config["SCREEN"]["DB_TABLE"])
    output:
        "data/processed/{pid}/screen_deltas.csv"
    script:
        "../src/features/screen_deltas.R"

rule google_activity_recognition_deltas:
    input:
        expand("data/raw/{{pid}}/{sensor}_with_datetime_unified.csv", sensor=config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["ANDROID"])
    output:
        expand("data/processed/{{pid}}/{sensor}_deltas.csv", sensor=config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["ANDROID"])
    script:
        "../src/features/activity_recognition_deltas.R"

rule ios_activity_recognition_deltas:
    input:
        expand("data/raw/{{pid}}/{sensor}_with_datetime_unified.csv", sensor=config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["IOS"])
    output:
        expand("data/processed/{{pid}}/{sensor}_deltas.csv", sensor=config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["IOS"])
    script:
        "../src/features/activity_recognition_deltas.R"

rule locations_python_features:
    input:
        location_data = expand("data/raw/{{pid}}/{sensor}_processed_{locations_to_use}.csv", sensor=config["LOCATIONS"]["DB_TABLE"], locations_to_use=config["LOCATIONS"]["LOCATIONS_TO_USE"]),
        day_segments_labels = "data/interim/day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["LOCATIONS"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
    output:
        "data/interim/{pid}/locations_features/locations_python_{provider_key}.csv"
    script:
        "../src/features/location/locations_entry.py"

rule locations_r_features:
    input:
        location_data = expand("data/raw/{{pid}}/{sensor}_processed_{locations_to_use}.csv", sensor=config["LOCATIONS"]["DB_TABLE"], locations_to_use=config["LOCATIONS"]["LOCATIONS_TO_USE"]),
        day_segments_labels = "data/interim/day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["LOCATIONS"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}"
    output:
        "data/interim/{pid}/locations_features/locations_r_{provider_key}.csv"
    script:
        "../src/features/location/locations_entry.R"

rule bluetooth_features:
    input: 
        expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["BLUETOOTH"]["DB_TABLE"]),
        day_segments = expand("data/interim/{sensor}_day_segments.csv", sensor=config["BLUETOOTH"]["DB_TABLE"])
    params:
        features = config["BLUETOOTH"]["FEATURES"]
    output:
        "data/processed/{pid}/bluetooth_features.csv"
    script:
        "../src/features/bluetooth_features.R"

rule activity_features:
    input:
        optional_ar_input
    params:
        segment = "{day_segment}",
        features = config["ACTIVITY_RECOGNITION"]["FEATURES"]
    output:
        "data/processed/{pid}/activity_recognition_{day_segment}.csv"
    script:
        "../src/features/activity_recognition.py"

rule battery_features:
    input:
        "data/processed/{pid}/battery_deltas.csv"
    params:
        day_segment = "{day_segment}",
        features = config["BATTERY"]["FEATURES"]
    output:
        "data/processed/{pid}/battery_{day_segment}.csv"
    script:
        "../src/features/battery_features.py"

rule screen_features:
    input:
        screen_deltas = "data/processed/{pid}/screen_deltas.csv",
        phone_sensed_bins = "data/interim/{pid}/phone_sensed_bins.csv"
    params:
        day_segment = "{day_segment}",
        reference_hour_first_use = config["SCREEN"]["REFERENCE_HOUR_FIRST_USE"],
        features_deltas = config["SCREEN"]["FEATURES_DELTAS"],
        episode_types = config["SCREEN"]["EPISODE_TYPES"],
        ignore_episodes_shorter_than = config["SCREEN"]["IGNORE_EPISODES_SHORTER_THAN"],
        ignore_episodes_longer_than = config["SCREEN"]["IGNORE_EPISODES_LONGER_THAN"],
        bin_size = config["PHONE_VALID_SENSED_BINS"]["BIN_SIZE"]
    output:
        "data/processed/{pid}/screen_{day_segment}.csv"
    script:
        "../src/features/screen_features.py"

rule light_features:
    input:
        expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["LIGHT"]["DB_TABLE"]),
    params:
        day_segment = "{day_segment}",
        features = config["LIGHT"]["FEATURES"],
    output:
        "data/processed/{pid}/light_{day_segment}.csv"
    script:
        "../src/features/light_features.py"

rule conversation_features:
    input:
        optional_conversation_input
    params:
        day_segment = "{day_segment}",
        features = config["CONVERSATION"]["FEATURES"],
        recordingMinutes = config["CONVERSATION"]["RECORDINGMINUTES"],
        pausedMinutes = config["CONVERSATION"]["PAUSEDMINUTES"],
    output:
        "data/processed/{pid}/conversation_{day_segment}.csv"
    script:
        "../src/features/conversation_features.py"

rule accelerometer_features:
    input:
        expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["ACCELEROMETER"]["DB_TABLE"]),
    params:
        day_segment = "{day_segment}",
        magnitude = config["ACCELEROMETER"]["FEATURES"]["MAGNITUDE"],
        exertional_activity_episode = config["ACCELEROMETER"]["FEATURES"]["EXERTIONAL_ACTIVITY_EPISODE"],
        nonexertional_activity_episode = config["ACCELEROMETER"]["FEATURES"]["NONEXERTIONAL_ACTIVITY_EPISODE"],
        valid_sensed_minutes = config["ACCELEROMETER"]["FEATURES"]["VALID_SENSED_MINUTES"],
    output:
        "data/processed/{pid}/accelerometer_{day_segment}.csv"
    script:
        "../src/features/accelerometer_features.py"

rule applications_foreground_features:
    input:
        expand("data/interim/{{pid}}/{sensor}_with_datetime_with_genre.csv", sensor=config["APPLICATIONS_FOREGROUND"]["DB_TABLE"])
    params:
        day_segment = "{day_segment}",
        single_categories = config["APPLICATIONS_FOREGROUND"]["SINGLE_CATEGORIES"],
        multiple_categories = config["APPLICATIONS_FOREGROUND"]["MULTIPLE_CATEGORIES"],
        single_apps = config["APPLICATIONS_FOREGROUND"]["SINGLE_APPS"],
        excluded_categories = config["APPLICATIONS_FOREGROUND"]["EXCLUDED_CATEGORIES"],
        excluded_apps = config["APPLICATIONS_FOREGROUND"]["EXCLUDED_APPS"],
        features = config["APPLICATIONS_FOREGROUND"]["FEATURES"],
    output:
        "data/processed/{pid}/applications_foreground_{day_segment}.csv"
    script:
        "../src/features/applications_foreground_features.py"

rule wifi_features:
    input: 
        expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["WIFI"]["DB_TABLE"]),
        day_segments = expand("data/interim/{sensor}_day_segments.csv", sensor=config["WIFI"]["DB_TABLE"])
    params:
        features = config["WIFI"]["FEATURES"]
    output:
        "data/processed/{pid}/wifi_features.csv"
    script:
        "../src/features/wifi_features.R"

rule fitbit_heartrate_features:
    input:
        heartrate_summary_data = "data/raw/{pid}/fitbit_heartrate_summary_with_datetime.csv",
        heartrate_intraday_data = "data/raw/{pid}/fitbit_heartrate_intraday_with_datetime.csv"
    params:
        day_segment = "{day_segment}",
        summary_features = config["HEARTRATE"]["SUMMARY_FEATURES"],
        intraday_features = config["HEARTRATE"]["INTRADAY_FEATURES"]
    output:
        "data/processed/{pid}/fitbit_heartrate_{day_segment}.csv"
    script:
        "../src/features/fitbit_heartrate_features.py"

rule fitbit_step_features:
    input:
        step_data = "data/raw/{pid}/fitbit_step_intraday_with_datetime.csv",
        sleep_data = optional_steps_sleep_input
    params:
        day_segment = "{day_segment}",
        features_all_steps = config["STEP"]["FEATURES"]["ALL_STEPS"],
        features_sedentary_bout = config["STEP"]["FEATURES"]["SEDENTARY_BOUT"],
        features_active_bout = config["STEP"]["FEATURES"]["ACTIVE_BOUT"],
        threshold_active_bout = config["STEP"]["THRESHOLD_ACTIVE_BOUT"],
        include_zero_step_rows = config["STEP"]["INCLUDE_ZERO_STEP_ROWS"],
        exclude_sleep = config["STEP"]["EXCLUDE_SLEEP"]["EXCLUDE"],
        exclude_sleep_type = config["STEP"]["EXCLUDE_SLEEP"]["TYPE"],
        exclude_sleep_fixed_start = config["STEP"]["EXCLUDE_SLEEP"]["FIXED"]["START"],
        exclude_sleep_fixed_end = config["STEP"]["EXCLUDE_SLEEP"]["FIXED"]["END"],
    output:
        "data/processed/{pid}/fitbit_step_{day_segment}.csv"
    script:
        "../src/features/fitbit_step_features.py"

rule fitbit_sleep_features:
    input:
        sleep_summary_data = "data/raw/{pid}/fitbit_sleep_summary_with_datetime.csv",
        sleep_intraday_data = "data/raw/{pid}/fitbit_sleep_intraday_with_datetime.csv"
    params:
        day_segment = "{day_segment}",
        summary_features = config["SLEEP"]["SUMMARY_FEATURES"],
        sleep_types = config["SLEEP"]["SLEEP_TYPES"]
    output:
        "data/processed/{pid}/fitbit_sleep_{day_segment}.csv"
    script:
        "../src/features/fitbit_sleep_features.py"
