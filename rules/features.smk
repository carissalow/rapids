rule join_features_from_providers:
    input:
        location_features = find_features_files
    output:
        "data/processed/features/{pid}/{sensor_key}.csv"
    script:
        "../src/features/join_features_from_providers.R"

rule messages_r_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["MESSAGES"]["DB_TABLE"]),
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["MESSAGES"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}"
    output:
        "data/interim/{pid}/messages_features/messages_r_{provider_key}.csv"
    script:
        "../src/features/messages/messages_entry.R"

rule messages_python_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["MESSAGES"]["DB_TABLE"]),
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["MESSAGES"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}"
    output:
        "data/interim/{pid}/messages_features/messages_python_{provider_key}.csv"
    script:
        "../src/features/messages/messages_entry.py"

rule calls_python_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor}_with_datetime_unified.csv", sensor=config["CALLS"]["DB_TABLE"]),
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["CALLS"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}"
    output:
        "data/interim/{pid}/calls_features/calls_python_{provider_key}.csv"
    script:
        "../src/features/calls/calls_entry.py"

rule calls_r_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor}_with_datetime_unified.csv", sensor=config["CALLS"]["DB_TABLE"]),
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["CALLS"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}"
    output:
        "data/interim/{pid}/calls_features/calls_r_{provider_key}.csv"
    script:
        "../src/features/calls/calls_entry.R"

rule battery_episodes:
    input:
        expand("data/raw/{{pid}}/{sensor}_raw.csv", sensor=config["BATTERY"]["DB_TABLE"])
    output:
        "data/interim/{pid}/battery_episodes.csv"
    script:
        "../src/features/battery/episodes/battery_episodes.R"

rule screen_episodes:
    input:
        screen = expand("data/raw/{{pid}}/{sensor}_with_datetime_unified.csv", sensor=config["SCREEN"]["DB_TABLE"])
    output:
        "data/interim/{pid}/screen_episodes.csv"
    script:
        "../src/features/screen/episodes/screen_episodes.R"

rule resample_episodes:
    input:
        "data/interim/{pid}/{sensor}_episodes.csv"
    params:
        sensor = "{sensor}"
    output:
        "data/interim/{pid}/{sensor}_episodes_resampled.csv"
    script:
        "../src/features/utils/resample_episodes.R"

rule resample_screen_episodes_with_datetime:
    input:
        sensor_input = "data/interim/{pid}/screen_episodes_resampled.csv",
        day_segments = "data/interim/day_segments/{pid}_day_segments.csv"
    params:
        timezones = None,
        fixed_timezone = config["READABLE_DATETIME"]["FIXED_TIMEZONE"],
        day_segments_type = config["DAY_SEGMENTS"]["TYPE"],
        include_past_periodic_segments = config["DAY_SEGMENTS"]["INCLUDE_PAST_PERIODIC_SEGMENTS"]
    output:
        "data/interim/{pid}/screen_episodes_resampled_with_datetime.csv"
    script:
        "../src/data/readable_datetime.R"


rule google_activity_recognition_deltas:
    input:
        expand("data/raw/{{pid}}/{sensor}_with_datetime_unified.csv", sensor=config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["ANDROID"])
    output:
        expand("data/interim/{{pid}}/{sensor}_episodes.csv", sensor=config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["ANDROID"])
    script:
        "../src/features/ar/episodes/activity_recognition_episodes.R"

rule ios_activity_recognition_deltas:
    input:
        expand("data/raw/{{pid}}/{sensor}_with_datetime_unified.csv", sensor=config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["IOS"])
    output:
        expand("data/interim/{{pid}}/{sensor}_episodes.csv", sensor=config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["IOS"])
    script:
        "../src/features/ar/episodes/activity_recognition_episodes.R"

rule locations_python_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor}_processed_{locations_to_use}.csv", sensor=config["LOCATIONS"]["DB_TABLE"], locations_to_use=config["LOCATIONS"]["LOCATIONS_TO_USE"]),
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["LOCATIONS"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
    output:
        "data/interim/{pid}/locations_features/locations_python_{provider_key}.csv"
    script:
        "../src/features/locations/locations_entry.py"

rule locations_r_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor}_processed_{locations_to_use}.csv", sensor=config["LOCATIONS"]["DB_TABLE"], locations_to_use=config["LOCATIONS"]["LOCATIONS_TO_USE"]),
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["LOCATIONS"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}"
    output:
        "data/interim/{pid}/locations_features/locations_r_{provider_key}.csv"
    script:
        "../src/features/locations/locations_entry.R"

rule bluetooth_r_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["BLUETOOTH"]["DB_TABLE"]),
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["BLUETOOTH"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}"
    output:
        "data/interim/{pid}/bluetooth_features/bluetooth_r_{provider_key}.csv"
    script:
        "../src/features/bluetooth/bluetooth_entry.R"

rule bluetooth_python_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["BLUETOOTH"]["DB_TABLE"]),
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["BLUETOOTH"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}"
    output:
        "data/interim/{pid}/bluetooth_features/bluetooth_python_{provider_key}.csv"
    script:
        "../src/features/bluetooth/bluetooth_entry.py"

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
        "data/interim/{pid}/battery_episodes.csv"
    params:
        day_segment = "{day_segment}",
        features = config["BATTERY"]["FEATURES"]
    output:
        "data/processed/{pid}/battery_{day_segment}.csv"
    script:
        "../src/features/battery_features.py"

rule screen_r_features:
    input:
        screen_episodes = "data/interim/{pid}/screen_episodes_resampled_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["SCREEN"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}"
    output:
        "data/interim/{pid}/screen_features/screen_r_{provider_key}.csv"
    script:
        "../src/features/screen/screen_entry.R"

rule screen_python_features:
    input:
        screen_episodes = "data/interim/{pid}/screen_episodes_resampled_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["SCREEN"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}"
    output:
        "data/interim/{pid}/screen_features/screen_python_{provider_key}.csv"
    script:
        "../src/features/screen/screen_entry.py"

rule light_r_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["LIGHT"]["DB_TABLE"]),
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["LIGHT"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}"
    output:
        "data/interim/{pid}/light_features/light_r_{provider_key}.csv"
    script:
        "../src/features/light/light_entry.R"

rule light_python_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["LIGHT"]["DB_TABLE"]),
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["LIGHT"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}"
    output:
        "data/interim/{pid}/light_features/light_python_{provider_key}.csv"
    script:
        "../src/features/light/light_entry.py"

rule conversation_r_features:
    input:
        sensor_data = optional_conversation_input,
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["CONVERSATION"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}"
    output:
        "data/interim/{pid}/conversation_features/conversation_r_{provider_key}.csv"
    script:
        "../src/features/conversation/conversation_entry.R"

rule conversation_python_features:
    input:
        sensor_data = optional_conversation_input,
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["CONVERSATION"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}"
    output:
        "data/interim/{pid}/conversation_features/conversation_python_{provider_key}.csv"
    script:
        "../src/features/conversation/conversation_entry.py"

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

rule applications_foreground_r_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor}_with_datetime_with_genre.csv", sensor=config["APPLICATIONS_FOREGROUND"]["DB_TABLE"]),
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["APPLICATIONS_FOREGROUND"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}"
    output:
        "data/interim/{pid}/applications_foreground_features/applications_foreground_r_{provider_key}.csv"
    script:
        "../src/features/applications_foreground/applications_foreground_entry.R"

rule applications_foreground_python_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor}_with_datetime_with_genre.csv", sensor=config["APPLICATIONS_FOREGROUND"]["DB_TABLE"]),
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["APPLICATIONS_FOREGROUND"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}"
    output:
        "data/interim/{pid}/applications_foreground_features/applications_foreground_python_{provider_key}.csv"
    script:
        "../src/features/applications_foreground/applications_foreground_entry.py"

rule wifi_r_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor_key}_with_datetime_visibleandconnected.csv", sensor_key="WIFI".lower()),
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["WIFI"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}"
    output:
        "data/interim/{pid}/wifi_features/wifi_r_{provider_key}.csv"
    script:
        "../src/features/wifi/wifi_entry.R"

rule wifi_python_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor_key}_with_datetime_visibleandconnected.csv", sensor_key="WIFI".lower()),
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["WIFI"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}"
    output:
        "data/interim/{pid}/wifi_features/wifi_python_{provider_key}.csv"
    script:
        "../src/features/wifi/wifi_entry.py"

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
