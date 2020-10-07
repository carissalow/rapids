rule join_features_from_providers:
    input:
        location_features = find_features_files
    output:
        "data/processed/features/{pid}/{sensor_key}.csv"
    script:
        "../src/features/join_features_from_providers.R"

rule resample_episodes:
    input:
        "data/interim/{pid}/{sensor}_episodes.csv"
    output:
        "data/interim/{pid}/{sensor}_episodes_resampled.csv"
    script:
        "../src/features/utils/resample_episodes.R"

rule resample_episodes_with_datetime:
    input:
        sensor_input = "data/interim/{pid}/{sensor}_episodes_resampled.csv",
        day_segments = "data/interim/day_segments/{pid}_day_segments.csv"
    params:
        timezones = None,
        fixed_timezone = config["READABLE_DATETIME"]["FIXED_TIMEZONE"],
        day_segments_type = config["DAY_SEGMENTS"]["TYPE"],
        include_past_periodic_segments = config["DAY_SEGMENTS"]["INCLUDE_PAST_PERIODIC_SEGMENTS"]
    output:
        "data/interim/{pid}/{sensor}_episodes_resampled_with_datetime.csv"
    script:
        "../src/data/readable_datetime.R"

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

rule activity_recognition_episodes:
    input:
        optional_ar_input
    output:
        "data/interim/{pid}/activity_recognition_episodes.csv"
    script:
        "../src/features/activity_recognition/episodes/activity_recognition_episodes.R"

rule activity_recognition_r_features:
    input:
        sensor_episodes = "data/interim/{pid}/activity_recognition_episodes_resampled_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["ACTIVITY_RECOGNITION"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
        sensor_key = "activity_recognition"
    output:
        "data/interim/{pid}/activity_recognition_features/activity_recognition_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule activity_recognition_python_features:
    input:
        sensor_episodes = "data/interim/{pid}/activity_recognition_episodes_resampled_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["ACTIVITY_RECOGNITION"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
        sensor_key = "activity_recognition"
    output:
        "data/interim/{pid}/activity_recognition_features/activity_recognition_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule applications_foreground_r_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor}_with_datetime_with_genre.csv", sensor=config["APPLICATIONS_FOREGROUND"]["DB_TABLE"])[0],
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["APPLICATIONS_FOREGROUND"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
        sensor_key = "applications_foreground"
    output:
        "data/interim/{pid}/applications_foreground_features/applications_foreground_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule applications_foreground_python_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor}_with_datetime_with_genre.csv", sensor=config["APPLICATIONS_FOREGROUND"]["DB_TABLE"])[0],
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["APPLICATIONS_FOREGROUND"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
        sensor_key = "applications_foreground"
    output:
        "data/interim/{pid}/applications_foreground_features/applications_foreground_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule battery_episodes:
    input:
        expand("data/raw/{{pid}}/{sensor}_raw.csv", sensor=config["BATTERY"]["DB_TABLE"])
    output:
        "data/interim/{pid}/battery_episodes.csv"
    script:
        "../src/features/battery/episodes/battery_episodes.R"

rule battery_r_features:
    input:
        sensor_episodes = "data/interim/{pid}/battery_episodes_resampled_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["BATTERY"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
        sensor_key = "battery"
    output:
        "data/interim/{pid}/battery_features/battery_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule battery_python_features:
    input:
        sensor_episodes = "data/interim/{pid}/battery_episodes_resampled_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["BATTERY"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
        sensor_key = "battery"
    output:
        "data/interim/{pid}/battery_features/battery_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule bluetooth_r_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["BLUETOOTH"]["DB_TABLE"])[0],
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["BLUETOOTH"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
        sensor_key = "bluetooth"
    output:
        "data/interim/{pid}/bluetooth_features/bluetooth_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule bluetooth_python_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["BLUETOOTH"]["DB_TABLE"])[0],
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["BLUETOOTH"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
        sensor_key = "bluetooth"
    output:
        "data/interim/{pid}/bluetooth_features/bluetooth_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule calls_r_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor}_with_datetime_unified.csv", sensor=config["CALLS"]["DB_TABLE"])[0],
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["CALLS"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
        sensor_key = "calls"
    output:
        "data/interim/{pid}/calls_features/calls_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule calls_python_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor}_with_datetime_unified.csv", sensor=config["CALLS"]["DB_TABLE"])[0],
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["CALLS"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
        sensor_key = "calls"
    output:
        "data/interim/{pid}/calls_features/calls_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule conversation_r_features:
    input:
        sensor_data = optional_conversation_input,
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["CONVERSATION"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
        sensor_key = "conversation"
    output:
        "data/interim/{pid}/conversation_features/conversation_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule conversation_python_features:
    input:
        sensor_data = optional_conversation_input,
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["CONVERSATION"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
        sensor_key = "conversation"
    output:
        "data/interim/{pid}/conversation_features/conversation_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule light_r_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["LIGHT"]["DB_TABLE"])[0],
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["LIGHT"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
        sensor_key = "light"
    output:
        "data/interim/{pid}/light_features/light_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule light_python_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["LIGHT"]["DB_TABLE"])[0],
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["LIGHT"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
        sensor_key = "light"
    output:
        "data/interim/{pid}/light_features/light_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule locations_r_features:
    input:
        sensor_data = expand("data/interim/{{pid}}/{sensor}_processed_{locations_to_use}_with_datetime.csv", sensor=config["LOCATIONS"]["DB_TABLE"], locations_to_use=config["LOCATIONS"]["LOCATIONS_TO_USE"])[0],
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["LOCATIONS"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
        sensor_key = "locations"
    output:
        "data/interim/{pid}/locations_features/locations_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule locations_python_features:
    input:
        sensor_data = expand("data/interim/{{pid}}/{sensor}_processed_{locations_to_use}_with_datetime.csv", sensor=config["LOCATIONS"]["DB_TABLE"], locations_to_use=config["LOCATIONS"]["LOCATIONS_TO_USE"])[0],
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["LOCATIONS"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
        sensor_key = "locations"
    output:
        "data/interim/{pid}/locations_features/locations_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule messages_r_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["MESSAGES"]["DB_TABLE"])[0],
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["MESSAGES"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
        sensor_key = "messages"
    output:
        "data/interim/{pid}/messages_features/messages_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule messages_python_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["MESSAGES"]["DB_TABLE"])[0],
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["MESSAGES"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
        sensor_key = "messages"
    output:
        "data/interim/{pid}/messages_features/messages_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule screen_episodes:
    input:
        screen = expand("data/raw/{{pid}}/{sensor}_with_datetime_unified.csv", sensor=config["SCREEN"]["DB_TABLE"])
    output:
        "data/interim/{pid}/screen_episodes.csv"
    script:
        "../src/features/screen/episodes/screen_episodes.R"

rule screen_r_features:
    input:
        sensor_episodes = "data/interim/{pid}/screen_episodes_resampled_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["SCREEN"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
        sensor_key = "screen"
    output:
        "data/interim/{pid}/screen_features/screen_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule screen_python_features:
    input:
        sensor_episodes = "data/interim/{pid}/screen_episodes_resampled_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["SCREEN"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
        sensor_key = "screen"
    output:
        "data/interim/{pid}/screen_features/screen_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule wifi_r_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor_key}_with_datetime_visibleandconnected.csv", sensor_key="WIFI".lower())[0],
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["WIFI"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
        sensor_key = "wifi"
    output:
        "data/interim/{pid}/wifi_features/wifi_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule wifi_python_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/{sensor_key}_with_datetime_visibleandconnected.csv", sensor_key="WIFI".lower())[0],
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["WIFI"]["PROVIDERS"][wildcards.provider_key],
        provider_key = "{provider_key}",
        sensor_key = "wifi"
    output:
        "data/interim/{pid}/wifi_features/wifi_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

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
