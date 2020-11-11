rule join_features_from_providers:
    input:
        location_features = find_features_files
    output:
        "data/processed/features/{pid}/{sensor_key}.csv"
    script:
        "../src/features/join_features_from_providers.R"

rule phone_accelerometer_python_features:
    input:
        sensor_data = "data/raw/{pid}/phone_accelerometer_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_ACCELEROMETER"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_accelerometer"
    output:
        "data/interim/{pid}/phone_accelerometer_features/phone_accelerometer_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule phone_accelerometer_r_features:
    input:
        sensor_data = "data/raw/{pid}/phone_accelerometer_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_ACCELEROMETER"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_accelerometer"
    output:
        "data/interim/{pid}/phone_accelerometer_features/phone_accelerometer_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule activity_recognition_episodes:
    input:
        sensor_data = "data/raw/{pid}/phone_activity_recognition_with_datetime_unified.csv"
    params:
        episode_threshold_between_rows = config["PHONE_BATTERY"]["EPISODE_THRESHOLD_BETWEEN_ROWS"]
    output:
        "data/interim/{pid}/phone_activity_recognition_episodes.csv"
    script:
        "../src/features/phone_activity_recognition/episodes/activity_recognition_episodes.R"

rule phone_activity_recognition_python_features:
    input:
        sensor_episodes = "data/interim/{pid}/phone_activity_recognition_episodes_resampled_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_ACTIVITY_RECOGNITION"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_activity_recognition"
    output:
        "data/interim/{pid}/phone_activity_recognition_features/phone_activity_recognition_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule phone_activity_recognition_r_features:
    input:
        sensor_episodes = "data/interim/{pid}/phone_activity_recognition_episodes_resampled_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_ACTIVITY_RECOGNITION"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_activity_recognition"
    output:
        "data/interim/{pid}/phone_activity_recognition_features/phone_activity_recognition_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule phone_applications_foreground_python_features:
    input:
        sensor_data = "data/raw/{pid}/phone_applications_foreground_with_datetime_with_categories.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_APPLICATIONS_FOREGROUND"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_applications_foreground"
    output:
        "data/interim/{pid}/phone_applications_foreground_features/phone_applications_foreground_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule phone_applications_foreground_r_features:
    input:
        sensor_data = "data/raw/{pid}/phone_applications_foreground_with_datetime_with_categories.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_APPLICATIONS_FOREGROUND"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_applications_foreground"
    output:
        "data/interim/{pid}/phone_applications_foreground_features/phone_applications_foreground_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule battery_episodes:
    input:
        "data/raw/{pid}/phone_battery_raw.csv"
    params:
        episode_threshold_between_rows = config["PHONE_BATTERY"]["EPISODE_THRESHOLD_BETWEEN_ROWS"]
    output:
        "data/interim/{pid}/phone_battery_episodes.csv"
    script:
        "../src/features/phone_battery/episodes/battery_episodes.R"

rule phone_battery_python_features:
    input:
        sensor_episodes = "data/interim/{pid}/phone_battery_episodes_resampled_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_BATTERY"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_battery"
    output:
        "data/interim/{pid}/phone_battery_features/phone_battery_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule phone_battery_r_features:
    input:
        sensor_episodes = "data/interim/{pid}/phone_battery_episodes_resampled_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_BATTERY"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_battery"
    output:
        "data/interim/{pid}/phone_battery_features/phone_battery_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule phone_bluetooth_python_features:
    input:
        sensor_data = "data/raw/{pid}/phone_bluetooth_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_BLUETOOTH"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_bluetooth"
    output:
        "data/interim/{pid}/phone_bluetooth_features/phone_bluetooth_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule phone_bluetooth_r_features:
    input:
        sensor_data = "data/raw/{pid}/phone_bluetooth_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_BLUETOOTH"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_bluetooth"
    output:
        "data/interim/{pid}/phone_bluetooth_features/phone_bluetooth_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule calls_python_features:
    input:
        sensor_data = "data/raw/{pid}/phone_calls_with_datetime_unified.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_CALLS"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_calls"
    output:
        "data/interim/{pid}/phone_calls_features/phone_calls_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule calls_r_features:
    input:
        sensor_data = "data/raw/{pid}/phone_calls_with_datetime_unified.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_CALLS"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_calls"
    output:
        "data/interim/{pid}/phone_calls_features/phone_calls_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule conversation_python_features:
    input:
        sensor_data = "data/raw/{pid}/phone_conversation_with_datetime_unified.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_CONVERSATION"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_conversation"
    output:
        "data/interim/{pid}/phone_conversation_features/phone_conversation_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule conversation_r_features:
    input:
        sensor_data = "data/raw/{pid}/phone_conversation_with_datetime_unified.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_CONVERSATION"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_conversation"
    output:
        "data/interim/{pid}/phone_conversation_features/phone_conversation_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule phone_light_python_features:
    input:
        sensor_data = "data/raw/{pid}/phone_light_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_LIGHT"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_light"
    output:
        "data/interim/{pid}/phone_light_features/phone_light_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule phone_light_r_features:
    input:
        sensor_data = "data/raw/{pid}/phone_light_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_LIGHT"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_light"
    output:
        "data/interim/{pid}/phone_light_features/phone_light_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule phone_locations_python_features:
    input:
        sensor_data = "data/interim/{pid}/phone_locations_processed_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_LOCATIONS"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_locations"
    output:
        "data/interim/{pid}/phone_locations_features/phone_locations_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule phone_locations_r_features:
    input:
        sensor_data = "data/interim/{pid}/phone_locations_processed_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_LOCATIONS"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_locations"
    output:
        "data/interim/{pid}/phone_locations_features/phone_locations_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule phone_messages_python_features:
    input:
        sensor_data = "data/raw/{pid}/phone_messages_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_MESSAGES"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_messages"
    output:
        "data/interim/{pid}/phone_messages_features/phone_messages_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule phone_messages_r_features:
    input:
        sensor_data = "data/raw/{pid}/phone_messages_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_MESSAGES"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_messages"
    output:
        "data/interim/{pid}/phone_messages_features/phone_messages_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule screen_episodes:
    input:
        screen = "data/raw/{pid}/phone_screen_with_datetime_unified.csv"
    output:
        "data/interim/{pid}/phone_screen_episodes.csv"
    script:
        "../src/features/phone_screen/episodes/screen_episodes.R"

rule phone_screen_python_features:
    input:
        sensor_episodes = "data/interim/{pid}/phone_screen_episodes_resampled_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_SCREEN"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_screen"
    output:
        "data/interim/{pid}/phone_screen_features/phone_screen_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule phone_screen_r_features:
    input:
        sensor_episodes = "data/interim/{pid}/phone_screen_episodes_resampled_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_SCREEN"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_screen"
    output:
        "data/interim/{pid}/phone_screen_features/phone_screen_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule phone_wifi_connected_python_features:
    input:
        sensor_data = "data/raw/{pid}/phone_wifi_connected_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_WIFI_CONNECTED"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_wifi_connected"
    output:
        "data/interim/{pid}/phone_wifi_connected_features/phone_wifi_connected_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule phone_wifi_connected_r_features:
    input:
        sensor_data = "data/raw/{pid}/phone_wifi_connected_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_WIFI_CONNECTED"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_wifi_connected"
    output:
        "data/interim/{pid}/phone_wifi_connected_features/phone_wifi_connected_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule phone_wifi_visible_python_features:
    input:
        sensor_data = "data/raw/{pid}/phone_wifi_visible_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_WIFI_VISIBLE"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_wifi_visible"
    output:
        "data/interim/{pid}/phone_wifi_visible_features/phone_wifi_visible_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule phone_wifi_visible_r_features:
    input:
        sensor_data = "data/raw/{pid}/phone_wifi_visible_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_WIFI_VISIBLE"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_wifi_visible"
    output:
        "data/interim/{pid}/phone_wifi_visible_features/phone_wifi_visible_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule fitbit_heartrate_summary_python_features:
    input:
        sensor_data = "data/raw/{pid}/fitbit_heartrate_summary_parsed_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["FITBIT_HEARTRATE_SUMMARY"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "fitbit_heartrate_summary"
    output:
        "data/interim/{pid}/fitbit_heartrate_summary_features/fitbit_heartrate_summary_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule fitbit_heartrate_summary_r_features:
    input:
        sensor_data = "data/raw/{pid}/fitbit_heartrate_summary_parsed_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["FITBIT_HEARTRATE_SUMMARY"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "fitbit_heartrate_summary"
    output:
        "data/interim/{pid}/fitbit_heartrate_summary_features/fitbit_heartrate_summary_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule fitbit_heartrate_intraday_python_features:
    input:
        sensor_data = "data/raw/{pid}/fitbit_heartrate_intraday_parsed_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["FITBIT_HEARTRATE_INTRADAY"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "fitbit_heartrate_intraday"
    output:
        "data/interim/{pid}/fitbit_heartrate_intraday_features/fitbit_heartrate_intraday_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule fitbit_heartrate_intraday_r_features:
    input:
        sensor_data = "data/raw/{pid}/fitbit_heartrate_intraday_parsed_with_datetime.csv",
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["FITBIT_HEARTRATE_INTRADAY"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "fitbit_heartrate_intraday"
    output:
        "data/interim/{pid}/fitbit_heartrate_intraday_features/fitbit_heartrate_intraday_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule fitbit_steps_python_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/fitbit_steps_{fitbit_data_type}_parsed_with_datetime.csv", fitbit_data_type=["summary", "intraday"]),
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["FITBIT_STEPS"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "fitbit_steps"
    output:
        "data/interim/{pid}/fitbit_steps_features/fitbit_steps_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule fitbit_steps_r_features:
    input:
        sensor_data = expand("data/raw/{{pid}}/fitbit_steps_{fitbit_data_type}_parsed_with_datetime.csv", fitbit_data_type=["summary", "intraday"]),
        day_segments_labels = "data/interim/day_segments/{pid}_day_segments_labels.csv"
    params:
        provider = lambda wildcards: config["FITBIT_STEPS"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "fitbit_steps"
    output:
        "data/interim/{pid}/fitbit_steps_features/fitbit_steps_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

# rule fitbit_step_features:
#     input:
#         step_data = "data/raw/{pid}/fitbit_step_intraday_with_datetime.csv",
#         sleep_data = optional_steps_sleep_input
#     params:
#         day_segment = "{day_segment}",
#         features_all_steps = config["STEP"]["FEATURES"]["ALL_STEPS"],
#         features_sedentary_bout = config["STEP"]["FEATURES"]["SEDENTARY_BOUT"],
#         features_active_bout = config["STEP"]["FEATURES"]["ACTIVE_BOUT"],
#         threshold_active_bout = config["STEP"]["THRESHOLD_ACTIVE_BOUT"],
#         include_zero_step_rows = config["STEP"]["INCLUDE_ZERO_STEP_ROWS"],
#         exclude_sleep = config["STEP"]["EXCLUDE_SLEEP"]["EXCLUDE"],
#         exclude_sleep_type = config["STEP"]["EXCLUDE_SLEEP"]["TYPE"],
#         exclude_sleep_fixed_start = config["STEP"]["EXCLUDE_SLEEP"]["FIXED"]["START"],
#         exclude_sleep_fixed_end = config["STEP"]["EXCLUDE_SLEEP"]["FIXED"]["END"],
#     output:
#         "data/processed/{pid}/fitbit_step_{day_segment}.csv"
#     script:
#         "../src/features/fitbit_step_features.py"

# rule fitbit_sleep_features:
#     input:
#         sleep_summary_data = "data/raw/{pid}/fitbit_sleep_summary_with_datetime.csv",
#         sleep_intraday_data = "data/raw/{pid}/fitbit_sleep_intraday_with_datetime.csv"
#     params:
#         day_segment = "{day_segment}",
#         summary_features = config["SLEEP"]["SUMMARY_FEATURES"],
#         sleep_types = config["SLEEP"]["SLEEP_TYPES"]
#     output:
#         "data/processed/{pid}/fitbit_sleep_{day_segment}.csv"
#     script:
#         "../src/features/fitbit_sleep_features.py"
