rule join_features_from_providers:
    input:
        sensor_features = find_features_files
    wildcard_constraints:
        sensor_key = '(phone|fitbit|empatica).*'
    output:
        "data/processed/features/{pid}/{sensor_key}.csv"
    script:
        "../src/features/utils/join_features_from_providers.R"

rule phone_data_yield_python_features:
    input:
        sensor_data = "data/interim/{pid}/phone_yielded_timestamps_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_DATA_YIELD"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_data_yield"
    output:
        "data/interim/{pid}/phone_data_yield_features/phone_data_yield_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule phone_data_yield_r_features:
    input:
        sensor_data = "data/interim/{pid}/phone_yielded_timestamps_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_DATA_YIELD"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_data_yield"
    output:
        "data/interim/{pid}/phone_data_yield_features/phone_data_yield_r_{provider_key}.csv"
    script:
        "../src/features/entry.R" 

rule phone_accelerometer_python_features:
    input:
        sensor_data = "data/raw/{pid}/phone_accelerometer_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
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
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
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
        sensor_data = "data/raw/{pid}/phone_activity_recognition_with_datetime.csv"
    params:
        episode_threshold_between_rows = config["PHONE_ACTIVITY_RECOGNITION"]["EPISODE_THRESHOLD_BETWEEN_ROWS"]
    output:
        "data/interim/{pid}/phone_activity_recognition_episodes.csv"
    script:
        "../src/features/phone_activity_recognition/episodes/activity_recognition_episodes.R"

rule phone_activity_recognition_python_features:
    input:
        sensor_episodes = "data/interim/{pid}/phone_activity_recognition_episodes_resampled_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
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
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_ACTIVITY_RECOGNITION"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_activity_recognition"
    output:
        "data/interim/{pid}/phone_activity_recognition_features/phone_activity_recognition_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule phone_applications_crashes_python_features:
    input:
        sensor_data = "data/raw/{pid}/phone_applications_crashes_with_datetime_with_categories.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_APPLICATIONS_CRASHES"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_applications_crashes"
    output:
        "data/interim/{pid}/phone_applications_crashes_features/phone_applications_crashes_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule phone_applications_crashes_r_features:
    input:
        sensor_data = "data/raw/{pid}/phone_applications_crashes_with_datetime_with_categories.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_APPLICATIONS_CRASHES"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_applications_crashes"
    output:
        "data/interim/{pid}/phone_applications_crashes_features/phone_applications_crashes_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule phone_applications_foreground_python_features:
    input:
        sensor_data = "data/raw/{pid}/phone_applications_foreground_with_datetime_with_categories.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
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
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_APPLICATIONS_FOREGROUND"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_applications_foreground"
    output:
        "data/interim/{pid}/phone_applications_foreground_features/phone_applications_foreground_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule phone_applications_notifications_python_features:
    input:
        sensor_data = "data/raw/{pid}/phone_applications_notifications_with_datetime_with_categories.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_APPLICATIONS_NOTIFICATIONS"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_applications_notifications"
    output:
        "data/interim/{pid}/phone_applications_notifications_features/phone_applications_notifications_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule phone_applications_notifications_r_features:
    input:
        sensor_data = "data/raw/{pid}/phone_applications_notifications_with_datetime_with_categories.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_APPLICATIONS_NOTIFICATIONS"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_applications_notifications"
    output:
        "data/interim/{pid}/phone_applications_notifications_features/phone_applications_notifications_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule phone_log_python_features:
    input:
        sensor_data = "data/raw/{pid}/phone_log_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_LOG"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_log"
    output:
        "data/interim/{pid}/phone_log_features/phone_log_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule phone_log_r_features:
    input:
        sensor_data = "data/raw/{pid}/phone_log_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_LOG"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_log"
    output:
        "data/interim/{pid}/phone_log_features/phone_log_r_{provider_key}.csv"
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
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
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
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
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
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
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
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
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
        sensor_data = "data/raw/{pid}/phone_calls_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
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
        sensor_data = "data/raw/{pid}/phone_calls_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
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
        sensor_data = "data/raw/{pid}/phone_conversation_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
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
        sensor_data = "data/raw/{pid}/phone_conversation_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_CONVERSATION"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_conversation"
    output:
        "data/interim/{pid}/phone_conversation_features/phone_conversation_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule phone_keyboard_python_features:
    input:
        sensor_data = "data/raw/{pid}/phone_keyboard_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_KEYBOARD"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_keyboard"
    output:
        "data/interim/{pid}/phone_keyboard_features/phone_keyboard_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule phone_keyboard_r_features:
    input:
        sensor_data = "data/raw/{pid}/phone_keyboard_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_KEYBOARD"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_keyboard"
    output:
        "data/interim/{pid}/phone_keyboard_features/phone_keyboard_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule phone_light_python_features:
    input:
        sensor_data = "data/raw/{pid}/phone_light_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
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
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
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
        sensor_data = "data/interim/{pid}/phone_locations_processed_with_datetime_with_home.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_LOCATIONS"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_locations"
    output:
        "data/interim/{pid}/phone_locations_features/phone_locations_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule phone_locations_barnett_daily_features:
    input:
        sensor_data = "data/interim/{pid}/phone_locations_processed_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv",
    params:
        provider = lambda wildcards: config["PHONE_LOCATIONS"]["PROVIDERS"]["BARNETT"],
    output:
        "data/interim/{pid}/phone_locations_barnett_daily.csv"
    script:
        "../src/features/phone_locations/barnett/daily_features.R"

rule phone_locations_r_features:
    input:
        sensor_data = "data/interim/{pid}/phone_locations_processed_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv",
        barnett_daily =  get_barnett_daily
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
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
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
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
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
        screen = "data/raw/{pid}/phone_screen_with_datetime.csv"
    output:
        "data/interim/{pid}/phone_screen_episodes.csv"
    script:
        "../src/features/phone_screen/episodes/screen_episodes.R"

rule phone_screen_python_features:
    input:
        sensor_episodes = "data/interim/{pid}/phone_screen_episodes_resampled_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
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
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
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
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
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
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
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
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
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
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["PHONE_WIFI_VISIBLE"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "phone_wifi_visible"
    output:
        "data/interim/{pid}/phone_wifi_visible_features/phone_wifi_visible_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule fitbit_data_yield_python_features:
    input:
        sensor_data = "data/raw/{pid}/fitbit_heartrate_intraday_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["FITBIT_DATA_YIELD"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "fitbit_data_yield"
    output:
        "data/interim/{pid}/fitbit_data_yield_features/fitbit_data_yield_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule fitbit_data_yield_r_features:
    input:
        sensor_data = "data/raw/{pid}/fitbit_heartrate_intraday_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["FITBIT_DATA_YIELD"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "fitbit_data_yield"
    output:
        "data/interim/{pid}/fitbit_data_yield_features/fitbit_data_yield_r_{provider_key}.csv"
    script:
        "../src/features/entry.R" 

rule fitbit_heartrate_summary_python_features:
    input:
        sensor_data = "data/raw/{pid}/fitbit_heartrate_summary_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
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
        sensor_data = "data/raw/{pid}/fitbit_heartrate_summary_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
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
        sensor_data = "data/raw/{pid}/fitbit_heartrate_intraday_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
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
        sensor_data = "data/raw/{pid}/fitbit_heartrate_intraday_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["FITBIT_HEARTRATE_INTRADAY"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "fitbit_heartrate_intraday"
    output:
        "data/interim/{pid}/fitbit_heartrate_intraday_features/fitbit_heartrate_intraday_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule fitbit_steps_summary_python_features:
    input:
        sensor_data = "data/raw/{pid}/fitbit_steps_summary_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["FITBIT_STEPS_SUMMARY"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "fitbit_steps_summary"
    output:
        "data/interim/{pid}/fitbit_steps_summary_features/fitbit_steps_summary_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule fitbit_steps_summary_r_features:
    input:
        sensor_data = "data/raw/{pid}/fitbit_steps_summary_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["FITBIT_STEPS_SUMMARY"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "fitbit_steps_summary"
    output:
        "data/interim/{pid}/fitbit_steps_summary_features/fitbit_steps_summary_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule fitbit_steps_intraday_python_features:
    input:
        sensor_data = "data/raw/{pid}/fitbit_steps_intraday_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["FITBIT_STEPS_INTRADAY"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "fitbit_steps_intraday"
    output:
        "data/interim/{pid}/fitbit_steps_intraday_features/fitbit_steps_intraday_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule fitbit_steps_intraday_r_features:
    input:
        sensor_data = "data/raw/{pid}/fitbit_steps_intraday_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["FITBIT_STEPS_INTRADAY"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "fitbit_steps_intraday"
    output:
        "data/interim/{pid}/fitbit_steps_intraday_features/fitbit_steps_intraday_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule fitbit_sleep_summary_python_features:
    input:
        sensor_data = "data/raw/{pid}/fitbit_sleep_summary_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["FITBIT_SLEEP_SUMMARY"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "fitbit_sleep_summary"
    output:
        "data/interim/{pid}/fitbit_sleep_summary_features/fitbit_sleep_summary_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule fitbit_sleep_summary_r_features:
    input:
        sensor_data = "data/raw/{pid}/fitbit_sleep_summary_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["FITBIT_SLEEP_SUMMARY"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "fitbit_sleep_summary"
    output:
        "data/interim/{pid}/fitbit_sleep_summary_features/fitbit_sleep_summary_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule sleep_intraday_episodes:
    input:
        sleep_intraday = "data/raw/{pid}/fitbit_sleep_intraday_with_datetime.csv"
    output:
        "data/interim/{pid}/fitbit_sleep_intraday_episodes.csv"
    script:
        "../src/features/fitbit_sleep_intraday/episodes/sleep_intraday_episodes.py"

rule fitbit_sleep_intraday_python_features:
    input:
        sensor_data = "data/interim/{pid}/fitbit_sleep_intraday_episodes_resampled_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["FITBIT_SLEEP_INTRADAY"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "fitbit_sleep_intraday"
    output:
        "data/interim/{pid}/fitbit_sleep_intraday_features/fitbit_sleep_intraday_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule fitbit_sleep_intraday_r_features:
    input:
        sensor_data = "data/interim/{pid}/fitbit_sleep_intraday_episodes_resampled_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["FITBIT_SLEEP_INTRADAY"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "fitbit_sleep_intraday"
    output:
        "data/interim/{pid}/fitbit_sleep_intraday_features/fitbit_sleep_intraday_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule merge_sensor_features_for_individual_participants:
    input:
        feature_files = input_merge_sensor_features_for_individual_participants
    output:
        "data/processed/features/{pid}/all_sensor_features.csv"
    script:
        "../src/features/utils/merge_sensor_features_for_individual_participants.R"

rule merge_sensor_features_for_all_participants:
    input:
        feature_files = expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"])
    output:
        "data/processed/features/all_participants/all_sensor_features.csv"
    script:
        "../src/features/utils/merge_sensor_features_for_all_participants.R"

rule empatica_accelerometer_python_features:
    input:
        sensor_data = "data/raw/{pid}/empatica_accelerometer_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["EMPATICA_ACCELEROMETER"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "empatica_accelerometer"
    output:
        "data/interim/{pid}/empatica_accelerometer_features/empatica_accelerometer_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule empatica_accelerometer_r_features:
    input:
        sensor_data = "data/raw/{pid}/empatica_accelerometer_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["EMPATICA_ACCELEROMETER"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "empatica_accelerometer"
    output:
        "data/interim/{pid}/empatica_accelerometer_features/empatica_accelerometer_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule empatica_heartrate_python_features:
    input:
        sensor_data = "data/raw/{pid}/empatica_heartrate_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["EMPATICA_HEARTRATE"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "empatica_heartrate"
    output:
        "data/interim/{pid}/empatica_heartrate_features/empatica_heartrate_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule empatica_heartrate_r_features:
    input:
        sensor_data = "data/raw/{pid}/empatica_heartrate_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["EMPATICA_HEARTRATE"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "empatica_heartrate"
    output:
        "data/interim/{pid}/empatica_heartrate_features/empatica_heartrate_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule empatica_temperature_python_features:
    input:
        sensor_data = "data/raw/{pid}/empatica_temperature_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["EMPATICA_TEMPERATURE"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "empatica_temperature"
    output:
        "data/interim/{pid}/empatica_temperature_features/empatica_temperature_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule empatica_temperature_r_features:
    input:
        sensor_data = "data/raw/{pid}/empatica_temperature_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["EMPATICA_TEMPERATURE"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "empatica_temperature"
    output:
        "data/interim/{pid}/empatica_temperature_features/empatica_temperature_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule empatica_electrodermal_activity_python_features:
    input:
        sensor_data = "data/raw/{pid}/empatica_electrodermal_activity_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["EMPATICA_ELECTRODERMAL_ACTIVITY"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "empatica_electrodermal_activity"
    output:
        "data/interim/{pid}/empatica_electrodermal_activity_features/empatica_electrodermal_activity_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule empatica_electrodermal_activity_r_features:
    input:
        sensor_data = "data/raw/{pid}/empatica_electrodermal_activity_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["EMPATICA_ELECTRODERMAL_ACTIVITY"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "empatica_electrodermal_activity"
    output:
        "data/interim/{pid}/empatica_electrodermal_activity_features/empatica_electrodermal_activity_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule empatica_blood_volume_pulse_python_features:
    input:
        sensor_data = "data/raw/{pid}/empatica_blood_volume_pulse_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["EMPATICA_BLOOD_VOLUME_PULSE"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "empatica_blood_volume_pulse"
    output:
        "data/interim/{pid}/empatica_blood_volume_pulse_features/empatica_blood_volume_pulse_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule empatica_blood_volume_pulse_r_features:
    input:
        sensor_data = "data/raw/{pid}/empatica_blood_volume_pulse_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["EMPATICA_BLOOD_VOLUME_PULSE"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "empatica_blood_volume_pulse"
    output:
        "data/interim/{pid}/empatica_blood_volume_pulse_features/empatica_blood_volume_pulse_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule empatica_inter_beat_interval_python_features:
    input:
        sensor_data = "data/raw/{pid}/empatica_inter_beat_interval_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["EMPATICA_INTER_BEAT_INTERVAL"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "empatica_inter_beat_interval"
    output:
        "data/interim/{pid}/empatica_inter_beat_interval_features/empatica_inter_beat_interval_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule empatica_inter_beat_interval_r_features:
    input:
        sensor_data = "data/raw/{pid}/empatica_inter_beat_interval_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["EMPATICA_INTER_BEAT_INTERVAL"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "empatica_inter_beat_interval"
    output:
        "data/interim/{pid}/empatica_inter_beat_interval_features/empatica_inter_beat_interval_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"

rule empatica_tags_python_features:
    input:
        sensor_data = "data/raw/{pid}/empatica_tags_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["EMPATICA_TAGS"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "empatica_tags"
    output:
        "data/interim/{pid}/empatica_tags_features/empatica_tags_python_{provider_key}.csv"
    script:
        "../src/features/entry.py"

rule empatica_tags_r_features:
    input:
        sensor_data = "data/raw/{pid}/empatica_tags_with_datetime.csv",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        provider = lambda wildcards: config["EMPATICA_TAGS"]["PROVIDERS"][wildcards.provider_key.upper()],
        provider_key = "{provider_key}",
        sensor_key = "empatica_tags"
    output:
        "data/interim/{pid}/empatica_tags_features/empatica_tags_r_{provider_key}.csv"
    script:
        "../src/features/entry.R"
