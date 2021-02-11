from snakemake.utils import validate
configfile: "config.yaml"
validate(config, "tools/config.schema.yaml")
include: "rules/common.smk"
include: "rules/renv.smk"
include: "rules/preprocessing.smk"
include: "rules/features.smk"
include: "rules/reports.smk"

import itertools

files_to_compute = []

if len(config["PIDS"]) == 0:
    raise ValueError("Add participants IDs to PIDS in config.yaml. Remember to create their participant files in data/external")

for provider in config["PHONE_DATA_YIELD"]["PROVIDERS"].keys():
    if config["PHONE_DATA_YIELD"]["PROVIDERS"][provider]["COMPUTE"]:
        
        allowed_phone_sensors = get_phone_sensor_names()
        if not (set(config["PHONE_DATA_YIELD"]["SENSORS"]) <= set(allowed_phone_sensors)):
            raise ValueError('\nInvalid sensor(s) for PHONE_DATA_YIELD. config["PHONE_DATA_YIELD"]["SENSORS"] can have '
                            'one or more of the following phone sensors: {}.\nInstead you provided "{}".\n'
                            'Keep in mind that the sensors\' TABLE attribute must point to a valid database table'\
                            .format(', '.join(allowed_phone_sensors),
                                    ', '.join(set(config["PHONE_DATA_YIELD"]["SENSORS"]) - set(allowed_phone_sensors))))
        
        files_to_compute.extend(expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=map(str.lower, config["PHONE_DATA_YIELD"]["SENSORS"])))
        files_to_compute.extend(expand("data/interim/{pid}/phone_yielded_timestamps.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_yielded_timestamps_with_datetime.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_data_yield_features/phone_data_yield_{language}_{provider_key}.csv", pid=config["PIDS"], language=config["PHONE_DATA_YIELD"]["PROVIDERS"][provider]["SRC_LANGUAGE"].lower(), provider_key=provider.lower()))
        files_to_compute.extend(expand("data/processed/features/{pid}/phone_data_yield.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
        files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

for provider in config["PHONE_MESSAGES"]["PROVIDERS"].keys():
    if config["PHONE_MESSAGES"]["PROVIDERS"][provider]["COMPUTE"]:
        files_to_compute.extend(expand("data/raw/{pid}/phone_messages_raw.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/phone_messages_with_datetime.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_messages_features/phone_messages_{language}_{provider_key}.csv", pid=config["PIDS"], language=config["PHONE_MESSAGES"]["PROVIDERS"][provider]["SRC_LANGUAGE"].lower(), provider_key=provider.lower()))
        files_to_compute.extend(expand("data/processed/features/{pid}/phone_messages.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
        files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

for provider in config["PHONE_CALLS"]["PROVIDERS"].keys():
    if config["PHONE_CALLS"]["PROVIDERS"][provider]["COMPUTE"]:
        files_to_compute.extend(expand("data/raw/{pid}/phone_calls_raw.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/phone_calls_with_datetime.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/phone_calls_with_datetime_unified.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_calls_features/phone_calls_{language}_{provider_key}.csv", pid=config["PIDS"], language=config["PHONE_CALLS"]["PROVIDERS"][provider]["SRC_LANGUAGE"].lower(), provider_key=provider.lower()))
        files_to_compute.extend(expand("data/processed/features/{pid}/phone_calls.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
        files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

for provider in config["PHONE_BLUETOOTH"]["PROVIDERS"].keys():
    if config["PHONE_BLUETOOTH"]["PROVIDERS"][provider]["COMPUTE"]:
        files_to_compute.extend(expand("data/raw/{pid}/phone_bluetooth_raw.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/phone_bluetooth_with_datetime.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_bluetooth_features/phone_bluetooth_{language}_{provider_key}.csv", pid=config["PIDS"], language=config["PHONE_BLUETOOTH"]["PROVIDERS"][provider]["SRC_LANGUAGE"].lower(), provider_key=provider.lower()))
        files_to_compute.extend(expand("data/processed/features/{pid}/phone_bluetooth.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
        files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

for provider in config["PHONE_ACTIVITY_RECOGNITION"]["PROVIDERS"].keys():
    if config["PHONE_ACTIVITY_RECOGNITION"]["PROVIDERS"][provider]["COMPUTE"]:
        files_to_compute.extend(expand("data/raw/{pid}/phone_activity_recognition_raw.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/phone_activity_recognition_with_datetime.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/phone_activity_recognition_with_datetime_unified.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_activity_recognition_episodes.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_activity_recognition_episodes_resampled.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_activity_recognition_episodes_resampled_with_datetime.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_activity_recognition_features/phone_activity_recognition_{language}_{provider_key}.csv", pid=config["PIDS"], language=config["PHONE_ACTIVITY_RECOGNITION"]["PROVIDERS"][provider]["SRC_LANGUAGE"].lower(), provider_key=provider.lower()))
        files_to_compute.extend(expand("data/processed/features/{pid}/phone_activity_recognition.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
        files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

for provider in config["PHONE_BATTERY"]["PROVIDERS"].keys():
    if config["PHONE_BATTERY"]["PROVIDERS"][provider]["COMPUTE"]:
        files_to_compute.extend(expand("data/raw/{pid}/phone_battery_raw.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_battery_episodes.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_battery_episodes_resampled.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_battery_episodes_resampled_with_datetime.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_battery_features/phone_battery_{language}_{provider_key}.csv", pid=config["PIDS"], language=config["PHONE_BATTERY"]["PROVIDERS"][provider]["SRC_LANGUAGE"].lower(), provider_key=provider.lower()))
        files_to_compute.extend(expand("data/processed/features/{pid}/phone_battery.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
        files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

for provider in config["PHONE_SCREEN"]["PROVIDERS"].keys():
    if config["PHONE_SCREEN"]["PROVIDERS"][provider]["COMPUTE"]:
        # if "PHONE_SCREEN" in config["PHONE_DATA_YIELD"]["SENSORS"]:# not used for now because we took episodepersensedminutes out of the list of supported features
        #     files_to_compute.extend(expand("data/interim/{pid}/phone_yielded_timestamps.csv", pid=config["PIDS"]))
        # else:
        #     raise ValueError("Error: Add PHONE_SCREEN (and as many PHONE_SENSORS as you have in your database) to [PHONE_DATA_YIELD][SENSORS] in config.yaml. This is necessary to compute phone_yielded_timestamps (time when the smartphone was sensing data)")
        files_to_compute.extend(expand("data/raw/{pid}/phone_screen_raw.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/phone_screen_with_datetime.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/phone_screen_with_datetime_unified.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_screen_episodes.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_screen_episodes_resampled.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_screen_episodes_resampled_with_datetime.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_screen_features/phone_screen_{language}_{provider_key}.csv", pid=config["PIDS"], language=config["PHONE_SCREEN"]["PROVIDERS"][provider]["SRC_LANGUAGE"].lower(), provider_key=provider.lower()))
        files_to_compute.extend(expand("data/processed/features/{pid}/phone_screen.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
        files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

for provider in config["PHONE_LIGHT"]["PROVIDERS"].keys():
    if config["PHONE_LIGHT"]["PROVIDERS"][provider]["COMPUTE"]:
        files_to_compute.extend(expand("data/raw/{pid}/phone_light_raw.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/phone_light_with_datetime.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_light_features/phone_light_{language}_{provider_key}.csv", pid=config["PIDS"], language=config["PHONE_LIGHT"]["PROVIDERS"][provider]["SRC_LANGUAGE"].lower(), provider_key=provider.lower()))
        files_to_compute.extend(expand("data/processed/features/{pid}/phone_light.csv", pid=config["PIDS"],))
        files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
        files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

for provider in config["PHONE_ACCELEROMETER"]["PROVIDERS"].keys():
    if config["PHONE_ACCELEROMETER"]["PROVIDERS"][provider]["COMPUTE"]:
        files_to_compute.extend(expand("data/raw/{pid}/phone_accelerometer_raw.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/phone_accelerometer_with_datetime.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_accelerometer_features/phone_accelerometer_{language}_{provider_key}.csv", pid=config["PIDS"], language=config["PHONE_ACCELEROMETER"]["PROVIDERS"][provider]["SRC_LANGUAGE"].lower(), provider_key=provider.lower()))
        files_to_compute.extend(expand("data/processed/features/{pid}/phone_accelerometer.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
        files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

for provider in config["PHONE_APPLICATIONS_FOREGROUND"]["PROVIDERS"].keys():
    if config["PHONE_APPLICATIONS_FOREGROUND"]["PROVIDERS"][provider]["COMPUTE"]:
        files_to_compute.extend(expand("data/raw/{pid}/phone_applications_foreground_raw.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/phone_applications_foreground_with_datetime.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/phone_applications_foreground_with_datetime_with_categories.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_applications_foreground_features/phone_applications_foreground_{language}_{provider_key}.csv", pid=config["PIDS"], language=config["PHONE_APPLICATIONS_FOREGROUND"]["PROVIDERS"][provider]["SRC_LANGUAGE"].lower(), provider_key=provider.lower()))
        files_to_compute.extend(expand("data/processed/features/{pid}/phone_applications_foreground.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
        files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

for provider in config["PHONE_WIFI_VISIBLE"]["PROVIDERS"].keys():
    if config["PHONE_WIFI_VISIBLE"]["PROVIDERS"][provider]["COMPUTE"]:
        files_to_compute.extend(expand("data/raw/{pid}/phone_wifi_visible_raw.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/phone_wifi_visible_with_datetime.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_wifi_visible_features/phone_wifi_visible_{language}_{provider_key}.csv", pid=config["PIDS"], language=config["PHONE_WIFI_VISIBLE"]["PROVIDERS"][provider]["SRC_LANGUAGE"].lower(), provider_key=provider.lower()))
        files_to_compute.extend(expand("data/processed/features/{pid}/phone_wifi_visible.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
        files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

for provider in config["PHONE_WIFI_CONNECTED"]["PROVIDERS"].keys():
    if config["PHONE_WIFI_CONNECTED"]["PROVIDERS"][provider]["COMPUTE"]:
        files_to_compute.extend(expand("data/raw/{pid}/phone_wifi_connected_raw.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/phone_wifi_connected_with_datetime.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_wifi_connected_features/phone_wifi_connected_{language}_{provider_key}.csv", pid=config["PIDS"], language=config["PHONE_WIFI_CONNECTED"]["PROVIDERS"][provider]["SRC_LANGUAGE"].lower(), provider_key=provider.lower()))
        files_to_compute.extend(expand("data/processed/features/{pid}/phone_wifi_connected.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
        files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

for provider in config["PHONE_CONVERSATION"]["PROVIDERS"].keys():    
    if config["PHONE_CONVERSATION"]["PROVIDERS"][provider]["COMPUTE"]:
        files_to_compute.extend(expand("data/raw/{pid}/phone_conversation_raw.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/phone_conversation_with_datetime.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/phone_conversation_with_datetime_unified.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_conversation_features/phone_conversation_{language}_{provider_key}.csv", pid=config["PIDS"], language=config["PHONE_CONVERSATION"]["PROVIDERS"][provider]["SRC_LANGUAGE"].lower(), provider_key=provider.lower()))
        files_to_compute.extend(expand("data/processed/features/{pid}/phone_conversation.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
        files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

# We can delete these if's as soon as we add feature PROVIDERS to any of these sensors
if isinstance(config["PHONE_APPLICATIONS_CRASHES"]["PROVIDERS"], dict):
    for provider in config["PHONE_APPLICATIONS_CRASHES"]["PROVIDERS"].keys():
        if config["PHONE_APPLICATIONS_CRASHES"]["PROVIDERS"][provider]["COMPUTE"]:
            files_to_compute.extend(expand("data/raw/{pid}/phone_applications_crashes_raw.csv", pid=config["PIDS"]))
            files_to_compute.extend(expand("data/raw/{pid}/phone_applications_crashes_with_datetime.csv", pid=config["PIDS"]))
            files_to_compute.extend(expand("data/raw/{pid}/phone_applications_crashes_with_datetime_with_categories.csv", pid=config["PIDS"]))
            files_to_compute.extend(expand("data/interim/{pid}/phone_applications_crashes_features/phone_applications_crashes_{language}_{provider_key}.csv", pid=config["PIDS"], language=config["PHONE_APPLICATIONS_CRASHES"]["PROVIDERS"][provider]["SRC_LANGUAGE"].lower(), provider_key=provider.lower()))
            files_to_compute.extend(expand("data/processed/features/{pid}/phone_applications_crashes.csv", pid=config["PIDS"]))
            files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
            files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

if isinstance(config["PHONE_APPLICATIONS_NOTIFICATIONS"]["PROVIDERS"], dict):
    for provider in config["PHONE_APPLICATIONS_NOTIFICATIONS"]["PROVIDERS"].keys():
        if config["PHONE_APPLICATIONS_NOTIFICATIONS"]["PROVIDERS"][provider]["COMPUTE"]:
            files_to_compute.extend(expand("data/raw/{pid}/phone_applications_notifications_raw.csv", pid=config["PIDS"]))
            files_to_compute.extend(expand("data/raw/{pid}/phone_applications_notifications_with_datetime.csv", pid=config["PIDS"]))
            files_to_compute.extend(expand("data/raw/{pid}/phone_applications_notifications_with_datetime_with_categories.csv", pid=config["PIDS"]))
            files_to_compute.extend(expand("data/interim/{pid}/phone_applications_notifications_features/phone_applications_notifications_{language}_{provider_key}.csv", pid=config["PIDS"], language=config["PHONE_APPLICATIONS_NOTIFICATIONS"]["PROVIDERS"][provider]["SRC_LANGUAGE"].lower(), provider_key=provider.lower()))
            files_to_compute.extend(expand("data/processed/features/{pid}/phone_applications_notifications.csv", pid=config["PIDS"]))
            files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
            files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

if isinstance(config["PHONE_KEYBOARD"]["PROVIDERS"], dict):
    for provider in config["PHONE_KEYBOARD"]["PROVIDERS"].keys():    
        if config["PHONE_KEYBOARD"]["PROVIDERS"][provider]["COMPUTE"]:
            files_to_compute.extend(expand("data/raw/{pid}/phone_keyboard_raw.csv", pid=config["PIDS"]))
            files_to_compute.extend(expand("data/raw/{pid}/phone_keyboard_with_datetime.csv", pid=config["PIDS"]))
            files_to_compute.extend(expand("data/interim/{pid}/phone_keyboard_features/phone_keyboard_{language}_{provider_key}.csv", pid=config["PIDS"], language=config["PHONE_KEYBOARD"]["PROVIDERS"][provider]["SRC_LANGUAGE"].lower(), provider_key=provider.lower()))
            files_to_compute.extend(expand("data/processed/features/{pid}/phone_keyboard.csv", pid=config["PIDS"]))
            files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
            files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

if isinstance(config["PHONE_AWARE_LOG"]["PROVIDERS"], dict):
    for provider in config["PHONE_AWARE_LOG"]["PROVIDERS"].keys():    
        if config["PHONE_AWARE_LOG"]["PROVIDERS"][provider]["COMPUTE"]:
            files_to_compute.extend(expand("data/raw/{pid}/phone_aware_log_raw.csv", pid=config["PIDS"]))
            files_to_compute.extend(expand("data/raw/{pid}/phone_aware_log_with_datetime.csv", pid=config["PIDS"]))
            files_to_compute.extend(expand("data/interim/{pid}/phone_aware_log_features/phone_aware_log_{language}_{provider_key}.csv", pid=config["PIDS"], language=config["PHONE_AWARE_LOG"]["PROVIDERS"][provider]["SRC_LANGUAGE"].lower(), provider_key=provider.lower()))
            files_to_compute.extend(expand("data/processed/features/{pid}/phone_aware_log.csv", pid=config["PIDS"]))
            files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
            files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

for provider in config["PHONE_LOCATIONS"]["PROVIDERS"].keys():
    if config["PHONE_LOCATIONS"]["PROVIDERS"][provider]["COMPUTE"]:
        if config["PHONE_LOCATIONS"]["LOCATIONS_TO_USE"] in ["FUSED_RESAMPLED","ALL_RESAMPLED"]:
            if "PHONE_LOCATIONS" in config["PHONE_DATA_YIELD"]["SENSORS"]:
                files_to_compute.extend(expand("data/interim/{pid}/phone_yielded_timestamps.csv", pid=config["PIDS"]))
            else:
                raise ValueError("Error: Add PHONE_LOCATIONS (and as many PHONE_SENSORS as you have) to [PHONE_DATA_YIELD][SENSORS] in config.yaml. This is necessary to compute phone_yielded_timestamps (time when the smartphone was sensing data) which is used to resample fused location data (ALL_RESAMPLED and RESAMPLED_FUSED)")

        files_to_compute.extend(expand("data/raw/{pid}/phone_locations_raw.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_locations_processed.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_locations_processed_with_datetime.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/phone_locations_features/phone_locations_{language}_{provider_key}.csv", pid=config["PIDS"], language=config["PHONE_LOCATIONS"]["PROVIDERS"][provider]["SRC_LANGUAGE"].lower(), provider_key=provider.lower()))
        files_to_compute.extend(expand("data/processed/features/{pid}/phone_locations.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
        files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

for provider in config["FITBIT_DATA_YIELD"]["PROVIDERS"].keys():
    if config["FITBIT_DATA_YIELD"]["PROVIDERS"][provider]["COMPUTE"]:
        files_to_compute.extend(expand("data/raw/{pid}/fitbit_heartrate_intraday_raw.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/fitbit_heartrate_intraday_parsed.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/fitbit_heartrate_intraday_parsed_with_datetime.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/processed/features/{pid}/fitbit_data_yield.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
        files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

for provider in config["FITBIT_HEARTRATE_SUMMARY"]["PROVIDERS"].keys():
    if config["FITBIT_HEARTRATE_SUMMARY"]["PROVIDERS"][provider]["COMPUTE"]:
        files_to_compute.extend(expand("data/raw/{pid}/fitbit_heartrate_summary_raw.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/fitbit_heartrate_summary_parsed.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/fitbit_heartrate_summary_parsed_with_datetime.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/fitbit_heartrate_summary_features/fitbit_heartrate_summary_{language}_{provider_key}.csv", pid=config["PIDS"], language=config["FITBIT_HEARTRATE_SUMMARY"]["PROVIDERS"][provider]["SRC_LANGUAGE"].lower(), provider_key=provider.lower()))
        files_to_compute.extend(expand("data/processed/features/{pid}/fitbit_heartrate_summary.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
        files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

for provider in config["FITBIT_HEARTRATE_INTRADAY"]["PROVIDERS"].keys():
    if config["FITBIT_HEARTRATE_INTRADAY"]["PROVIDERS"][provider]["COMPUTE"]:
        files_to_compute.extend(expand("data/raw/{pid}/fitbit_heartrate_intraday_raw.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/fitbit_heartrate_intraday_parsed.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/fitbit_heartrate_intraday_parsed_with_datetime.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/fitbit_heartrate_intraday_features/fitbit_heartrate_intraday_{language}_{provider_key}.csv", pid=config["PIDS"], language=config["FITBIT_HEARTRATE_INTRADAY"]["PROVIDERS"][provider]["SRC_LANGUAGE"].lower(), provider_key=provider.lower()))
        files_to_compute.extend(expand("data/processed/features/{pid}/fitbit_heartrate_intraday.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
        files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

for provider in config["FITBIT_SLEEP_SUMMARY"]["PROVIDERS"].keys():
    if config["FITBIT_SLEEP_SUMMARY"]["PROVIDERS"][provider]["COMPUTE"]:
        files_to_compute.extend(expand("data/raw/{pid}/fitbit_sleep_summary_raw.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/fitbit_sleep_summary_parsed.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/fitbit_sleep_summary_parsed_with_datetime.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/fitbit_sleep_summary_features/fitbit_sleep_summary_{language}_{provider_key}.csv", pid=config["PIDS"], language=config["FITBIT_SLEEP_SUMMARY"]["PROVIDERS"][provider]["SRC_LANGUAGE"].lower(), provider_key=provider.lower()))
        files_to_compute.extend(expand("data/processed/features/{pid}/fitbit_sleep_summary.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
        files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

# for provider in config["FITBIT_SLEEP_INTRADAY"]["PROVIDERS"].keys():
#     if config["FITBIT_SLEEP_INTRADAY"]["PROVIDERS"][provider]["COMPUTE"]:
#         files_to_compute.extend(expand("data/raw/{pid}/fitbit_sleep_intraday_raw.csv", pid=config["PIDS"]))
#         files_to_compute.extend(expand("data/raw/{pid}/fitbit_sleep_intraday_parsed.csv", pid=config["PIDS"]))
#         files_to_compute.extend(expand("data/raw/{pid}/fitbit_sleep_intraday_parsed_with_datetime.csv", pid=config["PIDS"]))

for provider in config["FITBIT_STEPS_SUMMARY"]["PROVIDERS"].keys():
    if config["FITBIT_STEPS_SUMMARY"]["PROVIDERS"][provider]["COMPUTE"]:
        files_to_compute.extend(expand("data/raw/{pid}/fitbit_steps_summary_raw.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/fitbit_steps_summary_parsed.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/fitbit_steps_summary_parsed_with_datetime.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/fitbit_steps_summary_features/fitbit_steps_summary_{language}_{provider_key}.csv", pid=config["PIDS"], language=config["FITBIT_STEPS_SUMMARY"]["PROVIDERS"][provider]["SRC_LANGUAGE"].lower(), provider_key=provider.lower()))
        files_to_compute.extend(expand("data/processed/features/{pid}/fitbit_steps_summary.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
        files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

for provider in config["FITBIT_STEPS_INTRADAY"]["PROVIDERS"].keys():
    if config["FITBIT_STEPS_INTRADAY"]["PROVIDERS"][provider]["COMPUTE"]:
        files_to_compute.extend(expand("data/raw/{pid}/fitbit_steps_intraday_raw.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/fitbit_steps_intraday_parsed.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/raw/{pid}/fitbit_steps_intraday_parsed_with_datetime.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/interim/{pid}/fitbit_steps_intraday_features/fitbit_steps_intraday_{language}_{provider_key}.csv", pid=config["PIDS"], language=config["FITBIT_STEPS_INTRADAY"]["PROVIDERS"][provider]["SRC_LANGUAGE"].lower(), provider_key=provider.lower()))
        files_to_compute.extend(expand("data/processed/features/{pid}/fitbit_steps_intraday.csv", pid=config["PIDS"]))
        files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
        files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

# for provider in config["FITBIT_CALORIES"]["PROVIDERS"].keys():
#     if config["FITBIT_CALORIES"]["PROVIDERS"][provider]["COMPUTE"]:
#         files_to_compute.extend(expand("data/raw/{pid}/fitbit_calories_{fitbit_data_type}_raw.csv", pid=config["PIDS"], fitbit_data_type=(["json"] if config["FITBIT_CALORIES"]["TABLE_FORMAT"] == "JSON" else ["summary", "intraday"])))
#         files_to_compute.extend(expand("data/raw/{pid}/fitbit_calories_{fitbit_data_type}_parsed.csv", pid=config["PIDS"], fitbit_data_type=["summary", "intraday"]))
#         files_to_compute.extend(expand("data/raw/{pid}/fitbit_calories_{fitbit_data_type}_parsed_with_datetime.csv", pid=config["PIDS"], fitbit_data_type=["summary", "intraday"]))
#         files_to_compute.extend(expand("data/processed/features/{pid}/all_sensor_features.csv", pid=config["PIDS"]))
#         files_to_compute.append("data/processed/features/all_participants/all_sensor_features.csv")

# Visualization for Data Exploration
if config["HISTOGRAM_PHONE_DATA_YIELD"]["PLOT"]:
    files_to_compute.append("reports/data_exploration/histogram_phone_data_yield.html")

if config["HEATMAP_SENSORS_PER_MINUTE_PER_TIME_SEGMENT"]["PLOT"]:
    files_to_compute.extend(expand("reports/interim/{pid}/heatmap_sensors_per_minute_per_time_segment.html", pid=config["PIDS"]))
    files_to_compute.append("reports/data_exploration/heatmap_sensors_per_minute_per_time_segment.html")

if config["HEATMAP_SENSOR_ROW_COUNT_PER_TIME_SEGMENT"]["PLOT"]:
    files_to_compute.extend(expand("reports/interim/{pid}/heatmap_sensor_row_count_per_time_segment.html", pid=config["PIDS"]))
    files_to_compute.append("reports/data_exploration/heatmap_sensor_row_count_per_time_segment.html")

if config["HEATMAP_PHONE_DATA_YIELD_PER_PARTICIPANT_PER_TIME_SEGMENT"]["PLOT"]:
    files_to_compute.append("reports/data_exploration/heatmap_phone_data_yield_per_participant_per_time_segment.html")

if config["HEATMAP_FEATURE_CORRELATION_MATRIX"]["PLOT"]:
    files_to_compute.append("reports/data_exploration/heatmap_feature_correlation_matrix.html")


rule all:
    input:
        files_to_compute

rule clean:
    shell:
        "rm -rf data/raw/* && rm -rf data/interim/* && rm -rf data/processed/* && rm -rf reports/figures/* && rm -rf reports/*.zip && rm -rf reports/compliance/*"