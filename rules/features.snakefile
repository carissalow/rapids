def optional_ar_input(wildcards):
    with open("data/external/"+wildcards.pid, encoding="ISO-8859-1") as external_file:
        external_file_content = external_file.readlines()
    platform = external_file_content[1].strip()
    if platform == "android":
        return ["data/raw/{pid}/plugin_google_activity_recognition_with_datetime_unified.csv",
                "data/processed/{pid}/plugin_google_activity_recognition_deltas.csv"]
    else:
        return ["data/raw/{pid}/plugin_ios_activity_recognition_with_datetime_unified.csv",
                "data/processed/{pid}/plugin_ios_activity_recognition_deltas.csv"]

rule sms_features:
    input: 
        "data/raw/{pid}/messages_with_datetime.csv"
    params:
        sms_type = "{sms_type}",
        day_segment = "{day_segment}",
        features = lambda wildcards: config["SMS"]["FEATURES"][wildcards.sms_type]
    output:
        "data/processed/{pid}/sms_{sms_type}_{day_segment}.csv"
    script:
        "../src/features/sms_features.R"

rule call_features:
    input: 
        "data/raw/{pid}/calls_with_datetime_unified.csv"
    params:
        call_type = "{call_type}",
        day_segment = "{day_segment}",
        features = lambda wildcards: config["CALLS"]["FEATURES"][wildcards.call_type]
    output:
        "data/processed/{pid}/call_{call_type}_{day_segment}.csv"
    script:
        "../src/features/call_features.R"

rule battery_deltas:
    input:
        "data/raw/{pid}/battery_with_datetime_unified.csv"
    output:
        "data/processed/{pid}/battery_deltas.csv"
    script:
        "../src/features/battery_deltas.R"

rule screen_deltas:
    input:
        screen = "data/raw/{pid}/screen_with_datetime.csv",
        participant_info = "data/external/{pid}"
    output:
        "data/processed/{pid}/screen_deltas.csv"
    script:
        "../src/features/screen_deltas.R"

rule google_activity_recognition_deltas:
    input:
        "data/raw/{pid}/plugin_google_activity_recognition_with_datetime_unified.csv"
    output:
        "data/processed/{pid}/plugin_google_activity_recognition_deltas.csv"
    script:
        "../src/features/activity_recognition_deltas.R"

rule ios_activity_recognition_deltas:
    input:
        "data/raw/{pid}/plugin_ios_activity_recognition_with_datetime_unified.csv"
    output:
        "data/processed/{pid}/plugin_ios_activity_recognition_deltas.csv"
    script:
        "../src/features/activity_recognition_deltas.R"

rule location_barnett_features:
    input:
        raw = "data/raw/{pid}/locations_raw.csv",
        fused = rules.resample_fused_location.output
    params:
        features = config["BARNETT_LOCATION"]["FEATURES"],
        locations_to_use = config["BARNETT_LOCATION"]["LOCATIONS_TO_USE"],
        accuracy_limit = config["BARNETT_LOCATION"]["ACCURACY_LIMIT"],
        timezone = config["BARNETT_LOCATION"]["TIMEZONE"],
        day_segment = "{day_segment}"
    output:
        "data/processed/{pid}/location_barnett_{day_segment}.csv"
    script:
        "../src/features/location_barnett_features.R"

rule bluetooth_features:
    input: 
        "data/raw/{pid}/bluetooth_with_datetime.csv"
    params:
        day_segment = "{day_segment}",
        features = config["BLUETOOTH"]["FEATURES"]
    output:
        "data/processed/{pid}/bluetooth_{day_segment}.csv"
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
        bin_size = config["PHONE_VALID_SENSED_DAYS"]["BIN_SIZE"]
    output:
        "data/processed/{pid}/screen_{day_segment}.csv"
    script:
        "../src/features/screen_features.py"

rule light_features:
    input:
        "data/raw/{pid}/light_with_datetime.csv",
    params:
        day_segment = "{day_segment}",
        features = config["LIGHT"]["FEATURES"],
    output:
        "data/processed/{pid}/light_{day_segment}.csv"
    script:
        "../src/features/light_features.py"

rule accelerometer_features:
    input:
        "data/raw/{pid}/accelerometer_with_datetime.csv",
    params:
        day_segment = "{day_segment}",
        features = config["ACCELEROMETER"]["FEATURES"],
    output:
        "data/processed/{pid}/accelerometer_{day_segment}.csv"
    script:
        "../src/features/accelerometer_features.py"

rule applications_foreground_features:
    input:
        "data/interim/{pid}/applications_foreground_with_datetime_with_genre.csv",
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
        "data/raw/{pid}/wifi_with_datetime.csv"
    params:
        day_segment = "{day_segment}",
        features = config["WIFI"]["FEATURES"]
    output:
        "data/processed/{pid}/wifi_{day_segment}.csv"
    script:
        "../src/features/wifi_features.R"

rule fitbit_heartrate_features:
    input:
        heartrate_summary_data = "data/raw/{pid}/fitbit_heartrate_summary_with_datetime.csv",
        heartrate_intraday_data = "data/raw/{pid}/fitbit_heartrate_intraday_with_datetime.csv"
    params:
        day_segment = "{day_segment}",
        features = config["HEARTRATE"]["FEATURES"],
        daily_features_from_summary_data = config["HEARTRATE"]["DAILY_FEATURES_FROM_SUMMARY_DATA"]
    output:
        "data/processed/{pid}/fitbit_heartrate_{day_segment}.csv"
    script:
        "../src/features/fitbit_heartrate_features.py"

rule fitbit_step_features:
    input:
        steps_data = "data/raw/{pid}/fitbit_steps_with_datetime.csv",
    params:
        day_segment = "{day_segment}",
        features_all_steps = config["STEP"]["FEATURES"]["ALL_STEPS"],
        features_sedentary_bout = config["STEP"]["FEATURES"]["SEDENTARY_BOUT"],
        features_active_bout = config["STEP"]["FEATURES"]["ACTIVE_BOUT"],
        threshold_active_bout = config["STEP"]["THRESHOLD_ACTIVE_BOUT"],
        include_zero_step_rows = config["STEP"]["INCLUDE_ZERO_STEP_ROWS"]
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
        sleep_types = config["SLEEP"]["SLEEP_TYPES"],
        daily_features_from_summary_data = config["SLEEP"]["DAILY_FEATURES_FROM_SUMMARY_DATA"]
    output:
        "data/processed/{pid}/fitbit_sleep_{day_segment}.csv"
    script:
        "../src/features/fitbit_sleep_features.py"
