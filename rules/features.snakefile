def optional_ar_input(wildcards):
    with open("data/external/"+wildcards.pid, encoding="ISO-8859-1") as external_file:
        external_file_content = external_file.readlines()
    platforms = external_file_content[1].strip().split(",")
    if platforms[0] == "multiple" or (len(platforms) > 1 and "android" in platforms and "ios" in platforms):
        platform = "android"
    else:
        platform = platforms[0]

    if platform == "android": 
        return ["data/raw/{pid}/" + config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["ANDROID"] + "_with_datetime_unified.csv",
                "data/processed/{pid}/" + config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["ANDROID"] + "_deltas.csv"]
    elif platform == "ios":
        return ["data/raw/{pid}/"+config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["IOS"]+"_with_datetime_unified.csv",
                "data/processed/{pid}/"+config["ACTIVITY_RECOGNITION"]["DB_TABLE"]["IOS"]+"_deltas.csv"]
    else:
        raise ValueError("Platform (line 2) in a participant file should be 'android', 'ios', or 'multiple'. You typed '" + platforms + "'")

def optional_conversation_input(wildcards):
    with open("data/external/"+wildcards.pid, encoding="ISO-8859-1") as external_file:
        external_file_content = external_file.readlines()
    platforms = external_file_content[1].strip().split(",")
    if platforms[0] == "multiple" or (len(platforms) > 1 and "android" in platforms and "ios" in platforms):
        platform = "android"
    else:
        platform = platforms[0]

    if platform == "android":
        return ["data/raw/{pid}/" + config["CONVERSATION"]["DB_TABLE"]["ANDROID"] + "_with_datetime.csv"]
    elif platform == "ios":
        return ["data/raw/{pid}/" + config["CONVERSATION"]["DB_TABLE"]["IOS"] + "_with_datetime.csv"]
    else:
        raise ValueError("Platform (line 2) in a participant file should be 'android' or 'ios', or 'multiple'. You typed '" + platforms + "'")

def optional_location_input(wildcards):
    if config["BARNETT_LOCATION"]["LOCATIONS_TO_USE"] == "RESAMPLE_FUSED":
        return expand("data/raw/{{pid}}/{sensor}_resampled.csv", sensor=config["BARNETT_LOCATION"]["DB_TABLE"])
    else:
        return expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["BARNETT_LOCATION"]["DB_TABLE"])

def optional_location_doryab_input(wildcards):
    if config["DORYAB_LOCATION"]["LOCATIONS_TO_USE"] == "RESAMPLE_FUSED":
        return expand("data/raw/{{pid}}/{sensor}_resampled.csv", sensor=config["DORYAB_LOCATION"]["DB_TABLE"])
    else:
        return expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["DORYAB_LOCATION"]["DB_TABLE"])

def optional_steps_sleep_input(wildcards):
    if config["STEP"]["EXCLUDE_SLEEP"]["EXCLUDE"] == True and config["STEP"]["EXCLUDE_SLEEP"]["TYPE"] == "FITBIT_BASED":
        return  "data/raw/{pid}/fitbit_sleep_summary_with_datetime.csv"
    else:
        return []

rule messages_features:
    input: 
        expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["MESSAGES"]["DB_TABLE"])
    params:
        messages_type = "{messages_type}",
        day_segment = "{day_segment}",
        features = lambda wildcards: config["MESSAGES"]["FEATURES"][wildcards.messages_type]
    output:
        "data/processed/{pid}/messages_{messages_type}_{day_segment}.csv"
    script:
        "../src/features/messages_features.R"

rule call_features:
    input: 
        expand("data/raw/{{pid}}/{sensor}_with_datetime_unified.csv", sensor=config["CALLS"]["DB_TABLE"])
    params:
        call_type = "{call_type}",
        day_segment = "{day_segment}",
        features = lambda wildcards: config["CALLS"]["FEATURES"][wildcards.call_type]
    output:
        "data/processed/{pid}/calls_{call_type}_{day_segment}.csv"
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

rule location_barnett_features:
    input:
        locations = optional_location_input
    params:
        features = config["BARNETT_LOCATION"]["FEATURES"],
        locations_to_use = config["BARNETT_LOCATION"]["LOCATIONS_TO_USE"],
        accuracy_limit = config["BARNETT_LOCATION"]["ACCURACY_LIMIT"],
        timezone = config["BARNETT_LOCATION"]["TIMEZONE"],
        minutes_data_used = config["BARNETT_LOCATION"]["MINUTES_DATA_USED"],
        day_segment = "{day_segment}"
    output:
        "data/processed/{pid}/location_barnett_{day_segment}.csv"
    script:
        "../src/features/location_barnett_features.R"

rule location_doryab_features:
    input:
        locations = optional_location_doryab_input
    params:
        features = config["DORYAB_LOCATION"]["FEATURES"],
        day_segment = "{day_segment}",
        dbscan_eps = config["DORYAB_LOCATION"]["DBSCAN_EPS"],
        dbscan_minsamples = config["DORYAB_LOCATION"]["DBSCAN_MINSAMPLES"],
        threshold_static = config["DORYAB_LOCATION"]["THRESHOLD_STATIC"],
        maximum_gap_allowed = config["DORYAB_LOCATION"]["MAXIMUM_GAP_ALLOWED"]
    output:
        "data/processed/{pid}/location_doryab_{day_segment}.csv"
    script:
        "../src/features/location_doryab_features.py"

rule bluetooth_features:
    input: 
        expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["BLUETOOTH"]["DB_TABLE"])
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
        expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["WIFI"]["DB_TABLE"])
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
