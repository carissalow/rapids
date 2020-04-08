rule sms_metrics:
    input: 
        "data/raw/{pid}/messages_with_datetime.csv"
    params:
        sms_type = "{sms_type}",
        day_segment = "{day_segment}",
        metrics = lambda wildcards: config["SMS"]["METRICS"][wildcards.sms_type]
    output:
        "data/processed/{pid}/sms_{sms_type}_{day_segment}.csv"
    script:
        "../src/features/sms_metrics.R"

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
        "data/raw/{pid}/plugin_google_activity_recognition_with_datetime.csv"
    output:
        "data/processed/{pid}/plugin_google_activity_recognition_deltas.csv"
    script:
        "../src/features/google_activity_recognition_deltas.R"

rule location_barnett_metrics:
    input:
        raw = "data/raw/{pid}/locations_raw.csv",
        fused = rules.resample_fused_location.output
    params:
        metrics = config["BARNETT_LOCATION"]["METRICS"],
        locations_to_use = config["BARNETT_LOCATION"]["LOCATIONS_TO_USE"],
        accuracy_limit = config["BARNETT_LOCATION"]["ACCURACY_LIMIT"],
        timezone = config["BARNETT_LOCATION"]["TIMEZONE"],
        day_segment = "{day_segment}"
    output:
        "data/processed/{pid}/location_barnett_{day_segment}.csv"
    script:
        "../src/features/location_barnett_metrics.R"

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
        
rule activity_metrics:
    input:
        gar_events = "data/raw/{pid}/plugin_google_activity_recognition_with_datetime.csv",
        gar_deltas = "data/processed/{pid}/plugin_google_activity_recognition_deltas.csv"
    params:
        segment = "{day_segment}",
        metrics = config["GOOGLE_ACTIVITY_RECOGNITION"]["METRICS"]
    output:
        "data/processed/{pid}/google_activity_recognition_{day_segment}.csv"
    script:
        "../src/features/google_activity_recognition.py"

rule battery_metrics:
    input:
        "data/processed/{pid}/battery_deltas.csv"
    params:
        day_segment = "{day_segment}",
        metrics = config["BATTERY"]["METRICS"]
    output:
        "data/processed/{pid}/battery_{day_segment}.csv"
    script:
        "../src/features/battery_metrics.py"

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

rule fitbit_heartrate_metrics:
    input:
        "data/raw/{pid}/fitbit_heartrate_with_datetime.csv",
    params:
        day_segment = "{day_segment}",
        metrics = config["HEARTRATE"]["METRICS"],
    output:
        "data/processed/{pid}/fitbit_heartrate_{day_segment}.csv"
    script:
        "../src/features/fitbit_heartrate_metrics.py"

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
