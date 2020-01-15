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

rule call_metrics:
    input: 
        "data/raw/{pid}/calls_with_datetime_unified.csv"
    params:
        call_type = "{call_type}",
        day_segment = "{day_segment}",
        metrics = lambda wildcards: config["CALLS"]["METRICS"][wildcards.call_type]
    output:
        "data/processed/{pid}/call_{call_type}_{day_segment}.csv"
    script:
        "../src/features/call_metrics.R"

rule battery_deltas:
    input:
        "data/raw/{pid}/battery_with_datetime_unified.csv"
    output:
        "data/processed/{pid}/battery_deltas.csv"
    script:
        "../src/features/battery_deltas.R"

rule screen_deltas:
    input:
        "data/raw/{pid}/screen_with_datetime.csv"
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
        locations_to_use = config["BARNETT_LOCATION"]["LOCATIONS_TO_USE"],
        accuracy_limit = config["BARNETT_LOCATION"]["ACCURACY_LIMIT"],
        timezone = config["BARNETT_LOCATION"]["TIMEZONE"]
    output:
        "data/processed/{pid}/location_barnett.csv"
    script:
        "../src/features/location_barnett_metrics.R"

rule bluetooth_metrics:
    input: 
        "data/raw/{pid}/bluetooth_with_datetime.csv"
    params:
        day_segment = "{day_segment}",
        metrics = config["BLUETOOTH"]["METRICS"]
    output:
        "data/processed/{pid}/bluetooth_{day_segment}.csv"
    script:
        "../src/features/bluetooth_metrics.R"
        
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

rule screen_metrics:
    input:
        screen_events = "data/raw/{pid}/screen_with_datetime.csv",
        screen_deltas = "data/processed/{pid}/screen_deltas.csv",
        phone_sensed_bins = "data/interim/{pid}/phone_sensed_bins.csv"
    params:
        day_segment = "{day_segment}",
        metrics_events = config["SCREEN"]["METRICS_EVENTS"],
        metrics_deltas = config["SCREEN"]["METRICS_DELTAS"],
        episodes = config["SCREEN"]["EPISODES"],
        bin_size = config["PHONE_VALID_SENSED_DAYS"]["BIN_SIZE"]
    output:
        "data/processed/{pid}/screen_{day_segment}.csv"
    script:
        "../src/features/screen_metrics.py"

rule light_metrics:
    input:
        "data/raw/{pid}/light_with_datetime.csv",
    params:
        day_segment = "{day_segment}",
        metrics = config["LIGHT"]["METRICS"],
    output:
        "data/processed/{pid}/light_{day_segment}.csv"
    script:
        "../src/features/light_metrics.py"

rule accelerometer_metrics:
    input:
        "data/raw/{pid}/accelerometer_with_datetime.csv",
    params:
        day_segment = "{day_segment}",
        metrics = config["ACCELEROMETER"]["METRICS"],
    output:
        "data/processed/{pid}/accelerometer_{day_segment}.csv"
    script:
        "../src/features/accelerometer_metrics.py"
