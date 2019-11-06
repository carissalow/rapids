rule communication_sms_metrics:
    input: 
        "data/raw/{pid}/messages_with_datetime.csv"
    params:
        sms_type = "{sms_type}",
        day_segment = "{day_segment}",
        metric = "{metric}"
    output:
        "data/processed/{pid}/com_sms_{sms_type}_{day_segment}_{metric}.csv"
    script:
        "../src/features/communication_sms_metrics.R"

rule communication_call_metrics:
    input: 
        "data/raw/{pid}/calls_with_datetime.csv"
    params:
        call_type = "{call_type}",
        day_segment = "{day_segment}",
        metric = "{metric}"
    output:
        "data/processed/{pid}/com_call_{call_type}_{day_segment}_{metric}.csv"
    script:
        "../src/features/communication_call_metrics.R"

rule battery_deltas:
    input:
        "data/raw/{pid}/battery_with_datetime.csv"
    output:
        "data/processed/{pid}/battery_deltas.csv"
    script:
        "../src/features/battery_deltas.R"

rule location_barnett_metrics:
    input:
        "data/raw/{pid}/locations_with_datetime.csv"
    params:
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
        "data/raw/{pid}/plugin_google_activity_recognition_with_datetime.csv"
    output:
        "data/processed/{pid}/google_activity_recognition.csv"
    script:
        "../src/features/google_activity_recognition.py"
