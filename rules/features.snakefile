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