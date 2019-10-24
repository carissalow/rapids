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