rule download_dataset:
    input:
        "data/external/{pid}"
    params:
        group = config["DOWNLOAD_DATASET"]["GROUP"],
        table = "{sensor}"
    output:
        "data/raw/{pid}/{sensor}_raw.csv"
    script:
        "../src/data/download_dataset.R"

rule readable_datetime:
    input:
        sensor_input = rules.download_dataset.output
    params:
        timezones = None,
        fixed_timezone = config["READABLE_DATETIME"]["FIXED_TIMEZONE"]
    output:
        "data/raw/{pid}/{sensor}_with_datetime.csv"
    script:
        "../src/data/readable_datetime.R"

rule phone_valid_sensed_days:
    input:
        all_sensors =  expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["SENSORS"])
    params:
        bin_size = config["PHONE_VALID_SENSED_DAYS"]["BIN_SIZE"],
        min_valid_hours = config["PHONE_VALID_SENSED_DAYS"]["MIN_VALID_HOURS"],
        min_bins_per_hour = config["PHONE_VALID_SENSED_DAYS"]["MIN_BINS_PER_HOUR"]
    output:
        "data/interim/{pid}/phone_valid_sensed_days.csv"
    script:
        "../src/data/phone_valid_sensed_days.R"

rule unify_ios_android:
    input:
        sensor_data = "data/raw/{pid}/{sensor}_with_datetime.csv",
        participant_info = "data/external/{pid}"
    params:
        sensor = "{sensor}"
    output:
        "data/raw/{pid}/{sensor}_with_datetime_unified.csv"
    script:
        "../src/data/unify_ios_android.R"