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
    wildcard_constraints:
        sensor = '(' + '|'.join([re.escape(x) for x in config["SENSORS"]]) + ')' # only process smartphone sensors, not fitbit
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

rule phone_sensed_bins:
    input:
        all_sensors =  expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["SENSORS"])
    params:
        bin_size = config["PHONE_VALID_SENSED_DAYS"]["BIN_SIZE"]
    output:
        "data/interim/{pid}/phone_sensed_bins.csv"
    script:
        "../src/data/phone_sensed_bins.R"

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

rule resample_fused_location:
    input:
        locations = "data/raw/{pid}/locations_raw.csv",
        phone_sensed_bins = rules.phone_sensed_bins.output
    params:
        bin_size = config["PHONE_VALID_SENSED_DAYS"]["BIN_SIZE"],
        timezone = config["RESAMPLE_FUSED_LOCATION"]["TIMEZONE"],
        consecutive_threshold = config["RESAMPLE_FUSED_LOCATION"]["CONSECUTIVE_THRESHOLD"],
        time_since_valid_location = config["RESAMPLE_FUSED_LOCATION"]["TIME_SINCE_VALID_LOCATION"]
    output:
        "data/raw/{pid}/locations_resampled.csv"
    script:
        "../src/data/resample_fused_location.R"

rule application_genres:
    input:
        "data/raw/{pid}/applications_foreground_with_datetime.csv"
    params:
        catalogue_source = config["APPLICATION_GENRES"]["CATALOGUE_SOURCE"],
        catalogue_file = config["APPLICATION_GENRES"]["CATALOGUE_FILE"],
        update_catalogue_file = config["APPLICATION_GENRES"]["UPDATE_CATALOGUE_FILE"],
        scrape_missing_genres = config["APPLICATION_GENRES"]["SCRAPE_MISSING_GENRES"]
    output:
        "data/interim/{pid}/applications_foreground_with_datetime_with_genre.csv"
    script:
        "../src/data/application_genres.R"

rule fitbit_with_datetime:
    input:
        "data/raw/{pid}/fitbit_data_raw.csv"
    params:
        local_timezone = config["READABLE_DATETIME"]["FIXED_TIMEZONE"],
        fitbit_sensor = "{fitbit_sensor}"
    output:
        "data/raw/{pid}/fitbit_{fitbit_sensor}_with_datetime.csv"
    script:
        "../src/data/fitbit_readable_datetime.py"

