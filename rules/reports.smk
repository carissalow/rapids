rule histogram_phone_data_yield:
    input:
        "data/processed/features/all_participants/all_sensor_features.csv"
    params:
        time_segments_type = config["TIME_SEGMENTS"]["TYPE"]
    output:
        "reports/data_exploration/histogram_phone_data_yield.html"
    script:
        "../src/visualization/histogram_phone_data_yield.py"

rule heatmap_sensors_per_minute_per_time_segment:
    input:
        phone_data_yield = "data/interim/{pid}/phone_yielded_timestamps_with_datetime.csv",
        participant_file = "data/external/participant_files/{pid}.yaml",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        pid = "{pid}",
        time_segments_type = config["TIME_SEGMENTS"]["TYPE"]
    output:
        "reports/interim/{pid}/heatmap_sensors_per_minute_per_time_segment.html"
    script:
        "../src/visualization/heatmap_sensors_per_minute_per_time_segment.py"

rule merge_heatmap_sensors_per_minute_per_time_segment:
    input:
        heatmap_sensors_per_minute_per_time_segment = expand("reports/interim/{pid}/heatmap_sensors_per_minute_per_time_segment.html", pid=config["PIDS"])
    output:
        "reports/data_exploration/heatmap_sensors_per_minute_per_time_segment.html"
    script:
        "../src/visualization/merge_heatmap_sensors_per_minute_per_time_segment.Rmd"

rule heatmap_sensor_row_count_per_time_segment:
    input:
        all_sensors = expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor = map(str.lower, config["HEATMAP_SENSOR_ROW_COUNT_PER_TIME_SEGMENT"]["SENSORS"])),
        phone_data_yield = "data/processed/features/{pid}/phone_data_yield.csv",
        participant_file = "data/external/participant_files/{pid}.yaml",
        time_segments_labels = "data/interim/time_segments/{pid}_time_segments_labels.csv"
    params:
        pid = "{pid}",
        sensor_names = config["HEATMAP_SENSOR_ROW_COUNT_PER_TIME_SEGMENT"]["SENSORS"],
        time_segments_type = config["TIME_SEGMENTS"]["TYPE"]
    output:
        "reports/interim/{pid}/heatmap_sensor_row_count_per_time_segment.html"
    script:
        "../src/visualization/heatmap_sensor_row_count_per_time_segment.py"

rule merge_heatmap_sensor_row_count_per_time_segment:
    input:
        heatmap_sensor_row_count_per_time_segment = expand("reports/interim/{pid}/heatmap_sensor_row_count_per_time_segment.html", pid=config["PIDS"])
    output:
        "reports/data_exploration/heatmap_sensor_row_count_per_time_segment.html"
    script:
        "../src/visualization/merge_heatmap_sensor_row_count_per_time_segment.Rmd"

rule heatmap_phone_data_yield_per_participant_per_time_segment:
    input:
        participant_files = expand("data/external/participant_files/{pid}.yaml", pid=config["PIDS"]),
        time_segments_file = config["TIME_SEGMENTS"]["FILE"],
        phone_data_yield = "data/processed/features/all_participants/all_sensor_features.csv",
    params:
        pids = config["PIDS"],
        time = config["HEATMAP_PHONE_DATA_YIELD_PER_PARTICIPANT_PER_TIME_SEGMENT"]["TIME"],
        time_segments_type = config["TIME_SEGMENTS"]["TYPE"]
    output:
        "reports/data_exploration/heatmap_phone_data_yield_per_participant_per_time_segment.html"
    script:
        "../src/visualization/heatmap_phone_data_yield_per_participant_per_time_segment.py"

rule heatmap_feature_correlation_matrix:
    input:
        all_sensor_features = "data/processed/features/all_participants/all_sensor_features.csv" # before data cleaning
    params:
        time_segments_type = config["TIME_SEGMENTS"]["TYPE"],
        min_rows_ratio = config["HEATMAP_FEATURE_CORRELATION_MATRIX"]["MIN_ROWS_RATIO"],
        corr_threshold = config["HEATMAP_FEATURE_CORRELATION_MATRIX"]["CORR_THRESHOLD"],
        corr_method = config["HEATMAP_FEATURE_CORRELATION_MATRIX"]["CORR_METHOD"]
    output:
        "reports/data_exploration/heatmap_feature_correlation_matrix.html"
    script:
        "../src/visualization/heatmap_feature_correlation_matrix.py"

