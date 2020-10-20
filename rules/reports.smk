rule heatmap_features_correlations:
    input:
        features = expand("data/processed/{pid}/{sensor}_{day_segment}.csv", pid=config["PIDS"], sensor=config["HEATMAP_FEATURES_CORRELATIONS"]["PHONE_FEATURES"]+config["HEATMAP_FEATURES_CORRELATIONS"]["FITBIT_FEATURES"], day_segment=config["DAY_SEGMENTS"]),
        phone_valid_sensed_days = expand("data/interim/{pid}/phone_valid_sensed_days_{{min_valid_hours_per_day}}hours_{{min_valid_bins_per_hour}}bins.csv", pid=config["PIDS"])
    params:
        min_rows_ratio = config["HEATMAP_FEATURES_CORRELATIONS"]["MIN_ROWS_RATIO"],
        corr_threshold = config["HEATMAP_FEATURES_CORRELATIONS"]["CORR_THRESHOLD"],
        corr_method = config["HEATMAP_FEATURES_CORRELATIONS"]["CORR_METHOD"]
    output:
        "reports/data_exploration/{min_valid_hours_per_day}hours_{min_valid_bins_per_hour}bins/heatmap_features_correlations.html"
    script:
        "../src/visualization/heatmap_features_correlations.py"

rule histogram_valid_sensed_hours:
    input:
        phone_valid_sensed_days = expand("data/interim/{pid}/phone_valid_sensed_days_{{min_valid_hours_per_day}}hours_{{min_valid_bins_per_hour}}bins.csv", pid=config["PIDS"])
    output:
        "reports/data_exploration/{min_valid_hours_per_day}hours_{min_valid_bins_per_hour}bins/histogram_valid_sensed_hours.html"
    script:
         "../src/visualization/histogram_valid_sensed_hours.py"

rule heatmap_days_by_sensors:
    input:
        sensors = optional_heatmap_days_by_sensors_input,
        phone_valid_sensed_days = "data/interim/{pid}/phone_valid_sensed_days_{min_valid_hours_per_day}hours_{min_valid_bins_per_hour}bins.csv"
    params:
        pid = "{pid}",
        expected_num_of_days = config["HEATMAP_DAYS_BY_SENSORS"]["EXPECTED_NUM_OF_DAYS"]
    output:
        "reports/interim/{min_valid_hours_per_day}hours_{min_valid_bins_per_hour}bins/{pid}/heatmap_days_by_sensors.html"
    script:
        "../src/visualization/heatmap_days_by_sensors.py"

rule heatmap_days_by_sensors_all_participants:
    input:
        heatmap_rows =  expand("reports/interim/{{min_valid_hours_per_day}}hours_{{min_valid_bins_per_hour}}bins/{pid}/heatmap_days_by_sensors.html", pid=config["PIDS"])
    output:
        "reports/data_exploration/{min_valid_hours_per_day}hours_{min_valid_bins_per_hour}bins/heatmap_days_by_sensors_all_participants.html"
    script:
        "../src/visualization/heatmap_days_by_sensors_all_participants.Rmd"

rule heatmap_sensed_bins:
    input:
        sensor = "data/interim/{pid}/phone_sensed_bins.csv",
        pid_file = "data/external/{pid}"
    params:
        pid = "{pid}",
        bin_size = config["HEATMAP_SENSED_BINS"]["BIN_SIZE"]
    output:
        "reports/interim/heatmap_sensed_bins/{pid}/heatmap_sensed_bins.html"
    script:
        "../src/visualization/heatmap_sensed_bins.py"

rule heatmap_sensed_bins_all_participants:
    input:
        heatmap_sensed_bins = expand("reports/interim/heatmap_sensed_bins/{pid}/heatmap_sensed_bins.html", pid=config["PIDS"])
    output:
        "reports/data_exploration/heatmap_sensed_bins_all_participants.html"
    script:
        "../src/visualization/heatmap_sensed_bins_all_participants.Rmd"

rule overall_compliance_heatmap:
    input:
        phone_sensed_bins =  expand("data/interim/{pid}/phone_sensed_bins.csv", pid=config["PIDS"]),
        phone_valid_sensed_days = expand("data/interim/{pid}/phone_valid_sensed_days_{{min_valid_hours_per_day}}hours_{{min_valid_bins_per_hour}}bins.csv", pid=config["PIDS"]),
        pid_files = expand("data/external/{pid}", pid=config["PIDS"])
    params:
        only_show_valid_days = config["OVERALL_COMPLIANCE_HEATMAP"]["ONLY_SHOW_VALID_DAYS"],
        local_timezone = config["SENSOR_DATA"]["PHONE"]["TIMEZONE"]["VALUE"],
        expected_num_of_days = config["OVERALL_COMPLIANCE_HEATMAP"]["EXPECTED_NUM_OF_DAYS"],
        bin_size = config["OVERALL_COMPLIANCE_HEATMAP"]["BIN_SIZE"],
        min_bins_per_hour = "{min_valid_bins_per_hour}"
    output:
        "reports/data_exploration/{min_valid_hours_per_day}hours_{min_valid_bins_per_hour}bins/overall_compliance_heatmap.html"
    script:
        "../src/visualization/overall_compliance_heatmap.py"
