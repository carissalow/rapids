def optional_heatmap_days_by_sensors_input(wildcards):
    with open("data/external/"+wildcards.pid, encoding="ISO-8859-1") as external_file:
        external_file_content = external_file.readlines()
    platforms = external_file_content[1].strip().split(",")
    if platforms[0] == "multiple" or (len(platforms) > 1 and "android" in platforms and "ios" in platforms):
        platform = "android"
    else:
        platform = platforms[0]
    
    if platform not in ["android", "ios"]:
        raise ValueError("Platform (line 2) in a participant file should be 'android', 'ios', or 'multiple'. You typed '" + platforms + "'")

    input_for_heatmap_days_by_sensors, tables = [], config["HEATMAP_DAYS_BY_SENSORS"]["DB_TABLES"]

    for sensor in ["ACTIVITY_RECOGNITION", "CONVERSATION"]:
        table = config[sensor]["DB_TABLE"][platform.upper()]
        if table in tables:
            input_for_heatmap_days_by_sensors.append("data/raw/{pid}/" + table + "_with_datetime.csv")
            tables = [x for x in tables if x not in config[sensor]["DB_TABLE"].values()]
    for table in tables:
        input_for_heatmap_days_by_sensors.append("data/raw/{pid}/" + table + "_with_datetime.csv")

    return input_for_heatmap_days_by_sensors

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
        local_timezone = config["READABLE_DATETIME"]["FIXED_TIMEZONE"],
        expected_num_of_days = config["OVERALL_COMPLIANCE_HEATMAP"]["EXPECTED_NUM_OF_DAYS"],
        bin_size = config["OVERALL_COMPLIANCE_HEATMAP"]["BIN_SIZE"],
        min_bins_per_hour = "{min_valid_bins_per_hour}"
    output:
        "reports/data_exploration/{min_valid_hours_per_day}hours_{min_valid_bins_per_hour}bins/overall_compliance_heatmap.html"
    script:
        "../src/visualization/overall_compliance_heatmap.py"

# rule heatmap_rows:
#     input:
#         sensor = "data/raw/{pid}/{sensor}_with_datetime.csv",
#         pid_file = "data/external/{pid}"
#     params:
#         table = "{sensor}",
#         pid = "{pid}",
#         bin_size = config["PHONE_VALID_SENSED_BINS"]["BIN_SIZE"]
#     output:
#         "reports/figures/{pid}/{sensor}_heatmap_rows.html"
#     script:
#         "../src/visualization/heatmap_rows.py"

# rule battery_consumption_rates_barchart:
#     input:
#         sensor = "data/processed/{pid}/battery_daily.csv",
#         pid_file = "data/external/{pid}"
#     params:
#         pid = "{pid}"
#     output:
#         "reports/figures/{pid}/battery_consumption_rates_barchart.html"
#     script:
#         "../src/visualization/battery_consumption_rates_barchart.py"

# rule compliance_report:
#     input:
#         sensor_heatmaps =  expand("reports/figures/{{pid}}/{sensor}_heatmap_rows.html", sensor=PHONE_SENSORS),
#         compliance_heatmap =  rules.compliance_heatmap.output
#     output:
#         "reports/compliance/{pid}/compliance_report.html",
#     script:
#         "../src/visualization/compliance_report.Rmd"