rule heatmap_rows:
    input:
        sensor = "data/raw/{pid}/{sensor}_with_datetime.csv",
        pid_file = "data/external/{pid}"
    params:
        table = "{sensor}",
        pid = "{pid}",
        bin_size = config["PHONE_VALID_SENSED_BINS"]["BIN_SIZE"]
    output:
        "reports/figures/{pid}/{sensor}_heatmap_rows.html"
    script:
        "../src/visualization/heatmap_rows.py"

rule compliance_heatmap:
    input:
        sensor = "data/interim/{pid}/phone_sensed_bins.csv",
        pid_file = "data/external/{pid}"
    params:
        pid = "{pid}",
        bin_size = config["PHONE_VALID_SENSED_BINS"]["BIN_SIZE"]
    output:
        "reports/figures/{pid}/compliance_heatmap.html"
    script:
        "../src/visualization/compliance_heatmap.py"

rule overall_compliance_heatmap:
    input:
        phone_sensed_bins =  expand("data/interim/{pid}/phone_sensed_bins.csv", pid=config["PIDS"]),
        phone_valid_sensed_days = expand("data/interim/{pid}/phone_valid_sensed_days.csv", pid=config["PIDS"]),
        pid_files = expand("data/external/{pid}", pid=config["PIDS"])
    params:
        local_timezone = config["READABLE_DATETIME"]["FIXED_TIMEZONE"],
        bin_size = config["PHONE_VALID_SENSED_BINS"]["BIN_SIZE"],
        min_bins_per_hour = config["PHONE_VALID_SENSED_DAYS"]["MIN_VALID_BINS_PER_HOUR"]
    output:
        "reports/figures/overall_compliance_heatmap.html"
    script:
        "../src/visualization/overall_compliance_heatmap.py"

rule battery_consumption_rates_barchart:
    input:
        sensor = "data/processed/{pid}/battery_daily.csv",
        pid_file = "data/external/{pid}"
    params:
        pid = "{pid}"
    output:
        "reports/figures/{pid}/battery_consumption_rates_barchart.html"
    script:
        "../src/visualization/battery_consumption_rates_barchart.py"

PHONE_SENSORS = []
PHONE_SENSORS.extend([config["MESSAGES"]["DB_TABLE"], config["CALLS"]["DB_TABLE"], config["BARNETT_LOCATION"]["DB_TABLE"], config["BLUETOOTH"]["DB_TABLE"], config["BATTERY"]["DB_TABLE"], config["SCREEN"]["DB_TABLE"], config["LIGHT"]["DB_TABLE"], config["ACCELEROMETER"]["DB_TABLE"], config["APPLICATIONS_FOREGROUND"]["DB_TABLE"],config["WIFI"]["DB_TABLE"]])

rule compliance_report:
    input:
        sensor_heatmaps =  expand("reports/figures/{{pid}}/{sensor}_heatmap_rows.html", sensor=PHONE_SENSORS),
        compliance_heatmap =  rules.compliance_heatmap.output
    output:
        "reports/compliance/{pid}/compliance_report.html",
    script:
        "../src/visualization/compliance_report.Rmd"