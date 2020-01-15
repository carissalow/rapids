configfile: "config.yaml"
include: "rules/packrat.snakefile"
include: "rules/preprocessing.snakefile"
include: "rules/features.snakefile"
include: "rules/reports.snakefile"

rule all:
    input:
        expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["SENSORS"]),
        expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["SENSORS"]),
        expand("data/processed/{pid}/battery_deltas.csv", pid=config["PIDS"]),
        expand("data/processed/{pid}/screen_deltas.csv", pid=config["PIDS"]),
        expand("data/processed/{pid}/plugin_google_activity_recognition_deltas.csv", pid=config["PIDS"]),
        expand("data/interim/{pid}/phone_valid_sensed_days.csv", pid=config["PIDS"]),
        expand("data/interim/{pid}/phone_sensed_bins.csv", pid=config["PIDS"]),
        expand("data/processed/{pid}/sms_{sms_type}_{day_segment}.csv",
                            pid=config["PIDS"],
                            sms_type = config["SMS"]["TYPES"],
                            day_segment = config["SMS"]["DAY_SEGMENTS"]),
        expand("data/processed/{pid}/call_{call_type}_{segment}.csv",
                            pid=config["PIDS"], 
                            call_type=config["CALLS"]["TYPES"],
                            segment = config["CALLS"]["DAY_SEGMENTS"]),
        expand("data/processed/{pid}/location_barnett.csv", pid=config["PIDS"]),
        expand("data/processed/{pid}/bluetooth_{segment}.csv",
                            pid=config["PIDS"], 
                            segment = config["BLUETOOTH"]["DAY_SEGMENTS"]),
        expand("data/processed/{pid}/google_activity_recognition_{segment}.csv",pid=config["PIDS"], 
                            segment = config["GOOGLE_ACTIVITY_RECOGNITION"]["DAY_SEGMENTS"]),
        expand("data/processed/{pid}/battery_{day_segment}.csv",
                            pid = config["PIDS"],
                            day_segment = config["BATTERY"]["DAY_SEGMENTS"]),
        expand("data/processed/{pid}/screen_{day_segment}.csv",
                            pid = config["PIDS"],
                            day_segment = config["SCREEN"]["DAY_SEGMENTS"]),
        expand("data/processed/{pid}/light_{day_segment}.csv",
                            pid = config["PIDS"],
                            day_segment = config["LIGHT"]["DAY_SEGMENTS"]),
        expand("data/processed/{pid}/accelerometer_{day_segment}.csv",
                    pid = config["PIDS"],
                    day_segment = config["ACCELEROMETER"]["DAY_SEGMENTS"]),
        # Reports
        expand("reports/figures/{pid}/{sensor}_heatmap_rows.html", pid=config["PIDS"], sensor=config["SENSORS"]),
        expand("reports/figures/{pid}/compliance_heatmap.html", pid=config["PIDS"]),
        expand("reports/figures/{pid}/battery_consumption_rates_barchart.html", pid=config["PIDS"]),
        expand("reports/compliance/{pid}/compliance_report.html", pid=config["PIDS"]),

rule clean:
    shell:
        "rm -rf data/raw/* && rm -rf data/interim/* && rm -rf data/processed/* && rm -rf reports/figures/* && rm -rf reports/*.zip && rm -rf reports/compliance/*"