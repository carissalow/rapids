configfile: "config.yaml"
include: "rules/packrat.snakefile"
include: "rules/preprocessing.snakefile"
include: "rules/features.snakefile"
include: "rules/models.snakefile"
include: "rules/reports.snakefile"

rule all:
    input:
        expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["SENSORS"]),
        expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["FITBIT_TABLE"]),
        expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["SENSORS"]),
        expand("data/processed/{pid}/battery_deltas.csv", pid=config["PIDS"]),
        expand("data/interim/{pid}/applications_foreground_with_datetime_with_genre.csv", pid=config["PIDS"]),
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
        expand("data/processed/{pid}/location_barnett_{segment}.csv", 
                            pid=config["PIDS"],
                            segment = config["BARNETT_LOCATION"]["DAY_SEGMENTS"]),
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
        expand("data/processed/{pid}/applications_foreground_{day_segment}.csv",
                            pid = config["PIDS"],
                            day_segment = config["APPLICATIONS_FOREGROUND"]["DAY_SEGMENTS"]),
        expand("data/raw/{pid}/fitbit_{fitbit_sensor}_with_datetime.csv",
                            pid=config["PIDS"],
                            fitbit_sensor=config["FITBIT_SENSORS"]),
        expand("data/processed/{pid}/fitbit_heartrate_{day_segment}.csv",
                            pid = config["PIDS"],
                            day_segment = config["HEARTRATE"]["DAY_SEGMENTS"]),
        expand("data/processed/{pid}/fitbit_step_{day_segment}.csv",
                            pid = config["PIDS"],
                            day_segment = config["STEP"]["DAY_SEGMENTS"]),
        # Models
        expand("models/input/merged_single_participant/{pid}/{source}_{day_segment}.csv",
                                pid = config["PIDS"],
                                source = config["METRICS_FOR_ANALYSIS"]["SOURCES"],
                                day_segment = config["METRICS_FOR_ANALYSIS"]["DAY_SEGMENTS"]),
        expand("models/input/merged_all_participants/{source}_{day_segment}.csv",
                                source = config["METRICS_FOR_ANALYSIS"]["SOURCES"],
                                day_segment = config["METRICS_FOR_ANALYSIS"]["DAY_SEGMENTS"]),
        # Reports
        expand("reports/figures/{pid}/{sensor}_heatmap_rows.html", pid=config["PIDS"], sensor=config["SENSORS"]),
        expand("reports/figures/{pid}/compliance_heatmap.html", pid=config["PIDS"]),
        expand("reports/figures/{pid}/battery_consumption_rates_barchart.html", pid=config["PIDS"]),
        expand("reports/compliance/{pid}/compliance_report.html", pid=config["PIDS"]),
        expand("reports/figures/overall_compliance_heatmap.html"),

rule clean:
    shell:
        "rm -rf data/raw/* && rm -rf data/interim/* && rm -rf data/processed/* && rm -rf reports/figures/* && rm -rf reports/*.zip && rm -rf reports/compliance/*"