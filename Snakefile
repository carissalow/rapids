configfile: "config.yaml"
include: "rules/preprocessing.snakefile"
include: "rules/features.snakefile"
include: "rules/reports.snakefile"

rule all:
    input:
        expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["SENSORS"]),
        expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["SENSORS"]),
        expand("data/processed/{pid}/battery_deltas.csv", pid=config["PIDS"]),
        expand("data/interim/{pid}/phone_valid_sensed_days.csv", pid=config["PIDS"]),
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
        # Reports
        expand("reports/figures/{pid}/{sensor}_heatmap_rows.html", pid=config["PIDS"], sensor=config["SENSORS"]),
        expand("reports/figures/{pid}/compliance_heatmap.html", pid=config["PIDS"], sensor=config["SENSORS"]),
        expand("reports/figures/{pid}/battery_consumption_rates_barchart.html", pid=config["PIDS"]),

# --- Packrat Rules --- #
## Taken from https://github.com/lachlandeer/snakemake-econ-r

## packrat_install: installs packrat onto machine
rule packrat_install:
    shell:
        "R -e 'install.packages(\"packrat\", repos=\"http://cran.us.r-project.org\")'"

## packrat_install: initialize a packrat environment for this project
rule packrat_init:
    shell:
        "R -e 'packrat::init()'"

## packrat_snap   : Look for new R packages in files & archives them
rule packrat_snap:
    shell:
        "R -e 'packrat::snapshot()'"

## packrat_restore: Installs archived packages onto a new machine
rule packrat_restore:
    shell:
        "R -e 'packrat::restore()'"