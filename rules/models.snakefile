def input_merge_metrics_of_single_participant(wildcards):
    if wildcards.source == "phone_fitbit_metrics":
        return expand("data/processed/{pid}/{metrics}_{day_segment}.csv", pid=wildcards.pid, metrics=config["METRICS_FOR_ANALYSIS"]["PHONE_METRICS"] + config["METRICS_FOR_ANALYSIS"]["FITBIT_METRICS"], day_segment=wildcards.day_segment)
    else:
        return expand("data/processed/{pid}/{metrics}_{day_segment}.csv", pid=wildcards.pid, metrics=config["METRICS_FOR_ANALYSIS"][wildcards.source.upper()], day_segment=wildcards.day_segment)

rule merge_metrics_for_individual_model:
    input:
        metric_files = input_merge_metrics_of_single_participant,
        phone_valid_sensed_days = "data/interim/{pid}/phone_valid_sensed_days.csv"
    params:
        drop_valid_sensed_days = config["METRICS_FOR_ANALYSIS"]["DROP_VALID_SENSED_DAYS"],
        source = "{source}"
    output:
        "data/processed/{pid}/metrics_for_individual_model/{source}_{day_segment}.csv"
    script:
        "../src/models/merge_metrics_for_individual_model.R"

rule merge_metrics_for_population_model:
    input:
        metric_files = expand("data/processed/{pid}/metrics_for_individual_model/{{source}}_{{day_segment}}.csv", pid=config["PIDS"])
    output:
        "data/processed/metrics_for_population_model/{source}_{day_segment}.csv" 
    script:
        "../src/models/merge_metrics_for_population_model.R"