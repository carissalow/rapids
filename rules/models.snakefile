def input_merge_metrics_of_single_participant(wildcards):
    if wildcards.source == "phone_fitbit_metrics":
        return expand("data/processed/{pid}/{metrics}_{day_segment}.csv", pid=wildcards.pid, metrics=config["METRICS_FOR_ANALYSIS"]["PHONE_METRICS"] + config["METRICS_FOR_ANALYSIS"]["FITBIT_METRICS"], day_segment=wildcards.day_segment)
    else:
        return expand("data/processed/{pid}/{metrics}_{day_segment}.csv", pid=wildcards.pid, metrics=config["METRICS_FOR_ANALYSIS"][wildcards.source.upper()], day_segment=wildcards.day_segment)

rule merge_metrics_of_single_participant:
    input:
        metric_files = input_merge_metrics_of_single_participant
    output:
        "models/input/merged_single_participant/{pid}/{source}_{day_segment}.csv"
    script:
        "../src/models/merge_metrics_of_single_participant.R"

rule merge_metrics_of_all_participants:
    input:
        metric_files = expand("models/input/merged_single_participant/{pid}/{{source}}_{{day_segment}}.csv", pid=config["PIDS"])
    output:
        "models/input/merged_all_participants/{source}_{day_segment}.csv" 
    script:
        "../src/models/merge_metrics_of_all_participants.R"