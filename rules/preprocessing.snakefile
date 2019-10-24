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
    output:
        "data/raw/{pid}/{sensor}_with_datetime.csv"
    script:
        "../src/data/readable_datetime.R"