rule download_dataset:
    input:
        "data/external/{pid}"
    params:
        group="AAPECS"
    output:
        "data/raw/{pid}/{sensor}.csv"
    script:
        "../src/data/download_dataset.R"