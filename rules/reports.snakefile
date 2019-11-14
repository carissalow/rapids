rule heatmap_rows:
    input:
        "data/raw/{pid}/{sensor}_with_datetime.csv"
    params:
        table = "{sensor}",
        pid = "{pid}"
    output:
        "reports/figures/{pid}/{sensor}_heatmap_rows.html"
    script:
        "../src/visualization/heatmap_rows.py"

rule compliance_heatmap:
    input:
        expand("data/raw/{{pid}}/{sensor}_with_datetime.csv", sensor=config["SENSORS"])
    params:
        pid = "{pid}"
    output:
        "reports/figures/{pid}/compliance_heatmap.html"
    script:
        "../src/visualization/compliance_heatmap.py"

rule battery_consumption_rates_barchart:
    input:
        "data/processed/{pid}/battery_daily.csv"
    params:
        pid = "{pid}"
    output:
        "reports/figures/{pid}/battery_consumption_rates_barchart.html"
    script:
        "../src/visualization/battery_consumption_rates_barchart.py"

rule zip_aapecs_features:
    output:
        "reports/features_{date}.zip",
    shell:
        "zip -r {output} data/processed"

rule zip_aapecs_interim:
    output:
        "reports/interim_{date}.zip",
    shell:
        "zip -r {output} data/interim"

rule zip_aapecs_figures:
    output:
        "reports/figures_{date}.zip",
    shell:
        "zip -r {output} reports/figures"