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
