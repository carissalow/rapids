configfile: "config.yaml"
include: "rules/preprocessing.snakefile"

rule all:
    input:
        expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["SENSORS"]),
        expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["SENSORS"])

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