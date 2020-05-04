## renv_install: installs renv onto machine
rule renv_install:
    shell:
        "R -e 'if (!requireNamespace(\"renv\", quietly = TRUE)) install.packages(\"renv\", repos=\"http://cran.us.r-project.org\")'"

## renv_install: initialize a renv environment for this project
rule renv_init:
    shell:
        "R -e 'renv::init()'"

## renv_snap   : Look for new R packages in files & archives them
rule renv_snap:
    shell:
        "R -e 'renv::snapshot()'"

## renv_restore: Installs archived packages onto a new machine
rule renv_restore:
    shell:
        "R -e 'renv::restore()'"
