# Welcome to RAPIDS documentation

Reproducible Analysis Pipeline for Data Streams (RAPIDS) allows you to process smartphone and wearable data to extract **behavioral features** (a.k.a. digital biomarkers/phenotypes).

RAPIDS is open source, documented, modular, tested, and reproducible. At the moment we support smartphone data collected with [AWARE](awareframework.com/) and wearable data from Fitbit devices.

:material-slack: Questions or feedback can be posted on \#rapids in AWARE Framework\'s [slack](http://awareframework.com:3000/). 

:material-github: Bugs should be reported on [Github](https://github.com/carissalow/rapids/issues). 

:fontawesome-solid-tasks: Join our discussions on our algorithms and assumptions for feature [processing](https://github.com/carissalow/rapids/issues?q=is%3Aissue+is%3Aopen+label%3Adiscussion).

## How does it work?

RAPIDS is formed by R and Python scripts orchestrated by [Snakemake](https://snakemake.readthedocs.io/en/stable/). We suggest you read Snakemake's docs but in short: every link in the analysis chain is atomic and has files as input and output. Behavioral features are processed per sensor and per participant.

## What are the benefits of using RAPIDS?

1. **Consistent analysis**. Every participant sensor dataset is analyzed in the exact same way and isolated from each other.
2. **Efficient analysis**. Every analysis step is executed only once. Whenever your data or configuration changes only the affected files are updated.
5. **Parallel execution**. Thanks to Snakemake, your analysis can be executed over multiple cores without changing your code.
6. **Extensible code**. You can easily add your own behavioral features in R or Python and keep authorship and citations.
3. **Timezone aware**. Your data is adjusted to the specified timezone (multiple timezones suport *coming soon*).
4. **Flexible day segments**. You can extract behavioral features on time windows of any length (e.g. 5 minutes, 3 hours, 2 days), on every day or particular days (e.g. weekends, Mondays, the 1st of each month, etc.) or around events of interest (e.g. surveys or clinical relapses).
7. **Tested code**. We are constantly adding tests to make sure our behavioral features are correct.
8. **Reproducible code**. You can be sure your code will run in other computers as intended thanks to R and Python virtual environments. You can share your analysis code along your publications without any overhead.
9. **Private**. All your data is processed locally.




## How is it organized?

The `config.yaml` file is the only file that you will have to modify. It includes parameters to manage participants, data sources, sensor data, visualizations and more.

All data is saved in `data/`. The `data/external/` folder stores any data imported by the user, `data/raw/` stores sensor data as imported from your database, `data/interim/` has intermediate files necessary to compute behavioral features from raw data, and `data/processed/` has all the final files with the behavioral features per sensor and participant.

All the source code is saved in `src/`. The `src/data/` folder stores scripts to download, clean and pre-process sensor data, `src/features` has scripts to extract behavioral features organized in their respective subfolders , `src/models/` can host any script to create models or statistical analyses with the behavioral features you extract, and `src/visualization/` has scripts to create plots of the raw and processed data.

There are other important files and folders but only relevant if you are interested in extending RAPIDS (e.g. virtual env files, docs, tests, Dockerfile, the Snakefile, etc.).