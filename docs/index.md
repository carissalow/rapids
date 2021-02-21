# Welcome to RAPIDS documentation

Reproducible Analysis Pipeline for Data Streams (RAPIDS) allows you to process smartphone and wearable data to [extract](features/feature-introduction.md) and [create](features/add-new-features.md) **behavioral features** (a.k.a. digital biomarkers), [visualize](visualizations/data-quality-visualizations.md) mobile sensor data and [structure](workflow-examples/analysis.md) your analysis into reproducible workflows.

RAPIDS is open source, documented, modular, tested, and reproducible. At the moment we support smartphone data collected with [AWARE](https://awareframework.com/), wearable data from Fitbit devices, and wearable data from Empatica devices (in collaboration with the [DBDP](https://dbdp.org/)).

!!! tip
    :material-slack: Questions or feedback can be posted on the \#rapids channel in AWARE Framework\'s [slack](http://awareframework.com:3000/). 

    :material-github: Bugs and feature requests should be posted on [Github](https://github.com/carissalow/rapids/issues). 

    :fontawesome-solid-tasks: Join our discussions on our algorithms and assumptions for feature [processing](https://github.com/carissalow/rapids/discussions).

    :fontawesome-solid-play: Ready to start? Go to [Installation](setup/installation/), then to [Configuration](setup/configuration/), and then to [Execution](setup/execution/)

    :fontawesome-solid-sync-alt: Are you upgrading from RAPIDS [beta](https://rapidspitt.readthedocs.io/en/latest/)? Follow this [guide](migrating-from-old-versions)

## How does it work?

RAPIDS is formed by R and Python scripts orchestrated by [Snakemake](https://snakemake.readthedocs.io/en/stable/). We suggest you read Snakemake's docs but in short: every link in the analysis chain is atomic and has files as input and output. Behavioral features are processed per sensor and per participant.

## What are the benefits of using RAPIDS?

1. **Consistent analysis**. Every participant sensor dataset is analyzed in the exact same way and isolated from each other.
2. **Efficient analysis**. Every analysis step is executed only once. Whenever your data or configuration changes only the affected files are updated.
5. **Parallel execution**. Thanks to Snakemake, your analysis can be executed over multiple cores without changing your code.
6. **Code-free features**. Extract any of the behavioral features offered by RAPIDS without writing any code.
7. **Extensible code**. You can easily add your own behavioral features in R or Python, share them with the community, and keep authorship and citations.
8. **Timezone aware**. Your data is adjusted to the specified timezone (multiple timezones suport *coming soon*).
9. **Flexible time segments**. You can extract behavioral features on time windows of any length (e.g. 5 minutes, 3 hours, 2 days), on every day or particular days (e.g. weekends, Mondays, the 1st of each month, etc.) or around events of interest (e.g. surveys or clinical relapses).
10. **Tested code**. We are constantly adding tests to make sure our behavioral features are correct.
11. **Reproducible code**. If you structure your analysis within RAPIDS, you can be sure your code will run in other computers as intended thanks to R and Python virtual environments. You can share your analysis code along your publications without any overhead.
12. **Private**. All your data is processed locally.

## How is it organized?

In broad terms the `config.yaml`, [`.env` file](setup/configuration/#database-credentials), [participants files](setup/configuration/#participant-files), and [time segment files](setup/configuration/#time-segments) are the only ones that you will have to modify. All data is stored in `data/` and all scripts are stored in `src/`. For more information see RAPIDS' [File Structure](file-structure.md).