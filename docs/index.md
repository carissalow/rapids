# Welcome to RAPIDS documentation

Reproducible Analysis Pipeline for Data Streams (RAPIDS) allows you to process smartphone and wearable data to [extract](features/feature-introduction.md) and [create](features/add-new-features.md) **behavioral features** (a.k.a. digital biomarkers), [visualize](visualizations/data-quality-visualizations.md) mobile sensor data, and [structure](workflow-examples/analysis.md) your analysis into reproducible workflows.

RAPIDS is open source, documented, modular, tested, and reproducible. At the moment, we support [data streams](../../datastreams/data-streams-introduction) logged by smartphones, Fitbit wearables, and, in collaboration with the [DBDP](https://dbdp.org/), Empatica wearables (but you can [add your own]((../../datastreams/add-new-data-streams)) too). 

**If you want to know more head over to [Where do I start?](../setup/where-do-i-start/)**

!!! tip
    :material-slack: Questions or feedback can be posted on the \#rapids channel in AWARE Framework\'s [slack](http://awareframework.com:3000/). 

    :material-github: Bugs and feature requests should be posted on [Github](https://github.com/carissalow/rapids/issues). 

    :fontawesome-solid-tasks: Join our discussions on our algorithms and assumptions for feature [processing](https://github.com/carissalow/rapids/discussions).

    :fontawesome-solid-play: Ready to start? Go to [Installation](setup/installation/), then to [Configuration](setup/configuration/), and then to [Execution](setup/execution/)

    :fontawesome-solid-sync-alt: Are you upgrading from RAPIDS `0.4.x` or older? Follow this [guide](migrating-from-old-versions)

## What are the benefits of using RAPIDS?

1. **Consistent analysis**. Every participant sensor dataset is analyzed in the same way and isolated from each other.
2. **Efficient analysis**. Every analysis step is executed only once. Whenever your data or configuration changes, only the affected files are updated.
5. **Parallel execution**. Thanks to Snakemake, your analysis can be executed over multiple cores without changing your code.
6. **Code-free features**. Extract any of the behavioral features offered by RAPIDS without writing any code.
7. **Extensible code**. You can easily add your own data streams or behavioral features in R or Python, share them with the community, and keep authorship and citations.
8. **Timezone aware**. Your data is adjusted to one or more time zones per participant.
9. **Flexible time segments**. You can extract behavioral features on time windows of any length (e.g., 5 minutes, 3 hours, 2 days), on every day or particular days (e.g., weekends, Mondays, the 1st of each month, etc.), or around events of interest (e.g., surveys or clinical relapses).
10. **Tested code**. We are continually adding tests to make sure our behavioral features are correct.
11. **Reproducible code**. If you structure your analysis within RAPIDS, you can be sure your code will run in other computers as intended, thanks to R and Python virtual environments. You can share your analysis code along with your publications without any overhead.
12. **Private**. All your data is processed locally.

