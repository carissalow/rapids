# File Structure

!!! tip
    - Read this page if you want to learn more about how RAPIDS is structured. If you want to start using it go to [Installation](../setup/installation/), then to [Configuration](../setup/configuration/), and then to [Execution](../setup/execution/)
    - All paths mentioned in this page are relative to RAPIDS' root folder.

If you want to extract the behavioral features that RAPIDS offers, you will only have to create or modify the [`.env` file](../setup/configuration/#database-credentials), [participants files](../setup/configuration/#participant-files), [time segment files](../setup/configuration/#time-segments), and the `config.yaml` file as instructed in the [Configuration page](../setup/configuration). The `config.yaml` file is the heart of RAPIDS and includes parameters to manage participants, data sources, sensor data, visualizations and more.


All data is saved in `data/`. The `data/external/` folder stores any data imported or created by the user, `data/raw/` stores sensor data as imported from your database, `data/interim/` has intermediate files necessary to compute behavioral features from raw data, and `data/processed/` has all the final files with the behavioral features in folders per participant and sensor.

RAPIDS source code is saved in `src/`. The `src/data/` folder stores scripts to download, clean and pre-process sensor data, `src/features` has scripts to extract behavioral features organized in their respective sensor subfolders , `src/models/` can host any script to create models or statistical analyses with the behavioral features you extract, and `src/visualization/` has scripts to create plots of the raw and processed data. There are other files and folders but only relevant if you are interested in extending RAPIDS (e.g. virtual env files, docs, tests, Dockerfile, the Snakefile, etc.). 

In the figure below, we represent the interactions between users and files. After a user modifies the configuration files mentioned above, the `Snakefile` file will search for and execute the Snakemake rules that contain the Python or R scripts necessary to generate or update the required output files (behavioral features, plots, etc.).

<figure>
  <img src="../img/files.png" width="600" />
  <figcaption>Interaction diagram between the user, and important files in RAPIDS</figcaption>
</figure>

