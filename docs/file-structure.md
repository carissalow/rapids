# File Structure

!!! tip
    Read this page if you want to learn more about how RAPIDS is structured. If you want to start using it go to [Installation](https://www.rapids.science/setup/installation/) and then to [Initial Configuration](https://www.rapids.science/setup/configuration/)

All paths mentioned in this page are relative to RAPIDS' root folder.

If you want to extract the behavioral features that RAPIDS offers, you will only have to create or modify the [`.env` file](https://www.rapids.science/setup/configuration/#database-credentials), [participants files](https://www.rapids.science/setup/configuration/#participant-files), [day segment files](https://www.rapids.science/setup/configuration/#day-segments), and the `config.yaml` file. The `config.yaml` file is the heart of RAPIDS and includes parameters to manage participants, data sources, sensor data, visualizations and more.


All data is saved in `data/`. The `data/external/` folder stores any data imported or created by the user, `data/raw/` stores sensor data as imported from your database, `data/interim/` has intermediate files necessary to compute behavioral features from raw data, and `data/processed/` has all the final files with the behavioral features in folders per participant and sensor.

All the source code is saved in `src/`. The `src/data/` folder stores scripts to download, clean and pre-process sensor data, `src/features` has scripts to extract behavioral features organized in their respective subfolders , `src/models/` can host any script to create models or statistical analyses with the behavioral features you extract, and `src/visualization/` has scripts to create plots of the raw and processed data.

There are other important files and folders but only relevant if you are interested in extending RAPIDS (e.g. virtual env files, docs, tests, Dockerfile, the Snakefile, etc.). In the figure below, we represent the interactions between users and files. After a user modifies `config.yaml` and `.env` the `Snakefile` file will decide what Snakemake rules have to be executed to produce the required output files (behavioral features) and what scripts are in charge of producing such files. In addition, users can add or modifiy files in the `data` folder (for example to configure the [participants files](https://www.rapids.science/setup/configuration/#participant-files) or the [day segment files](https://www.rapids.science/setup/configuration/#day-segments)).

<figure>
  <img src="/img/files.png" width="600" />
  <figcaption>Interaction diagram between the user, and important files in RAPIDS</figcaption>
</figure>

