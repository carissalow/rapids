# Execution

After you have [installed](../installation) and [configured](../configuration) RAPIDS, use the following command to execute it.

```bash
./rapids -j1
```

!!! info
    The script `#!bash ./rapids` is a wrapper around Snakemake so you can pass any parameters that Snakemake accepts (e.g. `-j1`). 
    
!!! hint "Updating RAPIDS output after modifying `config.yaml`"
    Any changes to the `config.yaml` file will be applied automatically and only the relevant files will be updated. This means that after modifying the features list for `PHONE_MESSAGE` for example, RAPIDS will update the output file with the correct features.

!!! hint "Multi-core"
    You can run RAPIDS over multiple cores by modifying the `-j` argument (e.g. use `-j8` to use 8 cores). **However**, take into account that this means multiple sensor datasets for different participants will be load in memory at the same time. If RAPIDS crashes because it ran out of memory reduce the number of cores and try again.

    As reference, we have run RAPIDS over 12 cores and 32 Gb of RAM without problems for a study with 200 participants with 14 days of low-frequency smartphone data (no accelerometer, gyroscope, or magnetometer).

!!! hint "Forcing a complete rerun"
    If you want to update your data from your database or rerun the whole pipeline from scratch run one or both of the following commands depending on the devices you are using:

    ```bash
    ./rapids -j1 -R download_phone_data
    ./rapids -j1 -R download_fitbit_data
    ```

!!! hint "Deleting RAPIDS output"
    If you  want to delete all the output files RAPIDS produces you can execute the following command (the content of these folders will be deleted: `data/raw`, `data/interim`, `data/processed`, `reports/figures`, and `reports/compliance`)

    ```bash
    ./rapids -j1 -R clean
    ```

!!! done "Ready to extract behavioral features"
    If you are ready to extract features head over to the [Behavioral Features Introduction](../../features/feature-introduction/)