# Execution

After you have [installed](../installation) and [configured](../configuration) RAPIDS, use the following command to execute it.

```bash
./rapids -j1
```

!!! done "Ready to extract behavioral features"
    If you are ready to extract features head over to the [Behavioral Features Introduction](../../features/feature-introduction/)

!!! hint "We wrap Snakemake" 
    The script `#!bash ./rapids` is a wrapper around Snakemake so you can pass any parameters that Snakemake accepts (e.g. `-j1`). 
    
!!! hint "Updating RAPIDS output after modifying `config.yaml`"
    Any changes to the `config.yaml` file will be applied automatically and only the relevant files will be updated. This means that after modifying the features list for `PHONE_MESSAGE` for example, RAPIDS will execute the script that computes `MESSAGES` features and update its output file.

!!! hint "Multi-core"
    You can run RAPIDS over multiple cores by modifying the `-j` argument (e.g. use `-j8` to use 8 cores). **However**, take into account that this means multiple sensor datasets for different participants will be loaded in memory at the same time. If RAPIDS crashes because it ran out of memory, reduce the number of cores and try again.

    As reference, we have run RAPIDS over 12 cores and 32 Gb of RAM without problems for a study with 200 participants with 14 days of low-frequency smartphone data (no accelerometer, gyroscope, or magnetometer).

!!! hint "Deleting RAPIDS output"
    If you  want to delete all the output files RAPIDS produces, you can execute the following command:

    ```bash
    ./rapids -j1 --delete-all-output
    ```

!!! hint "Forcing a complete rerun or updating your raw data in RAPIDS"
    If you want to update your raw data or rerun the whole pipeline from scratch, run the following commands:

    ```bash
    ./rapids -j1 --delete-all-output
    ./rapids -j1
    ```