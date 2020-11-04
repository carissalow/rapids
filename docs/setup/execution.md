# Execution

After you have [installed](/setup/installation) and [configured](/setup/configuration) RAPIDS, use the following command to execute it.

```bash
./rapids -j1
```

!!! info
    The script `#!bash ./rapids` is a wrapper around Snakemake so you can pass any parameters that Snakemake accepts (e.g. `-j1`). 
    
    Any changes to the `config.yaml` file will be applied automatically and only the relevant files will be updated.

!!! hint "Multi-core"
    You can run RAPIDS over multiple cores by modifying the `-j` argument (e.g. use `-j8` to use 8 cores). **However**, take into account that this means multiple sensor datasets for different participants will be load in memory at the same time. If RAPIDS crashes because it ran out of memory reduce the number of cores and try again.

    As reference, we have run RAPIDS over 12 cores and 32 Gb of RAM without problems for a study with 200 participants with 14 days of low-frequency smartphone data (no accelerometer, gyroscope, or magnetometer).