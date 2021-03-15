# Fitbit Steps Summary

Sensor parameters description for `[FITBIT_STEPS_SUMMARY]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Container where your steps summary data is stored, depending on the data stream you are using this can be a database table, a CSV file, etc. |


## RAPIDS provider

!!! info "Available time segments"
    - Only available for segments that span 1 or more complete days (e.g. Jan 1st 00:00 to Jan 3rd 23:59)

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/fitbit_steps_summary_raw.csv
    - data/raw/{pid}/fitbit_steps_summary_with_datetime.csv
    - data/interim/{pid}/fitbit_steps_summary_features/fitbit_steps_summary_{language}_{provider_key}.csv
    - data/processed/features/{pid}/fitbit_steps_summary.csv
    ```


Parameters description for `[FITBIT_STEPS_SUMMARY][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`  | Set to `True` to extract `FITBIT_STEPS_SUMMARY` features from the `RAPIDS` provider|
|`[FEATURES]` |         Features to be computed from steps summary data, see table below          |


Features description for `[FITBIT_STEPS_SUMMARY][PROVIDERS][RAPIDS]`:

|Feature                    |Units      |Description                                  |
|-------------------------- |---------- |-------------------------------------------- |
|maxsumsteps                |steps      |The maximum daily step count during a time segment.
|minsumsteps                |steps      |The minimum daily step count during a time segment.
|avgsumsteps                |steps      |The average daily step count during a time segment.
|mediansumsteps             |steps      |The median of daily step count during a time segment.
|stdsumsteps                |steps      |The standard deviation of daily step count during a time segment.

!!! note "Assumptions/Observations"
    
    NA
