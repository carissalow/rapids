# Fitbit Steps Summary

Sensor parameters description for `[FITBIT_STEPS_SUMMARY]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[TABLE]`| Database table name or file path where the steps summary data is stored. The configuration keys in [Device Data Source Configuration](../../setup/configuration/#device-data-source-configuration) control whether this parameter is interpreted as table or file.

The format of the column(s) containing the Fitbit sensor data can be `JSON` or `PLAIN_TEXT`. The data in `JSON` format is obtained directly from the Fitbit API. We support `PLAIN_TEXT` in case you already parsed your data and don't have access to your participants' Fitbit accounts anymore. If your data is in `JSON` format then summary and intraday data come packed together. 

We provide examples of the input format that RAPIDS expects, note that both examples for `JSON` and `PLAIN_TEXT` are tabular and the actual format difference comes in the `fitbit_data` column (we truncate the `JSON` example for brevity).

??? example "Example of the structure of source data"

    === "JSON"

        |device_id                                |fitbit_data                                               |
        |---------------------------------------- |--------------------------------------------------------- |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |"activities-steps":[{"dateTime":"2020-10-07","value":"1775"}],"activities-steps-intraday":{"dataset":[{"time":"00:00:00","value":5},{"time":"00:01:00","value":3},{"time":"00:02:00","value":0},...],"datasetInterval":1,"datasetType":"minute"}}
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |"activities-steps":[{"dateTime":"2020-10-08","value":"3201"}],"activities-steps-intraday":{"dataset":[{"time":"00:00:00","value":14},{"time":"00:01:00","value":11},{"time":"00:02:00","value":10},...],"datasetInterval":1,"datasetType":"minute"}}
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |"activities-steps":[{"dateTime":"2020-10-09","value":"998"}],"activities-steps-intraday":{"dataset":[{"time":"00:00:00","value":0},{"time":"00:01:00","value":0},{"time":"00:02:00","value":0},...],"datasetInterval":1,"datasetType":"minute"}}
    
    === "PLAIN_TEXT"
        All columns are mandatory.

        |device_id                              |local_date_time        |steps     |
        |-------------------------------------- |---------------------- |--------- |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-07             |1775      |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-08             |3201      |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-09             |998       |


## RAPIDS provider

!!! info "Available time segments"
    - Only available for segments that span 1 or more complete days (e.g. Jan 1st 00:00 to Jan 3rd 23:59)

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/fitbit_steps_summary_raw.csv
    - data/raw/{pid}/fitbit_steps_summary_parsed.csv
    - data/raw/{pid}/fitbit_steps_summary_parsed_with_datetime.csv
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
