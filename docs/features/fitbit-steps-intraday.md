# Fitbit Steps Intraday

Sensor parameters description for `[FITBIT_STEPS_INTRADAY]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[TABLE]`| Database table name or file path where the steps intraday data is stored. The configuration keys in [Device Data Source Configuration](../../setup/configuration/#device-data-source-configuration) control whether this parameter is interpreted as table or file.

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

        |device_id                              |local_date_time        |steps     |
        |-------------------------------------- |---------------------- |--------- |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-07 00:00:00    |5         |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-07 00:01:00    |3         |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-07 00:02:00    |0         |


## RAPIDS provider

!!! info "Available time segments"
    - Available for all time segments

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/fitbit_steps_intraday_raw.csv
    - data/raw/{pid}/fitbit_steps_intraday_parsed.csv
    - data/raw/{pid}/fitbit_steps_intraday_parsed_with_datetime.csv
    - data/interim/{pid}/fitbit_steps_intraday_features/fitbit_steps_intraday_{language}_{provider_key}.csv
    - data/processed/features/{pid}/fitbit_steps_intraday.csv
    ```


Parameters description for `[FITBIT_STEPS_INTRADAY][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`                | Set to `True` to extract `FITBIT_STEPS_INTRADAY` features from the `RAPIDS` provider|
|`[FEATURES]`               |         Features to be computed from steps intraday data, see table below           |
|`[THRESHOLD_ACTIVE_BOUT]`  | Every minute with Fitbit steps data wil be labelled as `sedentary` if its step count is below this threshold, otherwise, `active`.    |
|`[INCLUDE_ZERO_STEP_ROWS]` | Whether or not to include time segments with a 0 step count during the whole day.                          |


Features description for `[FITBIT_STEPS_INTRADAY][PROVIDERS][RAPIDS]`:

|Feature                    |Units          |Description                                                  |
|-------------------------- |-------------- |-------------------------------------------------------------|
|sumsteps                   |steps          |The total step count during a time segment.
|maxsteps                   |steps          |The maximum step count during a time segment.
|minsteps                   |steps          |The minimum step count during a time segment.
|avgsteps                   |steps          |The average step count during a time segment.
|stdsteps                   |steps          |The standard deviation of step count during a time segment.
|countepisodesedentarybout  |bouts          |Number of sedentary bouts during a time segment.
|sumdurationsedentarybout   |minutes        |Total duration of all sedentary bouts during a time segment.
|maxdurationsedentarybout   |minutes        |The maximum duration of any sedentary bout during a time segment.
|mindurationsedentarybout   |minutes        |The minimum duration of any sedentary bout during a time segment.
|avgdurationsedentarybout   |minutes        |The average duration of sedentary bouts during a time segment.
|stddurationsedentarybout   |minutes        |The standard deviation of the duration of sedentary bouts during a time segment.
|countepisodeactivebout     |bouts          |Number of active bouts during a time segment.
|sumdurationactivebout      |minutes        |Total duration of all active bouts during a time segment.
|maxdurationactivebout      |minutes        |The maximum duration of any active bout during a time segment.
|mindurationactivebout      |minutes        |The minimum duration of any active bout during a time segment.
|avgdurationactivebout      |minutes        |The average duration of active bouts during a time segment.
|stddurationactivebout      |minutes        |The standard deviation of the duration of active bouts during a time segment.

!!! note "Assumptions/Observations"
    
    1. _Active and sedentary bouts_. If the step count per minute is smaller than `THRESHOLD_ACTIVE_BOUT` (default value is 10), that minute is labelled as sedentary, otherwise, is labelled as active. Active and sedentary bouts are periods of consecutive minutes labelled as `active` or `sedentary`.

