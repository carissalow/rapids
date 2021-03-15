# Fitbit Heart Rate Summary

Sensor parameters description for `[FITBIT_HEARTRATE_SUMMARY]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Container where your heart rate summary data is stored, depending on the data stream you are using this can be a database table, a CSV file, etc. |


## RAPIDS provider

!!! info "Available time segments"
    - Only available for segments that span 1 or more complete days (e.g. Jan 1st 00:00 to Jan 3rd 23:59)

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/fitbit_heartrate_summary_raw.csv
    - data/raw/{pid}/fitbit_heartrate_summary_with_datetime.csv
    - data/interim/{pid}/fitbit_heartrate_summary_features/fitbit_heartrate_summary_{language}_{provider_key}.csv
    - data/processed/features/{pid}/fitbit_heartrate_summary.csv
    ```


Parameters description for `[FITBIT_HEARTRATE_SUMMARY][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`  | Set to `True` to extract `FITBIT_HEARTRATE_SUMMARY` features from the `RAPIDS` provider|
|`[FEATURES]` |         Features to be computed from heart rate summary data, see table below          |


Features description for `[FITBIT_HEARTRATE_SUMMARY][PROVIDERS][RAPIDS]`:

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
|maxrestinghr               |beats/mins     |The maximum daily resting heart rate during a time segment.
|minrestinghr               |beats/mins     |The minimum daily resting heart rate during a time segment.
|avgrestinghr               |beats/mins     |The average daily resting heart rate during a time segment.
|medianrestinghr            |beats/mins     |The median of daily resting heart rate during a time segment.
|moderestinghr              |beats/mins     |The mode of daily resting heart rate during a time segment.
|stdrestinghr               |beats/mins     |The standard deviation of daily resting heart rate during a time segment.
|diffmaxmoderestinghr       |beats/mins     |The difference between the maximum and mode daily resting heart rate during a time segment.
|diffminmoderestinghr       |beats/mins     |The difference between the mode and minimum daily resting heart rate during a time segment.
|entropyrestinghr           |nats           |Shannon’s entropy measurement based on daily resting heart rate during a time segment.
|sumcaloriesZONE            |cals           |The total daily calories burned within `heartrate_zone` during a time segment.
|maxcaloriesZONE            |cals           |The maximum daily calories burned within `heartrate_zone` during a time segment.
|mincaloriesZONE            |cals           |The minimum daily calories burned within `heartrate_zone` during a time segment.
|avgcaloriesZONE            |cals           |The average daily calories burned within `heartrate_zone` during a time segment.
|mediancaloriesZONE         |cals           |The median of daily calories burned within `heartrate_zone` during a time segment.
|stdcaloriesZONE            |cals           |The standard deviation of daily calories burned within `heartrate_zone` during a time segment.
|entropycaloriesZONE        |nats           |Shannon’s entropy measurement based on daily calories burned within `heartrate_zone` during a time segment.

!!! note "Assumptions/Observations"
    
    1. There are four heart rate zones (ZONE): ``outofrange``, ``fatburn``, ``cardio``, and ``peak``. Please refer to [Fitbit documentation](https://help.fitbit.com/articles/en_US/Help_article/1565.htm) for more information about the way they are computed.

    2. Calories' accuracy depends on the users’ Fitbit profile (weight, height, etc.).
