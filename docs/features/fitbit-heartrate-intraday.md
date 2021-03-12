# Fitbit Heart Rate Intraday

Sensor parameters description for `[FITBIT_HEARTRATE_INTRADAY]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Container where your heart rate intraday data is stored, depending on the data stream you are using this can be a database table, a CSV file, etc. |


## RAPIDS provider

!!! info "Available time segments"
    - Available for all time segments

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/fitbit_heartrate_intraday_raw.csv
    - data/raw/{pid}/fitbit_heartrate_intraday_with_datetime.csv
    - data/interim/{pid}/fitbit_heartrate_intraday_features/fitbit_heartrate_intraday_{language}_{provider_key}.csv
    - data/processed/features/{pid}/fitbit_heartrate_intraday.csv
    ```


Parameters description for `[FITBIT_HEARTRATE_INTRADAY][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`  | Set to `True` to extract `FITBIT_HEARTRATE_INTRADAY` features from the `RAPIDS` provider|
|`[FEATURES]` |         Features to be computed from heart rate intraday data, see table below          |


Features description for `[FITBIT_HEARTRATE_INTRADAY][PROVIDERS][RAPIDS]`:

|Feature                    |Units          |Description|
|-------------------------- |-------------- |---------------------------|
|maxhr                      |beats/mins     |The maximum heart rate during a time segment.
|minhr                      |beats/mins     |The minimum heart rate during a time segment.
|avghr                      |beats/mins     |The average heart rate during a time segment.
|medianhr                   |beats/mins     |The median of heart rate during a time segment.
|modehr                     |beats/mins     |The mode of heart rate during a time segment.
|stdhr                      |beats/mins     |The standard deviation of heart rate during a time segment.
|diffmaxmodehr              |beats/mins     |The difference between the maximum and mode heart rate during a time segment.
|diffminmodehr              |beats/mins     |The difference between the mode and minimum heart rate during a time segment.
|entropyhr                  |nats           |Shannon’s entropy measurement based on heart rate during a time segment.
|minutesonZONE              |minutes        |Number of minutes the user’s heart rate fell within each `heartrate_zone` during a time segment.

!!! note "Assumptions/Observations"
    
    1. There are four heart rate zones (ZONE): ``outofrange``, ``fatburn``, ``cardio``, and ``peak``. Please refer to [Fitbit documentation](https://help.fitbit.com/articles/en_US/Help_article/1565.htm) for more information about the way they are computed.
