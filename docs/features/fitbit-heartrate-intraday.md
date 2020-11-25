# Fitbit Heart Rate Intraday

Sensor parameters description for `[FITBIT_HEARTRATE_INTRADAY]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[TABLE]`| Database table name or file path where the heart rate intraday data is stored. The configuration keys in [Device Data Source Configuration](../../setup/configuration/#device-data-source-configuration) control whether this parameter is interpreted as table or file.

The format of the column(s) containing the Fitbit sensor data can be `JSON` or `PLAIN_TEXT`. The data in `JSON` format is obtained directly from the Fitbit API. We support `PLAIN_TEXT` in case you already parsed your data and don't have access to your participants' Fitbit accounts anymore. If your data is in `JSON` format then summary and intraday data come packed together. 

We provide examples of the input format that RAPIDS expects, note that both examples for `JSON` and `PLAIN_TEXT` are tabular and the actual format difference comes in the `fitbit_data` column (we truncate the `JSON` example for brevity).

??? example "Example of the structure of source data"

    === "JSON"

        |device_id                                |fitbit_data                                               |
        |---------------------------------------- |--------------------------------------------------------- |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"activities-heart":[{"dateTime":"2020-10-07","value":{"customHeartRateZones":[],"heartRateZones":[{"caloriesOut":1200.6102,"max":88,"min":31,"minutes":1058,"name":"Out of Range"},{"caloriesOut":760.3020,"max":120,"min":86,"minutes":366,"name":"Fat Burn"},{"caloriesOut":15.2048,"max":146,"min":120,"minutes":2,"name":"Cardio"},{"caloriesOut":0,"max":221,"min":148,"minutes":0,"name":"Peak"}],"restingHeartRate":72}}],"activities-heart-intraday":{"dataset":[{"time":"00:00:00","value":68},{"time":"00:01:00","value":67},{"time":"00:02:00","value":67},...],"datasetInterval":1,"datasetType":"minute"}}
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"activities-heart":[{"dateTime":"2020-10-08","value":{"customHeartRateZones":[],"heartRateZones":[{"caloriesOut":1100.1120,"max":89,"min":30,"minutes":921,"name":"Out of Range"},{"caloriesOut":660.0012,"max":118,"min":82,"minutes":361,"name":"Fat Burn"},{"caloriesOut":23.7088,"max":142,"min":108,"minutes":3,"name":"Cardio"},{"caloriesOut":0,"max":221,"min":148,"minutes":0,"name":"Peak"}],"restingHeartRate":70}}],"activities-heart-intraday":{"dataset":[{"time":"00:00:00","value":77},{"time":"00:01:00","value":75},{"time":"00:02:00","value":73},...],"datasetInterval":1,"datasetType":"minute"}}
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"activities-heart":[{"dateTime":"2020-10-09","value":{"customHeartRateZones":[],"heartRateZones":[{"caloriesOut":750.3615,"max":77,"min":30,"minutes":851,"name":"Out of Range"},{"caloriesOut":734.1516,"max":107,"min":77,"minutes":550,"name":"Fat Burn"},{"caloriesOut":131.8579,"max":130,"min":107,"minutes":29,"name":"Cardio"},{"caloriesOut":0,"max":220,"min":130,"minutes":0,"name":"Peak"}],"restingHeartRate":69}}],"activities-heart-intraday":{"dataset":[{"time":"00:00:00","value":90},{"time":"00:01:00","value":89},{"time":"00:02:00","value":88},...],"datasetInterval":1,"datasetType":"minute"}}
    
    === "PLAIN_TEXT"

        |device_id                              |local_date_time        |heartrate |heartrate_zone  |
        |-------------------------------------- |---------------------- |--------- |--------------- |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-07 00:00:00    |68        |outofrange      |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-07 00:01:00    |67        |outofrange      |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-07 00:02:00    |67        |outofrange      |


## RAPIDS provider

!!! info "Available day segments"
    - Available for all day segments

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/fitbit_heartrate_intraday_raw.csv
    - data/raw/{pid}/fitbit_heartrate_intraday_parsed.csv
    - data/raw/{pid}/fitbit_heartrate_intraday_parsed_with_datetime.csv
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
|maxhr                      |beats/mins     |The maximum heart rate during a day segment.
|minhr                      |beats/mins     |The minimum heart rate during a day segment.
|avghr                      |beats/mins     |The average heart rate during a day segment.
|medianhr                   |beats/mins     |The median of heart rate during a day segment.
|modehr                     |beats/mins     |The mode of heart rate during a day segment.
|stdhr                      |beats/mins     |The standard deviation of heart rate during a day segment.
|diffmaxmodehr              |beats/mins     |The difference between the maximum and mode heart rate during a day segment.
|diffminmodehr              |beats/mins     |The difference between the mode and minimum heart rate during a day segment.
|entropyhr                  |nats           |Shannon’s entropy measurement based on heart rate during a day segment.
|minutesonZONE              |minutes        |Number of minutes the user’s heart rate fell within each `heartrate_zone` during a day segment.

!!! note "Assumptions/Observations"
    
    1. There are four heart rate zones (ZONE): ``outofrange``, ``fatburn``, ``cardio``, and ``peak``. Please refer to [Fitbit documentation](https://help.fitbit.com/articles/en_US/Help_article/1565.htm) for more information about the way they are computed.
