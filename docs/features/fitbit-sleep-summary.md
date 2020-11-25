# Fitbit Sleep Summary

Sensor parameters description for `[FITBIT_SLEEP_SUMMARY]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[TABLE]`| Database table name or file path where the sleep summary data is stored. The configuration keys in [Device Data Source Configuration](../../setup/configuration/#device-data-source-configuration) control whether this parameter is interpreted as table or file.

The format of the column(s) containing the Fitbit sensor data can be `JSON` or `PLAIN_TEXT`. The data in `JSON` format is obtained directly from the Fitbit API. We support `PLAIN_TEXT` in case you already parsed your data and don't have access to your participants' Fitbit accounts anymore. If your data is in `JSON` format then summary and intraday data come packed together. 

We provide examples of the input format that RAPIDS expects, note that both examples for `JSON` and `PLAIN_TEXT` are tabular and the actual format difference comes in the `fitbit_data` column (we truncate the `JSON` example for brevity).

??? example "Example of the structure of source data with Fitbit’s sleep API Version 1"

    === "JSON"

        |device_id                                |fitbit_data                                               |
        |---------------------------------------- |--------------------------------------------------------- |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"sleep": [{"awakeCount": 2, "awakeDuration": 3, "awakeningsCount": 10, "dateOfSleep": "2020-10-07", "duration": 8100000, "efficiency": 91, "endTime": "2020-10-07T18:10:00.000", "isMainSleep": true, "logId": 14147921940, "minuteData": [{"dateTime": "15:55:00", "value": "3"}, {"dateTime": "15:56:00", "value": "3"}, {"dateTime": "15:57:00", "value": "2"},...], "minutesAfterWakeup": 0, "minutesAsleep": 123, "minutesAwake": 12, "minutesToFallAsleep": 0, "restlessCount": 8, "restlessDuration": 9, "startTime": "2020-10-07T15:55:00.000", "timeInBed": 135}, {"awakeCount": 0, "awakeDuration": 0, "awakeningsCount": 1, "dateOfSleep": "2020-10-07", "duration": 3780000, "efficiency": 100, "endTime": "2020-10-07T10:52:30.000", "isMainSleep": false, "logId": 14144903977, "minuteData": [{"dateTime": "09:49:00", "value": "1"}, {"dateTime": "09:50:00", "value": "1"}, {"dateTime": "09:51:00", "value": "1"},...], "minutesAfterWakeup": 1, "minutesAsleep": 62, "minutesAwake": 0, "minutesToFallAsleep": 0, "restlessCount": 1, "restlessDuration": 1, "startTime": "2020-10-07T09:49:00.000", "timeInBed": 63}], "summary": {"totalMinutesAsleep": 185, "totalSleepRecords": 2, "totalTimeInBed": 198}}
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"sleep": [{"awakeCount": 3, "awakeDuration": 21, "awakeningsCount": 16, "dateOfSleep": "2020-10-08", "duration": 19260000, "efficiency": 89, "endTime": "2020-10-08T06:01:30.000", "isMainSleep": true, "logId": 14150613895, "minuteData": [{"dateTime": "00:40:00", "value": "3"}, {"dateTime": "00:41:00", "value": "3"}, {"dateTime": "00:42:00", "value": "3"},...], "minutesAfterWakeup": 0, "minutesAsleep": 275, "minutesAwake": 33, "minutesToFallAsleep": 0, "restlessCount": 13, "restlessDuration": 25, "startTime": "2020-10-08T00:40:00.000", "timeInBed": 321}], "summary": {"totalMinutesAsleep": 275, "totalSleepRecords": 1, "totalTimeInBed": 321}}
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"sleep": [{"awakeCount": 1, "awakeDuration": 3, "awakeningsCount": 8, "dateOfSleep": "2020-10-09", "duration": 19320000, "efficiency": 96, "endTime": "2020-10-09T05:57:30.000", "isMainSleep": true, "logId": 14161136803, "minuteData": [{"dateTime": "00:35:30", "value": "2"}, {"dateTime": "00:36:30", "value": "1"}, {"dateTime": "00:37:30", "value": "1"},...], "minutesAfterWakeup": 0, "minutesAsleep": 309, "minutesAwake": 13, "minutesToFallAsleep": 0, "restlessCount": 7, "restlessDuration": 10, "startTime": "2020-10-09T00:35:30.000", "timeInBed": 322}], "summary": {"totalMinutesAsleep": 309, "totalSleepRecords": 1, "totalTimeInBed": 322}}
    
    === "PLAIN_TEXT"

        |device_id                              |local_start_date_time  |local_end_date_time    |efficiency  |minutes_after_wakeup  |minutes_asleep  |minutes_awake  |minutes_to_fall_asleep  |minutes_in_bed  |is_main_sleep  |type     |count_awake |duration_awake  |count_awakenings  |count_restless  |duration_restless  |
        |-------------------------------------- |---------------------- |---------------------- |----------- |--------------------- |--------------- |-------------- |----------------------- |--------------- |-------------- |-------- |----------- |--------------- |----------------- |--------------- |------------------ |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-07 15:55:00    |2020-10-07 18:10:00    |91          |0                     |123             |12             |0                       |135             |1              |classic  |2           |3               |10                |8               |9                  |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-07 09:49:00    |2020-10-07 10:52:30    |100         |1                     |62              |0              |0                       |63              |0              |classic  |0           |0               |1                 |1               |1                  |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-08 00:40:00    |2020-10-08 06:01:30    |89          |0                     |275             |33             |0                       |321             |1              |classic  |3           |21              |16                |13              |25                 |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-09 00:35:30    |2020-10-09 05:57:30    |96          |0                     |309             |13             |0                       |322             |1              |classic  |1           |3               |8                 |7               |10                 |

??? example "Example of the structure of source data with Fitbit’s sleep API Version 1.2"

    === "JSON"

        |device_id                                |fitbit_data                                               |
        |---------------------------------------- |--------------------------------------------------------- |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"sleep":[{"dateOfSleep":"2020-10-10","duration":3600000,"efficiency":92,"endTime":"2020-10-10T16:37:00.000","infoCode":2,"isMainSleep":false,"levels":{"data":[{"dateTime":"2020-10-10T15:36:30.000","level":"restless","seconds":60},{"dateTime":"2020-10-10T15:37:30.000","level":"asleep","seconds":660},{"dateTime":"2020-10-10T15:48:30.000","level":"restless","seconds":60},...], "summary":{"asleep":{"count":0,"minutes":56},"awake":{"count":0,"minutes":0},"restless":{"count":3,"minutes":4}}},"logId":26315914306,"minutesAfterWakeup":0,"minutesAsleep":55,"minutesAwake":5,"minutesToFallAsleep":0,"startTime":"2020-10-10T15:36:30.000","timeInBed":60,"type":"classic"},{"dateOfSleep":"2020-10-10","duration":22980000,"efficiency":88,"endTime":"2020-10-10T08:10:00.000","infoCode":0,"isMainSleep":true,"levels":{"data":[{"dateTime":"2020-10-10T01:46:30.000","level":"light","seconds":420},{"dateTime":"2020-10-10T01:53:30.000","level":"deep","seconds":1230},{"dateTime":"2020-10-10T02:14:00.000","level":"light","seconds":360},...], "summary":{"deep":{"count":3,"minutes":92,"thirtyDayAvgMinutes":0},"light":{"count":29,"minutes":193,"thirtyDayAvgMinutes":0},"rem":{"count":4,"minutes":33,"thirtyDayAvgMinutes":0},"wake":{"count":28,"minutes":65,"thirtyDayAvgMinutes":0}}},"logId":26311786557,"minutesAfterWakeup":0,"minutesAsleep":318,"minutesAwake":65,"minutesToFallAsleep":0,"startTime":"2020-10-10T01:46:30.000","timeInBed":383,"type":"stages"}],"summary":{"stages":{"deep":92,"light":193,"rem":33,"wake":65},"totalMinutesAsleep":373,"totalSleepRecords":2,"totalTimeInBed":443}}
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"sleep":[{"dateOfSleep":"2020-10-11","duration":41640000,"efficiency":89,"endTime":"2020-10-11T11:47:00.000","infoCode":0,"isMainSleep":true,"levels":{"data":[{"dateTime":"2020-10-11T00:12:30.000","level":"wake","seconds":450},{"dateTime":"2020-10-11T00:20:00.000","level":"light","seconds":870},{"dateTime":"2020-10-11T00:34:30.000","level":"wake","seconds":780},...], "summary":{"deep":{"count":4,"minutes":52,"thirtyDayAvgMinutes":62},"light":{"count":32,"minutes":442,"thirtyDayAvgMinutes":364},"rem":{"count":6,"minutes":68,"thirtyDayAvgMinutes":58},"wake":{"count":29,"minutes":132,"thirtyDayAvgMinutes":94}}},"logId":26589710670,"minutesAfterWakeup":1,"minutesAsleep":562,"minutesAwake":132,"minutesToFallAsleep":0,"startTime":"2020-10-11T00:12:30.000","timeInBed":694,"type":"stages"}],"summary":{"stages":{"deep":52,"light":442,"rem":68,"wake":132},"totalMinutesAsleep":562,"totalSleepRecords":1,"totalTimeInBed":694}}
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"sleep":[{"dateOfSleep":"2020-10-12","duration":28980000,"efficiency":93,"endTime":"2020-10-12T09:34:30.000","infoCode":0,"isMainSleep":true,"levels":{"data":[{"dateTime":"2020-10-12T01:31:00.000","level":"wake","seconds":600},{"dateTime":"2020-10-12T01:41:00.000","level":"light","seconds":60},{"dateTime":"2020-10-12T01:42:00.000","level":"deep","seconds":2340},...], "summary":{"deep":{"count":4,"minutes":63,"thirtyDayAvgMinutes":59},"light":{"count":27,"minutes":257,"thirtyDayAvgMinutes":364},"rem":{"count":5,"minutes":94,"thirtyDayAvgMinutes":58},"wake":{"count":24,"minutes":69,"thirtyDayAvgMinutes":95}}},"logId":26589710673,"minutesAfterWakeup":0,"minutesAsleep":415,"minutesAwake":68,"minutesToFallAsleep":0,"startTime":"2020-10-12T01:31:00.000","timeInBed":483,"type":"stages"}],"summary":{"stages":{"deep":63,"light":257,"rem":94,"wake":69},"totalMinutesAsleep":415,"totalSleepRecords":1,"totalTimeInBed":483}}
    
    === "PLAIN_TEXT"

        |device_id                              |local_start_date_time  |local_end_date_time    |efficiency  |minutes_after_wakeup  |minutes_asleep  |minutes_awake  |minutes_to_fall_asleep  |minutes_in_bed  |is_main_sleep  |type     |
        |-------------------------------------- |---------------------- |---------------------- |----------- |--------------------- |--------------- |-------------- |----------------------- |--------------- |-------------- |-------- |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-10 15:36:30    |2020-10-10 16:37:00    |92          |0                     |55              |5              |0                       |60              |0              |classic  |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-10 01:46:30    |2020-10-10 08:10:00    |88          |0                     |318             |65             |0                       |383             |1              |stages   |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-11 00:12:30    |2020-10-11 11:47:00    |89          |1                     |562             |132            |0                       |694             |1              |stages   |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-12 01:31:00    |2020-10-12 09:34:30    |93          |0                     |415             |68             |0                       |483             |1              |stages   |


## RAPIDS provider

!!! info "Available day segments"
    - Only available for segments that span 1 or more complete days (e.g. Jan 1st 00:00 to Jan 3rd 23:59)

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/fitbit_sleep_summary_raw.csv
    - data/raw/{pid}/fitbit_sleep_summary_parsed.csv
    - data/raw/{pid}/fitbit_sleep_summary_parsed_with_datetime.csv
    - data/interim/{pid}/fitbit_sleep_summary_features/fitbit_sleep_summary_{language}_{provider_key}.csv
    - data/processed/features/{pid}/fitbit_sleep_summary.csv
    ```


Parameters description for `[FITBIT_SLEEP_SUMMARY][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`     | Set to `True` to extract `FITBIT_SLEEP_SUMMARY` features from the `RAPIDS` provider                                                |
|`[SLEEP_TYPES]` | Types of sleep to be included in the feature extraction computation. Fitbit provides 3 types of sleep: `main`, `nap`, `all`.       |
|`[FEATURES]`    |         Features to be computed from sleep summary data, see table below                                                           |


Features description for `[FITBIT_SLEEP_SUMMARY][PROVIDERS][RAPIDS]`:

|Feature                        |Units      |Description                                  |
|------------------------------ |---------- |-------------------------------------------- |
|countepisodeTYPE               |episodes   |Number of sleep episodes for a certain sleep type during a day segment.
|avgefficiencyTYPE              |scores     |Average sleep efficiency for a certain sleep type during a day segment.
|sumdurationafterwakeupTYPE     |minutes    |Total duration the user stayed in bed after waking up for a certain sleep type during a day segment.
|sumdurationasleepTYPE          |minutes    |Total sleep duration for a certain sleep type during a day segment.
|sumdurationawakeTYPE           |minutes    |Total duration the user stayed awake but still in bed for a certain sleep type during a day segment.
|sumdurationtofallasleepTYPE    |minutes    |Total duration the user spent to fall asleep for a certain sleep type during a day segment.
|sumdurationinbedTYPE           |minutes    |Total duration the user stayed in bed (sumdurationtofallasleep + sumdurationawake + sumdurationasleep + sumdurationafterwakeup) for a certain sleep type during a day segment.
|avgdurationafterwakeupTYPE     |minutes    |Average duration the user stayed in bed after waking up for a certain sleep type during a day segment.
|avgdurationasleepTYPE          |minutes    |Average sleep duration for a certain sleep type during a day segment.
|avgdurationawakeTYPE           |minutes    |Average duration the user stayed awake but still in bed for a certain sleep type during a day segment.
|avgdurationtofallasleepTYPE    |minutes    |Average duration the user spent to fall asleep for a certain sleep type during a day segment.
|avgdurationinbedTYPE           |minutes    |Average duration the user stayed in bed (sumdurationtofallasleep + sumdurationawake + sumdurationasleep + sumdurationafterwakeup) for a certain sleep type during a day segment.



!!! note "Assumptions/Observations"
    
    1. There are three sleep types (TYPE): `main`, `nap`, `all`. The `all` type contains both main sleep and naps.

    2. There are two versions of Fitbit’s sleep API ([version 1](https://dev.fitbit.com/build/reference/web-api/sleep-v1/) and [version 1.2](https://dev.fitbit.com/build/reference/web-api/sleep/)), and each provides raw sleep data in a different format:
        - _Count & duration summaries_. `v1` contains `count_awake`, `duration_awake`, `count_awakenings`, `count_restless`, and `duration_restless` fields for every sleep record but `v1.2` does not.
    
    3. _API columns_. Features are computed based on the values provided by Fitbit’s API: `efficiency`, `minutes_after_wakeup`, `minutes_asleep`, `minutes_awake`, `minutes_to_fall_asleep`, `minutes_in_bed`, `is_main_sleep` and `type`.