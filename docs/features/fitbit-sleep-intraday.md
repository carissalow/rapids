# Fitbit Sleep Intraday

Sensor parameters description for `[FITBIT_SLEEP_INTRADAY]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[TABLE]`| Database table name or file path where the sleep intraday data is stored. The configuration keys in [Device Data Source Configuration](../../setup/configuration/#device-data-source-configuration) control whether this parameter is interpreted as table or file.
|`[INCLUDE_SLEEP_LATER_THAN]`| All resampled sleep rows (bin interval: one minute) that started after this time will be included in the feature computation. It is a number ranging from 0 (midnight) to 1439 (23:59) which denotes the number of minutes after midnight. If a segment is longer than one day, this value is for every day.
|`[REFERENCE_TIME]`| The reference point from which the `[ROUTINE]` features are to be computed. Chosen from `MIDNIGHT` and `START_OF_THE_SEGMENT`, default is `MIDNIGHT`. If you have multiple time segments per day it might be more informative to set this flag to `START_OF_THE_SEGMENT`.

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

        All columns are mandatory, however, all except `local_date_time` and `duration` can be empty if you don't have that data. Just have in mind that some features might be inaccurate or empty as `type_episode_id`, `level`, `is_main_sleep`, and `type` are used for sleep episodes extraction. `type_episode_id` is based on where it is extracted: if it is extracted from the 1st "minutesData" block, the `type_episode_id` field will be 0. Similarly, the kth block will be k-1.

        |type_episode_id  |local_date_time     |duration  |level      |is_main_sleep  |type           |
        |---------------- |------------------- |--------- |---------- |-------------- |-------------- |
        |0                |2020-10-07 15:55:00 |60        |awake      |0              |classic        |
        |0                |2020-10-07 15:56:00 |60        |awake      |0              |classic        |
        |0                |2020-10-07 15:57:00 |60        |restless   |0              |classic        |

??? example "Example of the structure of source data with Fitbit’s sleep API Version 1.2"

    === "JSON"

        |device_id                                |fitbit_data                                               |
        |---------------------------------------- |--------------------------------------------------------- |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"sleep":[{"dateOfSleep":"2020-10-10","duration":3600000,"efficiency":92,"endTime":"2020-10-10T16:37:00.000","infoCode":2,"isMainSleep":false,"levels":{"data":[{"dateTime":"2020-10-10T15:36:30.000","level":"restless","seconds":60},{"dateTime":"2020-10-10T15:37:30.000","level":"asleep","seconds":660},{"dateTime":"2020-10-10T15:48:30.000","level":"restless","seconds":60},...], "summary":{"asleep":{"count":0,"minutes":56},"awake":{"count":0,"minutes":0},"restless":{"count":3,"minutes":4}}},"logId":26315914306,"minutesAfterWakeup":0,"minutesAsleep":55,"minutesAwake":5,"minutesToFallAsleep":0,"startTime":"2020-10-10T15:36:30.000","timeInBed":60,"type":"classic"},{"dateOfSleep":"2020-10-10","duration":22980000,"efficiency":88,"endTime":"2020-10-10T08:10:00.000","infoCode":0,"isMainSleep":true,"levels":{"data":[{"dateTime":"2020-10-10T01:46:30.000","level":"light","seconds":420},{"dateTime":"2020-10-10T01:53:30.000","level":"deep","seconds":1230},{"dateTime":"2020-10-10T02:14:00.000","level":"light","seconds":360},...], "summary":{"deep":{"count":3,"minutes":92,"thirtyDayAvgMinutes":0},"light":{"count":29,"minutes":193,"thirtyDayAvgMinutes":0},"rem":{"count":4,"minutes":33,"thirtyDayAvgMinutes":0},"wake":{"count":28,"minutes":65,"thirtyDayAvgMinutes":0}}},"logId":26311786557,"minutesAfterWakeup":0,"minutesAsleep":318,"minutesAwake":65,"minutesToFallAsleep":0,"startTime":"2020-10-10T01:46:30.000","timeInBed":383,"type":"stages"}],"summary":{"stages":{"deep":92,"light":193,"rem":33,"wake":65},"totalMinutesAsleep":373,"totalSleepRecords":2,"totalTimeInBed":443}}
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"sleep":[{"dateOfSleep":"2020-10-11","duration":41640000,"efficiency":89,"endTime":"2020-10-11T11:47:00.000","infoCode":0,"isMainSleep":true,"levels":{"data":[{"dateTime":"2020-10-11T00:12:30.000","level":"wake","seconds":450},{"dateTime":"2020-10-11T00:20:00.000","level":"light","seconds":870},{"dateTime":"2020-10-11T00:34:30.000","level":"wake","seconds":780},...], "summary":{"deep":{"count":4,"minutes":52,"thirtyDayAvgMinutes":62},"light":{"count":32,"minutes":442,"thirtyDayAvgMinutes":364},"rem":{"count":6,"minutes":68,"thirtyDayAvgMinutes":58},"wake":{"count":29,"minutes":132,"thirtyDayAvgMinutes":94}}},"logId":26589710670,"minutesAfterWakeup":1,"minutesAsleep":562,"minutesAwake":132,"minutesToFallAsleep":0,"startTime":"2020-10-11T00:12:30.000","timeInBed":694,"type":"stages"}],"summary":{"stages":{"deep":52,"light":442,"rem":68,"wake":132},"totalMinutesAsleep":562,"totalSleepRecords":1,"totalTimeInBed":694}}
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"sleep":[{"dateOfSleep":"2020-10-12","duration":28980000,"efficiency":93,"endTime":"2020-10-12T09:34:30.000","infoCode":0,"isMainSleep":true,"levels":{"data":[{"dateTime":"2020-10-12T01:31:00.000","level":"wake","seconds":600},{"dateTime":"2020-10-12T01:41:00.000","level":"light","seconds":60},{"dateTime":"2020-10-12T01:42:00.000","level":"deep","seconds":2340},...], "summary":{"deep":{"count":4,"minutes":63,"thirtyDayAvgMinutes":59},"light":{"count":27,"minutes":257,"thirtyDayAvgMinutes":364},"rem":{"count":5,"minutes":94,"thirtyDayAvgMinutes":58},"wake":{"count":24,"minutes":69,"thirtyDayAvgMinutes":95}}},"logId":26589710673,"minutesAfterWakeup":0,"minutesAsleep":415,"minutesAwake":68,"minutesToFallAsleep":0,"startTime":"2020-10-12T01:31:00.000","timeInBed":483,"type":"stages"}],"summary":{"stages":{"deep":63,"light":257,"rem":94,"wake":69},"totalMinutesAsleep":415,"totalSleepRecords":1,"totalTimeInBed":483}}
    
    === "PLAIN_TEXT"

        All columns are mandatory, however, all except `local_date_time` and `duration` can be empty if you don't have that data. Just have in mind that some features might be inaccurate or empty as `type_episode_id`, `level`, `is_main_sleep`, and `type` are used for sleep episodes extraction. `type_episode_id` is based on where it is extracted: if it is extracted from the 1st "data" and "shortData" block, the `type_episode_id` field will be 0. Similarly, the kth block will be k-1.

        |type_episode_id  |local_date_time     |duration  |level      |is_main_sleep  |type           |
        |---------------- |------------------- |--------- |---------- |-------------- |-------------- |
        |0                |2020-10-10 15:36:30 |60        |restless   |1              |stages         |
        |0                |2020-10-10 15:37:30 |660       |asleep     |1              |stages         |
        |0                |2020-10-10 15:48:30 |60        |restless   |1              |stages         |

## RAPIDS provider

!!! info "Available time segments"
    - Available for all time segments

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/fitbit_sleep_intraday_raw.csv
    - data/raw/{pid}/fitbit_sleep_intraday_parsed.csv
    - data/interim/{pid}/fitbit_sleep_intraday_episodes_resampled.csv
    - data/interim/{pid}/fitbit_sleep_intraday_episodes_resampled_with_datetime.csv
    - data/interim/{pid}/fitbit_sleep_intraday_features/fitbit_sleep_intraday_{language}_{provider_key}.csv
    - data/processed/features/{pid}/fitbit_sleep_intraday.csv
    ```


Parameters description for `[FITBIT_SLEEP_INTRADAY][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`                                  | Set to `True` to extract `FITBIT_SLEEP_INTRADAY` features from the `RAPIDS` provider|
|`[FEATURES]`                                 |         Features to be computed from sleep intraday data, see table below           |
|`[SLEEP_LEVELS]`                             | Fitbit’s sleep API Version 1 only provides `CLASSIC` records. However, Version 1.2 provides 2 types of records: `CLASSIC` and `STAGES`. `STAGES` is only available in devices with a heart rate sensor and even those devices will fail to report it if the battery is low or the device is not tight enough. While `CLASSIC` contains 3 sleep levels (`awake`, `restless`, and `asleep`), `STAGES` contains 4 sleep levels (`wake`, `deep`, `light`, `rem`). To make it consistent, RAPIDS grouped them into 2 `UNIFIED` sleep levels: `awake` (`CLASSIC`: `awake` and `restless`; `STAGES`: `wake`) and `asleep` (`CLASSIC`: `asleep`; `STAGES`: `deep`, `light`, and `rem`).
|`[SLEEP_TYPES]`                              | Types of sleep to be included in the feature extraction computation. Fitbit provides 2 types of sleep: `main`, `nap`.


Features description for `[FITBIT_STEPS_INTRADAY][PROVIDERS][RAPIDS][LEVELS_AND_TYPES]`:

|Feature&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                    |Units          |Description                                                  |
|------------------------------- |-------------- |-------------------------------------------------------------|
|countepisode`[LEVEL][TYPE]`     |episodes       |Number of `[LEVEL][TYPE]`sleep episodes. `[LEVEL]`is one of `[SLEEP_LEVELS]` (e.g. awake-classic or rem-stages) and `[TYPE]` is one of `[SLEEP_TYPES]` (e.g. main). Both `[LEVEL]`and `[TYPE]` can also be `all` when ``LEVELS_AND_TYPES_COMBINING_ALL`` is True, which ignores the levels and groups by sleep types.
|sumduration`[LEVEL][TYPE]`      |minutes        |Total duration of all `[LEVEL][TYPE]`sleep episodes. `[LEVEL]`is one of `[SLEEP_LEVELS]` (e.g. awake-classic or rem-stages) and `[TYPE]` is one of `[SLEEP_TYPES]` (e.g. main). Both `[LEVEL]` and `[TYPE]`can also be `all` when `LEVELS_AND_TYPES_COMBINING_ALL` is True, which ignores the levels and groups by sleep types.
|maxduration`[LEVEL][TYPE]`      |minutes        | Longest duration of any `[LEVEL][TYPE]`sleep episode. `[LEVEL]`is one of `[SLEEP_LEVELS]` (e.g. awake-classic or rem-stages) and `[TYPE]` is one of `[SLEEP_TYPES]` (e.g. main). Both `[LEVEL]` and `[TYPE]`can also be `all` when `LEVELS_AND_TYPES_COMBINING_ALL` is True, which ignores the levels and groups by sleep types.
|minduration`[LEVEL][TYPE]`      |minutes        | Shortest duration of any `[LEVEL][TYPE]`sleep episode. `[LEVEL]`is one of `[SLEEP_LEVELS]` (e.g. awake-classic or rem-stages) and `[TYPE]` is one of `[SLEEP_TYPES]` (e.g. main). Both `[LEVEL]` and `[TYPE]`can also be `all` when `LEVELS_AND_TYPES_COMBINING_ALL` is True, which ignores the levels and groups by sleep types.
|avgduration`[LEVEL][TYPE]`      |minutes        | Average duration of all `[LEVEL][TYPE]`sleep episodes. `[LEVEL]`is one of `[SLEEP_LEVELS]` (e.g. awake-classic or rem-stages) and `[TYPE]` is one of `[SLEEP_TYPES]` (e.g. main). Both `[LEVEL]` and `[TYPE]`can also be `all` when `LEVELS_AND_TYPES_COMBINING_ALL` is True, which ignores the levels and groups by sleep types.
|medianduration`[LEVEL][TYPE]`   |minutes        | Median duration of all `[LEVEL][TYPE]`sleep episodes. `[LEVEL]`is one of `[SLEEP_LEVELS]` (e.g. awake-classic or rem-stages) and `[TYPE]` is one of `[SLEEP_TYPES]` (e.g. main). Both `[LEVEL]` and `[TYPE]`can also be `all` when `LEVELS_AND_TYPES_COMBINING_ALL` is True, which ignores the levels and groups by sleep types.
|stdduration`[LEVEL][TYPE]`      |minutes        | Standard deviation duration of all `[LEVEL][TYPE]`sleep episodes. `[LEVEL]`is one of `[SLEEP_LEVELS]` (e.g. awake-classic or rem-stages) and `[TYPE]` is one of `[SLEEP_TYPES]` (e.g. main). Both `[LEVEL]` and `[TYPE]`can also be `all` when `LEVELS_AND_TYPES_COMBINING_ALL` is True, which ignores the levels and groups by sleep types.


Features description for `[FITBIT_STEPS_INTRADAY][PROVIDERS][RAPIDS]` RATIOS `[ACROSS_LEVELS]`:

|Feature&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                    |Units |Description                                                        |
|-------------------------- |-------------- |-------------------------------------------------------------|
|ratiocount`[LEVEL]`         |-|Ratio between the **count** of episodes of a single sleep `[LEVEL]` and the **count** of all episodes of all levels during both `main` and `nap` sleep types. This answers the question: what percentage of all `wake`, `deep`, `light`, and `rem` episodes were `rem`? (e.g., $countepisode[remstages][all] / countepisode[all][all]$)
|ratioduration`[LEVEL]`      |-|Ratio between the **duration** of episodes of a single sleep `[LEVEL]` and the **duration** of all episodes of all levels during both `main` and `nap` sleep types. This answers the question: what percentage of all `wake`, `deep`, `light`, and `rem` time was `rem`? (e.g., $sumduration[remstages][all] / sumduration[all][all]$)


Features description for `[FITBIT_STEPS_INTRADAY][PROVIDERS][RAPIDS]` RATIOS `[ACROSS_TYPES]`:

|Feature&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                    |Units          |Description                                                  |
|-------------------------- |-------------- |-------------------------------------------------------------|
|ratiocountmain             |-              |Ratio between the **count** of all `main` episodes (independently of the levels inside) divided by the **count** of all `main` and `nap` episodes. This answers the question: what percentage of all sleep episodes (`main` and `nap`) were `main`? We do not provide the ratio for `nap` because is complementary. ($countepisode[all][main] / countepisode[all][all]$)
|ratiodurationmain          |-              |Ratio between the **duration** of all `main` episodes (independently of the levels inside) divided by the **duration** of all `main` and `nap` episodes. This answers the question: what percentage of all sleep time (`main` and `nap`) was `main`? We do not provide the ratio for `nap` because is complementary. ($sumduration[all][main] / sumduration[all][all]$)


Features description for `[FITBIT_STEPS_INTRADAY][PROVIDERS][RAPIDS]` RATIOS `[WITHIN_LEVELS]`:

|Feature&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                           |Units          |Description                                                  |
|--------------------------------- |-------------- |-------------------------------------------------------------|
|ratiocount`[TYPE]`within`[LEVEL]`    |-              |Ratio between the **count** of episodes of a single sleep `[LEVEL]` during `main` sleep divided by the **count** of episodes of a single sleep `[LEVEL]` during `main` **and** `nap`. This answers the question: are `rem` episodes more frequent during `main` than `nap` sleep? We do not provide the ratio for `nap` because is complementary. ($countepisode[remstages][main] / countepisode[remstages][all]$)
|ratioduration`[TYPE]`within`[LEVEL]` |-              |Ratio between the **duration** of episodes of a single sleep `[LEVEL]` during `main` sleep divided by the **duration** of episodes of a single sleep `[LEVEL]` during `main` **and** `nap`. This answers the question: is `rem` time more frequent during `main` than `nap` sleep? We do not provide the ratio for `nap` because is complementary. ($countepisode[remstages][main] / countepisode[remstages][all]$)


Features description for `[FITBIT_STEPS_INTRADAY][PROVIDERS][RAPIDS]` RATIOS `[WITHIN_TYPES]`:

|Feature&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|Units|Description|
| - |- | - |
|ratiocount`[LEVEL]`within`[TYPE]`    |-|Ratio between the **count** of episodes of a single sleep `[LEVEL]` and the **count** of all episodes of all levels during either `main` or `nap` sleep types. This answers the question: what percentage of all `wake`, `deep`, `light`, and `rem` episodes were `rem` during `main`/`nap` sleep time? (e.g., $countepisode[remstages][main] / countepisode[all][main]$)
|ratioduration`[LEVEL]`within`[TYPE]` |-|Ratio between the **duration** of episodes of a single sleep `[LEVEL]` and the **duration** of all episodes of all levels during either `main` or `nap` sleep types. This answers the question: what percentage of all `wake`, `deep`, `light`, and `rem` time was `rem` during `main`/`nap` sleep time? (e.g., $sumduration[remstages][main] / sumduration[all][main]$)


Features description for `[FITBIT_STEPS_INTRADAY][PROVIDERS][RAPIDS][ROUTINE]`:

|Feature                           |Units          |Description                                                  |
|--------------------------------- |-------------- |-------------------------------------------------------------|
|starttimefirstmainsleep           |minutes        |Start time (in minutes since `REFERENCE_TIME`) of the first main sleep episode after `INCLUDE_EPISODES_LATER_THAN`.
|endtimelastmainsleep              |minutes        |End time (in minutes since `REFERENCE_TIME`) of the last main sleep episode after `INCLUDE_EPISODES_LATER_THAN`.
|starttimefirstnap                 |minutes        |Start time (in minutes since `REFERENCE_TIME`) of the first nap episode after `INCLUDE_EPISODES_LATER_THAN`.
|endtimelastnap                    |minutes        |End time (in minutes since `REFERENCE_TIME`) of the last nap episode after `INCLUDE_EPISODES_LATER_THAN`.





!!! note "Assumptions/Observations"
    1. Deleting values from `[SLEEP_LEVELS]` or `[SLEEP_TYPES]` will only change the features you receive from `[LEVELS_AND_TYPES]`. For example if `STAGES` only contains `[rem, light]` you will not receive `countepisode[wake|deep][TYPE]` or sum, max, min, avg, median, or std `duration`. These values will not influence `RATIOS` or `ROUTINE` features.
    2. Any `[LEVEL]` grouping is done within the elements of each class `CLASSIC`, `STAGES`, and `UNIFIED`. That is, we never combine `CLASSIC` or `STAGES` types to compute features when `LEVELS_AND_TYPES_COMBINING_ALL` is True or when computing `RATIOS`.
    

## PRICE provider

!!! info "Available time segments"
    - Available for any time segments larger or equal to one day

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/fitbit_sleep_intraday_raw.csv
    - data/raw/{pid}/fitbit_sleep_intraday_parsed.csv
    - data/interim/{pid}/fitbit_sleep_intraday_episodes_resampled.csv
    - data/interim/{pid}/fitbit_sleep_intraday_episodes_resampled_with_datetime.csv
    - data/interim/{pid}/fitbit_sleep_intraday_features/fitbit_sleep_intraday_{language}_{provider_key}.csv
    - data/processed/features/{pid}/fitbit_sleep_intraday.csv
    ```


Parameters description for `[FITBIT_SLEEP_INTRADAY][PROVIDERS][PRICE]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`                                  | Set to `True` to extract `FITBIT_SLEEP_INTRADAY` features from the `PRICE` provider                                      |
|`[FEATURES]`                                 |         Features to be computed from sleep intraday data, see table below   
|`[DAY_TYPE]`                                             | The features of this provider can be computed using daily averages/medians that were extracted on `WEEKEND` days only, `WEEK` days only, or `ALL` days|
|`[INCLUDE_EPISODES_INTERSECTING]`            | This parameter contains 2 values: `[START_TIME]` and `[LENGTH]`. Only `main` sleep episodes that ended later than `[START_TIME]` or start before `[START_TIME]` + `[LENGTH]` are taken into account to compute the features described below. Both `[START_TIME]` and `[LENGTH]` are in minutes. `[START_TIME]` is a number ranging from 0 (midnight) to 1439 (23:59) which denotes the number of minutes after midnight. `[LENGTH]` is a number smaller than 1440 (24 hours).  |


Features description for `[FITBIT_STEPS_INTRADAY][PROVIDERS][PRICE]`:

|Feature&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                                  |Units          |Description                                                  |
|------------------------------------- |----------------- |-------------------------------------------------------------|
|avgratioawaketoinbedmain`[DAY_TYPE]`  |-                 | Average daily ratio between daily awake time and in-bed time inferred from `main` sleep episodes. An sleep episode can have `CLASSIC` or `STAGES` levels, awake time is formed by the `[awake, restless]` `CLASSIC` levels and the `[wake]` `STAGES` level. In turn, in-bed time is the total duration of all `main` sleep episodes for each day. You can include daily ratios that were computed on weekend days, week days or both depending on the value of the `DAY_TYPE` flag.
|avgstarttimeofepisodemain`[DAY_TYPE]` |minutes           | Average start time of the first `main` sleep episode of each day in a time segment. You can include daily start times from episodes detected on weekend days, week days or both depending on the value of the `DAY_TYPE` flag.
|avgendtimeofepisodemain`[DAY_TYPE]`   |minutes           |  Average end time of the last `main` sleep episode of each day in a time segment. You can include daily end times from episodes detected on weekend days, week days or both depending on the value of the `DAY_TYPE` flag.
|avgmidpointofepisodemain`[DAY_TYPE]`  |minutes           |  Average mid time between the start of the first `main` sleep episode and the end of the last `main` sleep episode of each day in a time segment. You can include episodes detected on weekend days, week days or both depending on the value of the `DAY_TYPE` flag.
|socialjetlag                          |minutes           | Difference in minutes between the avgstarttimeofepisodemain (bed time) of weekends and weekdays. 
|meanssdstarttimeofepisodemain         |minutes squared   | Same as `avgstarttimeofepisodemain[DAY_TYPE]` but the average is computed over the squared differences of each pair of consecutive start times.
|meanssdendtimeofepisodemain           |minutes squared   | Same as `avgendtimeofepisodemain[DAY_TYPE]` but the average is computed over the squared differences of each pair of consecutive end times.
|meanssdmidpointofepisodemain          |minutes squared   | Same as `avgmidpointofepisodemain[DAY_TYPE]` but the average is computed over the squared differences of each pair of consecutive mid times.
|medianssdstarttimeofepisodemain       |minutes squared   | Same as `avgstarttimeofepisodemain[DAY_TYPE]` but the median is computed over the squared differences of each pair of consecutive start times.
|medianssdendtimeofepisodemain         |minutes squared   | Same as `avgendtimeofepisodemain[DAY_TYPE]` but the median is computed over the squared differences of each pair of consecutive end times.
|medianssdmidpointofepisodemain        |minutes squared   | Same as `avgmidpointofepisodemain[DAY_TYPE]` but the median is computed over the squared differences of each pair of consecutive mid times.



!!! note "Assumptions/Observations"
    1. These features are based on descriptive statistics computed across daily values (start/end/mid times of sleep episodes). This is the reason why they are only available on time segments that are longer than 24 hours (we need at least 1 day to get the average).
    1. Even though Fitbit provides 2 types of sleep episodes (`main` and `nap`), only `main` sleep episodes are considered.
    2. We sum 24 to any start or bed times that happen *after* midnight so they can be averaged with other values happening *before* midnight.
    2. `main` sleep episodes belong to the day they start at and are only included in the feature computation if they start or end (overlap) between `[START_TIME]` and `[START_TIME]` + `[LENGTH]`. For example: 
        
        If we process the following `main` sleep episodes: 

        | episode |start|end|
        |-|-|-|
        |1|2021-02-01 12:00|2021-02-01 15:00|
        |2|2021-02-01 21:00|2021-02-02 03:00|
        |3|2021-02-02 05:00|2021-02-02 08:00|
        |4|2021-02-02 11:00|2021-02-02 14:00|
        |5|2021-02-02 19:00|2021-02-03 06:00|

        And our parameters:

        - `[INCLUDE_EPISODES_INTERSECTING][START_TIME]` = 1320 (today's 22:00)

        - `[INCLUDE_EPISODES_INTERSECTING][LENGTH]` = 720 (tomorrow's 10:00, or 22:00 + 12 hours)
        
        Only sleep episodes 2, 3,and 5 would be considered.

    3. All `main` sleep episodes are chunked within the requested [time segments](../../setup/configuration/#time-segments) which need to be at least 24 hours or more long (1, 2, 3, 7 days, etc.). Then, daily features will be extracted and averaged across the length of the time segment, for example:
    
        The daily features extracted on 2021-02-01 will be:

        - starttimeofepisodemain (bedtime) is `21 * 60` (episode 2 start time 2021-02-01 21:00)

        - endtimeofepisodemain (wake time) is `32 * 60 `(episode 3 end time 2021-02-02 08:00 + 24)

        - midpointofepisodemain (midpoint sleep) is `[(32 * 60) - (21 * 60)] / 2`


        The daily features extracted on 2021-02-02 will be:

        - starttimeofepisodemain (bedtime) is `19 * 60` (episode 5 start time 2021-02-01 19:00)

        - endtimeofepisodemain (wake time) is `30 * 60 `(episode 5 end time 2021-02-03 06:00 + 24)

        - midpointofepisodemain (midpoint sleep) is `[(30 * 60) - (19 * 60)] / 2`

        And `avgstarttimeofepisodemain[DAY_TYPE]` will be `([21 * 60] + [19 * 60]) / 2` 


    



