# Fitbit Sleep Summary

Sensor parameters description for `[FITBIT_SLEEP_SUMMARY]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Container where your sleep summary data is stored, depending on the data stream you are using this can be a database table, a CSV file, etc. |


## RAPIDS provider

!!! info "Available time segments"
    - Only available for segments that span 1 or more complete days (e.g. Jan 1st 00:00 to Jan 3rd 23:59)

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/fitbit_sleep_summary_raw.csv
    - data/raw/{pid}/fitbit_sleep_summary_with_datetime.csv
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
|countepisodeTYPE               |episodes   |Number of sleep episodes for a certain sleep type during a time segment.
|avgefficiencyTYPE              |scores     |Average sleep efficiency for a certain sleep type during a time segment.
|sumdurationafterwakeupTYPE     |minutes    |Total duration the user stayed in bed after waking up for a certain sleep type during a time segment.
|sumdurationasleepTYPE          |minutes    |Total sleep duration for a certain sleep type during a time segment.
|sumdurationawakeTYPE           |minutes    |Total duration the user stayed awake but still in bed for a certain sleep type during a time segment.
|sumdurationtofallasleepTYPE    |minutes    |Total duration the user spent to fall asleep for a certain sleep type during a time segment.
|sumdurationinbedTYPE           |minutes    |Total duration the user stayed in bed (sumdurationtofallasleep + sumdurationawake + sumdurationasleep + sumdurationafterwakeup) for a certain sleep type during a time segment.
|avgdurationafterwakeupTYPE     |minutes    |Average duration the user stayed in bed after waking up for a certain sleep type during a time segment.
|avgdurationasleepTYPE          |minutes    |Average sleep duration for a certain sleep type during a time segment.
|avgdurationawakeTYPE           |minutes    |Average duration the user stayed awake but still in bed for a certain sleep type during a time segment.
|avgdurationtofallasleepTYPE    |minutes    |Average duration the user spent to fall asleep for a certain sleep type during a time segment.
|avgdurationinbedTYPE           |minutes    |Average duration the user stayed in bed (sumdurationtofallasleep + sumdurationawake + sumdurationasleep + sumdurationafterwakeup) for a certain sleep type during a time segment.



!!! note "Assumptions/Observations"
    
    1. There are three sleep types (TYPE): `main`, `nap`, `all`. The `all` type contains both main sleep and naps.

    2. There are two versions of Fitbit’s sleep API ([version 1](https://dev.fitbit.com/build/reference/web-api/sleep-v1/) and [version 1.2](https://dev.fitbit.com/build/reference/web-api/sleep/)), and each provides raw sleep data in a different format:
        - _Count & duration summaries_. `v1` contains `count_awake`, `duration_awake`, `count_awakenings`, `count_restless`, and `duration_restless` fields for every sleep record but `v1.2` does not.
    
    3. _API columns_. Features are computed based on the values provided by Fitbit’s API: `efficiency`, `minutes_after_wakeup`, `minutes_asleep`, `minutes_awake`, `minutes_to_fall_asleep`, `minutes_in_bed`, `is_main_sleep` and `type`.