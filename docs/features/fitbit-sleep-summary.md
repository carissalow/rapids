# Fitbit Sleep Summary

Sensor parameters description for `[FITBIT_SLEEP_SUMMARY]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Container where your sleep summary data is stored, depending on the data stream you are using this can be a database table, a CSV file, etc. |


## RAPIDS provider

!!! hint "Understanding RAPIDS features"
    [This diagram](../../img/sleep_summary_rapids.png) will help you understand how sleep episodes are chunked and grouped within time segments using `SLEEP_SUMMARY_LAST_NIGHT_END` for the RAPIDS provider.

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
|`[SLEEP_TYPES]` | Types of sleep to be included in the feature extraction computation. There are three sleep types: `main`, `nap`, and `all`. The `all` type means both main sleep and naps are considered.       |
|`[FEATURES]`    |         Features to be computed from sleep summary data, see table below                                                           |
|`[FITBIT_DATA_STREAMS][data stream][SLEEP_SUMMARY_LAST_NIGHT_END]`    |  As an exception, the `LAST_NIGHT_END` parameter for this provider is in the data stream configuration section. This parameter controls how sleep episodes are assigned to different days and affects wake and bedtimes.|


Features description for `[FITBIT_SLEEP_SUMMARY][PROVIDERS][RAPIDS]`:

|Feature                        |Units      |Description                                  |
|------------------------------ |---------- |-------------------------------------------- |
|firstwaketimeTYPE              |minutes    |First wake time for a certain sleep type during a time segment. Wake time is number of minutes after midnight of a sleep episode's end time.
|lastwaketimeTYPE               |minutes    |Last wake time for a certain sleep type during a time segment. Wake time is number of minutes after midnight of a sleep episode's end time.
|firstbedtimeTYPE                   |minutes    |First bedtime for a certain sleep type during a time segment. Bedtime is number of minutes after midnight of a sleep episode's start time.
|lastbedtimeTYPE                    |minutes    |Last bedtime for a certain sleep type during a time segment. Bedtime is number of minutes after midnight of a sleep episode's start time.
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
    1. [This diagram](../../img/sleep_summary_rapids.png) will help you understand how sleep episodes are chunked and grouped within time segments using `LNE` for the RAPIDS provider.
    1. There are three sleep types (TYPE): `main`, `nap`, `all`. The `all` type groups both `main` sleep and `naps`. All types are based on Fitbit's labels.
    2. There are two versions of Fitbit’s sleep API ([version 1](https://dev.fitbit.com/build/reference/web-api/sleep-v1/) and [version 1.2](https://dev.fitbit.com/build/reference/web-api/sleep/)), and each provides raw sleep data in a different format:
        - _Count & duration summaries_. `v1` contains `count_awake`, `duration_awake`, `count_awakenings`, `count_restless`, and `duration_restless` fields for every sleep record but `v1.2` does not.
    3. _API columns_. Most features are computed based on the values provided by Fitbit’s API: `efficiency`, `minutes_after_wakeup`, `minutes_asleep`, `minutes_awake`, `minutes_to_fall_asleep`, `minutes_in_bed`, `is_main_sleep` and `type`.
    4. Bed time and sleep duration are based on episodes that started between today’s LNE and tomorrow’s LNE while awake time is based on the episodes that started between yesterday’s LNE and today’s LNE
    5. The reference point for bed/awake times is today’s 00:00. You can have bedtimes larger than 24 and awake times smaller than 0
    6. These features are only available for time segments that span midnight to midnight of the same or different day.
    7. We include first and last wake and bedtimes because, when `LAST_NIGHT_END` is 10 am, the first bedtime could match a nap at 2 pm, and the last bedtime could match a main overnight sleep episode that starts at 10pm.
    5. Set the value for `SLEEP_SUMMARY_LAST_NIGHT_END` int the config parameter [FITBIT_DATA_STREAMS][data stream][SLEEP_SUMMARY_LAST_NIGHT_END].