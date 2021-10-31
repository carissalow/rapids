# Fitbit Steps Intraday

Sensor parameters description for `[FITBIT_STEPS_INTRADAY]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Container where your steps intraday data is stored, depending on the data stream you are using this can be a database table, a CSV file, etc. |
|`[EXCLUDE_SLEEP]` | Step data will be excluded if it was logged during sleep periods when at least one `[EXCLUDE]` flag is set to `True`. Sleep can be delimited by (1) a fixed period that repeats on every day if `[TIME_BASED][EXCLUDE]` is True or (2) by Fitbit summary sleep episodes if `[FITBIT_BASED][EXCLUDE]` is True. If both are True (3), we use all Fitbit sleep episodes as well as the time-based episodes that do not overlap with any Fitbit episodes. If `[TIME_BASED][EXCLUDE]` is True, make sure Fitbit sleep summary container points to a valid table or file.

## RAPIDS provider

!!! info "Available time segments"
    - Available for all time segments

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/fitbit_steps_intraday_raw.csv
    - data/raw/{pid}/fitbit_steps_intraday_with_datetime.csv
    - data/raw/{pid}/fitbit_sleep_summary_raw.csv (Only when [EXCLUDE_SLEEP][EXCLUDE]=True and [EXCLUDE_SLEEP][TYPE]=FITBIT_BASED)
    - data/interim/{pid}/fitbit_steps_intraday_with_datetime_exclude_sleep.csv (Only when [EXCLUDE_SLEEP][EXCLUDE]=True)
    - data/interim/{pid}/fitbit_steps_intraday_features/fitbit_steps_intraday_{language}_{provider_key}.csv
    - data/processed/features/{pid}/fitbit_steps_intraday.csv
    ```


Parameters description for `[FITBIT_STEPS_INTRADAY][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`                | Set to `True` to extract `FITBIT_STEPS_INTRADAY` features from the `RAPIDS` provider|
|`[FEATURES]`               |         Features to be computed from steps intraday data, see table below           |
|`[REFERENCE_HOUR]`         | The reference point from which `firststeptime` or `laststeptime` is to be computed, default is midnight |
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
|firststeptime              |minutes        |Minutes until the first non-zero step count.
|laststeptime               |minutes        |Minutes until the last non-zero step count.
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

