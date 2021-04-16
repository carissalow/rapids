# Fitbit Calories Intraday

Sensor parameters description for `[FITBIT_CALORIES_INTRADAY]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Container where your calories intraday data is stored, depending on the data stream you are using this can be a database table, a CSV file, etc. |


## RAPIDS provider

!!! info "Available time segments"
    - Available for all time segments

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/fitbit_calories_intraday_raw.csv
    - data/raw/{pid}/fitbit_calories_intraday_with_datetime.csv
    - data/interim/{pid}/fitbit_calories_intraday_features/fitbit_calories_intraday_{language}_{provider_key}.csv
    - data/processed/features/{pid}/fitbit_calories_intraday.csv
    ```


Parameters description for `[FITBIT_CALORIES_INTRADAY][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`  | Set to `True` to extract `FITBIT_CALORIES_INTRADAY` features from the `RAPIDS` provider|
|`[FEATURES]` |         Features to be computed from calories intraday data, see table below          |
|`[EPISODE_TYPE]` |    RAPIDS will compute features for any episodes in this list. There are seven types of episodes defined as consecutive appearances of a label. Four are based on the activity level labels provided by Fitbit: `sedentary`, `lightly active`, `fairly active`, and `very active`. One is defined by RAPIDS as moderate to vigorous physical activity `MVPA` episodes that are based on all `fairly active`, and `very active`  labels. Two are defined by the user based on a threshold that divides low or high MET (metabolic equivalent) episodes.        |
|`EPISODE_TIME_THRESHOLD` | Any consecutive rows of the same `[EPISODE_TYPE]` will be considered a single episode if the time difference between them is less or equal than this threshold in minutes|
|`[EPISODE_MET_THRESHOLD]` |    Any 1-minute calorie data chunk with a MET value equal or higher than this threshold will be considered a high MET episode and low MET otherwise.  The default value is 3|
|`[EPISODE_MVPA_CATEGORIES]` |    The Fitbit level labels that are considered part of a moderate to vigorous physical activity episode. One or more of `sedentary`, `lightly active`, `fairly active`, and `very active`. The default are `fairly active` and `very active`|
|`[EPISODE_REFERENCE_TIME]` |   Reference time for the start/end time features. `MIDNIGHT` sets the reference time to 00:00 of each day, `START_OF_THE_SEGMENT` sets the reference time to the start of the time segment (useful when a segment is shorter than a day or spans multiple days)|


Features description for `[FITBIT_CALORIES_INTRADAY][PROVIDERS][RAPIDS]`:

|Feature&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
|starttimefirstepisode`EPISODE_TYPE`               |minutes     |Start time of the first episode of type `[EPISODE_TYPE]`
|endtimefirstepisode`EPISODE_TYPE`               |minutes     |End time of the first episode of type `[EPISODE_TYPE]`
|starttimelastepisode`EPISODE_TYPE`               |minutes     |Start time of the last episode of type `[EPISODE_TYPE]`
|endtimelastepisode`EPISODE_TYPE`               |minutes     |End time of the last episode of type `[EPISODE_TYPE]`
|starttimelongestepisode`EPISODE_TYPE`               |minutes     |Start time of the longest episode of type `[EPISODE_TYPE]`
|endtimelongestepisode`EPISODE_TYPE`               |minutes     |End time of the longest episode of type `[EPISODE_TYPE]`
|countepisode`EPISODE_TYPE`               |episodes     |The number of episodes of type `[EPISODE_TYPE]`
|sumdurationepisode`EPISODE_TYPE`               |minutes     |The sum of the duration of episodes of type `[EPISODE_TYPE]`
|avgdurationepisode`EPISODE_TYPE`               |minutes     |The average of the duration of episodes of type `[EPISODE_TYPE]`
|maxdurationepisode`EPISODE_TYPE`               |minutes     |The maximum of the duration of episodes of type `[EPISODE_TYPE]`
|mindurationepisode`EPISODE_TYPE`               |minutes     |The minimum of the duration of episodes of type `[EPISODE_TYPE]`
|stddurationepisode`EPISODE_TYPE`               |minutes     |The standard deviation of the duration of episodes of type `[EPISODE_TYPE]`
|summet`EPISODE_TYPE`               |METs     |The sum of all METs during episodes of type `[EPISODE_TYPE]`
|avgmet`EPISODE_TYPE`               |METs     |The average of all METs during episodes of type `[EPISODE_TYPE]`
|maxmet`EPISODE_TYPE`               |METs     |The maximum of all METs during episodes of type `[EPISODE_TYPE]`
|minmet`EPISODE_TYPE`               |METs     |The minimum of all METs during episodes of type `[EPISODE_TYPE]`
|stdmet`EPISODE_TYPE`               |METs     |The standard deviation of all METs during episodes of type `[EPISODE_TYPE]`
|sumcalories`EPISODE_TYPE`               |calories     |The sum of all calories during episodes of type `[EPISODE_TYPE]`
|avgcalories`EPISODE_TYPE`               |calories     |The average of all calories during episodes of type `[EPISODE_TYPE]`
|maxcalories`EPISODE_TYPE`               |calories     |The maximum of all calories during episodes of type `[EPISODE_TYPE]`
|mincalories`EPISODE_TYPE`               |calories     |The minimum of all calories during episodes of type `[EPISODE_TYPE]`
|stdcalories`EPISODE_TYPE`               |calories     |The standard deviation of all calories during episodes of type `[EPISODE_TYPE]`


!!! note "Assumptions/Observations"
    - These features are based on intraday calories data that is usually obtained in 1-minute chunks from Fitbit's API.
    - The MET value returned by Fitbit is divided by 10
    - Take into account that the [intraday data returned by Fitbit](https://dev.fitbit.com/build/reference/web-api/activity/#get-activity-intraday-time-series) can contain time series for calories burned inclusive of BMR, tracked activity, and manually logged activities.
