# Phone Activity Recognition

Sensor parameters description for `[PHONE_ACTIVITY_RECOGNITION]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER][ANDROID]`| Data stream [container](../../datastreams/data-streams-introduction/) (database table, CSV file, etc.) where the activity data from Android devices is stored (the AWARE client saves this data on different tables for Android and iOS)
|`[CONTAINER][IOS]`| Data stream [container](../../datastreams/data-streams-introduction/) (database table, CSV file, etc.) where the activity data from iOS devices is stored (the AWARE client saves this data on different tables for Android and iOS)
|`[EPISODE_THRESHOLD_BETWEEN_ROWS]` | Difference in minutes between any two rows for them to be considered part of the same activity episode

## RAPIDS provider

!!! info "Available time segments and platforms"
    - Available for all time segments
    - Available for Android and iOS

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/phone_activity_recognition_raw.csv
    - data/raw/{pid}/phone_activity_recognition_with_datetime.csv
    - data/interim/{pid}/phone_activity_recognition_episodes.csv
    - data/interim/{pid}/phone_activity_recognition_episodes_resampled.csv
    - data/interim/{pid}/phone_activity_recognition_episodes_resampled_with_datetime.csv
    - data/interim/{pid}/phone_activity_recognition_features/phone_activity_recognition_{language}_{provider_key}.csv
    - data/processed/features/{pid}/phone_activity_recognition.csv
    ```


Parameters description for `[PHONE_ACTIVITY_RECOGNITION][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`| Set to `True` to extract `PHONE_ACTIVITY_RECOGNITION` features from the `RAPIDS` provider|
|`[FEATURES]` |         Features to be computed, see table below
|`[ACTIVITY_CLASSES][STATIONARY]` | An array of the activity labels to be considered in the `STATIONARY` category choose any of `still`, `tilting`
|`[ACTIVITY_CLASSES][MOBILE]` | An array of the activity labels to be considered in the `MOBILE` category choose any of `on_foot`, `walking`, `running`, `on_bicycle`
|`[ACTIVITY_CLASSES][VEHICLE]` | An array of the activity labels to be considered in the `VEHICLE` category choose any of `in_vehicule`


Features description for `[PHONE_ACTIVITY_RECOGNITION][PROVIDERS][RAPIDS]`:

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
|count                   |rows             | Number of episodes.
|mostcommonactivity      |activity type   | The most common activity type (e.g. `still`, `on_foot`, etc.). If there is a tie, the first one is chosen.
|countuniqueactivities   |activity type   | Number of unique activities.
|durationstationary      |minutes          | The total duration of `[ACTIVITY_CLASSES][STATIONARY]` episodes of still and tilting activities
|durationmobile          |minutes          | The total duration of `[ACTIVITY_CLASSES][MOBILE]` episodes of on foot, running, and on bicycle activities
|durationvehicle         |minutes          | The total duration of `[ACTIVITY_CLASSES][VEHICLE]` episodes of on vehicle activity

!!! note "Assumptions/Observations"
    1. iOS Activity Recognition names and types are unified with Android labels: 

        | iOS Activity Name | Android Activity Name | Android Activity Type |
        |----|----|----|
        |`walking`|  `walking` |  `7`
        |`running`|  `running` |  `8`
        |`cycling`|  `on_bicycle` |  `1`
        |`automotive`|  `in_vehicle` |  `0`
        |`stationary`|  `still` |  `3`
        |`unknown`|  `unknown` |  `4`

    2. In AWARE, Activity Recognition data for Android and iOS are stored in two different database tables, RAPIDS automatically infers what platform each participant belongs to based on their [participant file](../../setup/configuration/#participant-files).