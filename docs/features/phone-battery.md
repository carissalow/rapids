# Phone Battery

Sensor parameters description for `[PHONE_BATTERY]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Data stream [container](../../datastreams/data-streams-introduction/) (database table, CSV file, etc.) where the battery data is stored
|`[EPISODE_THRESHOLD_BETWEEN_ROWS]` | Difference in minutes between any two rows for them to be considered part of the same battery charge or discharge episode

## RAPIDS provider

!!! info "Available time segments and platforms"
    - Available for all time segments
    - Available for Android and iOS

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/phone_battery_raw.csv
    - data/interim/{pid}/phone_battery_episodes.csv
    - data/interim/{pid}/phone_battery_episodes_resampled.csv
    - data/interim/{pid}/phone_battery_episodes_resampled_with_datetime.csv
    - data/interim/{pid}/phone_battery_features/phone_battery_{language}_{provider_key}.csv
    - data/processed/features/{pid}/phone_battery.csv
    ```


Parameters description for `[PHONE_BATTERY][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`| Set to `True` to extract `PHONE_BATTERY` features from the `RAPIDS` provider|
|`[FEATURES]` |         Features to be computed, see table below


Features description for `[PHONE_BATTERY][PROVIDERS][RAPIDS]`:

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
|countdischarge         |episodes           | Number of discharging episodes.
|sumdurationdischarge   |minutes            | The total duration of all discharging episodes.
|countcharge            |episodes           | Number of battery charging episodes.
|sumdurationcharge      |minutes            | The total duration of all charging episodes.
|avgconsumptionrate     |episodes/minutes   | The average of all episodes' consumption rates. An episode's consumption rate is defined as the ratio between its battery delta and duration
|maxconsumptionrate     |episodes/minutes   | The highest of all episodes' consumption rates. An episode's consumption rate is defined as the ratio between its battery delta and duration

!!! note "Assumptions/Observations"
    1. We convert battery data collected with iOS client v1 (autodetected because battery status `4` do not exist) to match Android battery format: we swap status `3` for `5` and `1` for `3`
    2. We group battery data into discharge or charge episodes considering any contiguous rows with consecutive reductions or increases of the battery level if they are logged within `[EPISODE_THRESHOLD_BETWEEN_ROWS]` minutes from each other.
