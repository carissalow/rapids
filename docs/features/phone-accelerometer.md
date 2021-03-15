# Phone Accelerometer

Sensor parameters description for `[PHONE_ACCELEROMETER]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Data stream [container](../../datastreams/data-streams-introduction/) (database table, CSV file, etc.) where the accelerometer data is stored

## RAPIDS provider

!!! info "Available time segments and platforms"
    - Available for all time segments
    - Available for Android and iOS

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/phone_accelerometer_raw.csv
    - data/raw/{pid}/phone_accelerometer_with_datetime.csv
    - data/interim/{pid}/phone_accelerometer_features/phone_accelerometer_{language}_{provider_key}.csv
    - data/processed/features/{pid}/phone_accelerometer.csv
    ```


Parameters description for `[PHONE_ACCELEROMETER][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`| Set to `True` to extract `PHONE_ACCELEROMETER` features from the `RAPIDS` provider|
|`[FEATURES]` |         Features to be computed, see table below


Features description for `[PHONE_ACCELEROMETER][PROVIDERS][RAPIDS]`:

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
|maxmagnitude      |m/s^2^    |The maximum magnitude of acceleration ($\|acceleration\| = \sqrt{x^2 + y^2 + z^2}$).
|minmagnitude      |m/s^2^    |The minimum magnitude of acceleration.
|avgmagnitude      |m/s^2^    |The average magnitude of acceleration.
|medianmagnitude   |m/s^2^    |The median magnitude of acceleration.
|stdmagnitude      |m/s^2^    |The standard deviation of acceleration.

!!! note "Assumptions/Observations"
    1. Analyzing accelerometer data is a memory intensive task. If RAPIDS crashes is likely because the accelerometer dataset for a participant is to big to fit in memory. We are considering different alternatives to overcome this problem.

## PANDA provider

These features are based on the work by [Panda et al](../../citation#panda-accelerometer).

!!! info "Available time segments and platforms"
    - Available for all time segments
    - Available for Android and iOS

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/phone_accelerometer_raw.csv
    - data/raw/{pid}/phone_accelerometer_with_datetime.csv
    - data/interim/{pid}/phone_accelerometer_features/phone_accelerometer_{language}_{provider_key}.csv
    - data/processed/features/{pid}/phone_accelerometer.csv
    ```


Parameters description for `[PHONE_ACCELEROMETER][PROVIDERS][PANDA]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`| Set to `True` to extract `PHONE_ACCELEROMETER` features from the `PANDA` provider|
|`[FEATURES]` |         Features to be computed for exertional and non-exertional activity episodes, see table below


Features description for `[PHONE_ACCELEROMETER][PROVIDERS][PANDA]`:

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
| sumduration    | minutes | Total duration of all exertional or non-exertional activity episodes.                     |
| maxduration    | minutes | Longest duration of any exertional or non-exertional activity episode.                    |
| minduration    | minutes | Shortest duration of any exertional or non-exertional activity episode.                   |
| avgduration    | minutes | Average duration of any exertional or non-exertional activity episode.                    |
| medianduration | minutes | Median duration of any exertional or non-exertional activity episode.                     |
| stdduration    | minutes | Standard deviation of the duration of all exertional or non-exertional activity episodes. |

!!! note "Assumptions/Observations"
    1. Analyzing accelerometer data is a memory intensive task. If RAPIDS crashes is likely because the accelerometer dataset for a participant is to big to fit in memory. We are considering different alternatives to overcome this problem.
    2. See [Panda et al](../../citation#panda-accelerometer) for a definition of exertional and non-exertional activity episodes
