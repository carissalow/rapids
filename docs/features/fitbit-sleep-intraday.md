# Fitbit Sleep Intraday

Sensor parameters description for `[FITBIT_SLEEP_INTRADAY]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Container where your sleep intraday data is stored, depending on the data stream you are using this can be a database table, a CSV file, etc. |

## RAPIDS provider

!!! hint "Understanding RAPIDS features"
    [This diagram](../../img/sleep_intraday_rapids.png) will help you understand how sleep episodes are chunked and grouped within time segments for the RAPIDS provider.


!!! info "Available time segments"
    - Available for all time segments

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/fitbit_sleep_intraday_raw.csv
    - data/raw/{pid}/fitbit_sleep_intraday_with_datetime.csv
    - data/interim/{pid}/fitbit_sleep_intraday_episodes.csv
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
|`[SLEEP_LEVELS]`                             | Fitbit’s sleep API Version 1 only provides `CLASSIC` records. However, Version 1.2 provides 2 types of records: `CLASSIC` and `STAGES`. `STAGES` is only available in devices with a heart rate sensor and even those devices will fail to report it if the battery is low or the device is not tight enough. While `CLASSIC` contains 3 sleep levels (`awake`, `restless`, and `asleep`), `STAGES` contains 4 sleep levels (`wake`, `deep`, `light`, `rem`). To make it consistent, RAPIDS groups them into 2 `UNIFIED` sleep levels: `awake` (`CLASSIC`: `awake` and `restless`; `STAGES`: `wake`) and `asleep` (`CLASSIC`: `asleep`; `STAGES`: `deep`, `light`, and `rem`). In this section, there is a boolean flag named `INCLUDE_ALL_GROUPS` that if set to TRUE, computes LEVELS_AND_TYPES features grouping all levels together in a single `all` category.
|`[SLEEP_TYPES]`                              | Types of sleep to be included in the feature extraction computation. There are three sleep types: `main`, `nap`, and `all`. The `all` type means both main sleep and naps are considered.


Features description for `[FITBIT_SLEEP_INTRADAY][PROVIDERS][RAPIDS][LEVELS_AND_TYPES]`:

|Feature&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                    |Units          |Description                                                  |
|------------------------------- |-------------- |-------------------------------------------------------------|
|countepisode`[LEVEL][TYPE]`     |episodes       |Number of `[LEVEL][TYPE]`sleep episodes. `[LEVEL]`is one of `[SLEEP_LEVELS]` (e.g. awake-classic or rem-stages) and `[TYPE]` is one of `[SLEEP_TYPES]` (e.g. main). `[LEVEL]` can also be `all` when `INCLUDE_ALL_GROUPS` is True, which ignores the levels and groups by sleep types.
|sumduration`[LEVEL][TYPE]`      |minutes        |Total duration of all `[LEVEL][TYPE]`sleep episodes. `[LEVEL]`is one of `[SLEEP_LEVELS]` (e.g. awake-classic or rem-stages) and `[TYPE]` is one of `[SLEEP_TYPES]` (e.g. main). `[LEVEL]` can also be `all` when `INCLUDE_ALL_GROUPS` is True, which ignores the levels and groups by sleep types.
|maxduration`[LEVEL][TYPE]`      |minutes        | Longest duration of any `[LEVEL][TYPE]`sleep episode. `[LEVEL]`is one of `[SLEEP_LEVELS]` (e.g. awake-classic or rem-stages) and `[TYPE]` is one of `[SLEEP_TYPES]` (e.g. main). `[LEVEL]` can also be `all` when `INCLUDE_ALL_GROUPS` is True, which ignores the levels and groups by sleep types.
|minduration`[LEVEL][TYPE]`      |minutes        | Shortest duration of any `[LEVEL][TYPE]`sleep episode. `[LEVEL]`is one of `[SLEEP_LEVELS]` (e.g. awake-classic or rem-stages) and `[TYPE]` is one of `[SLEEP_TYPES]` (e.g. main). `[LEVEL]` can also be `all` when `INCLUDE_ALL_GROUPS` is True, which ignores the levels and groups by sleep types.
|avgduration`[LEVEL][TYPE]`      |minutes        | Average duration of all `[LEVEL][TYPE]`sleep episodes. `[LEVEL]`is one of `[SLEEP_LEVELS]` (e.g. awake-classic or rem-stages) and `[TYPE]` is one of `[SLEEP_TYPES]` (e.g. main). `[LEVEL]` can also be `all` when `INCLUDE_ALL_GROUPS` is True, which ignores the levels and groups by sleep types.
|medianduration`[LEVEL][TYPE]`   |minutes        | Median duration of all `[LEVEL][TYPE]`sleep episodes. `[LEVEL]`is one of `[SLEEP_LEVELS]` (e.g. awake-classic or rem-stages) and `[TYPE]` is one of `[SLEEP_TYPES]` (e.g. main). `[LEVEL]` can also be `all` when `INCLUDE_ALL_GROUPS` is True, which ignores the levels and groups by sleep types.
|stdduration`[LEVEL][TYPE]`      |minutes        | Standard deviation duration of all `[LEVEL][TYPE]`sleep episodes. `[LEVEL]`is one of `[SLEEP_LEVELS]` (e.g. awake-classic or rem-stages) and `[TYPE]` is one of `[SLEEP_TYPES]` (e.g. main). `[LEVEL]` can also be `all` when `INCLUDE_ALL_GROUPS` is True, which ignores the levels and groups by sleep types.


Features description for `[FITBIT_SLEEP_INTRADAY][PROVIDERS][RAPIDS]` RATIOS `[ACROSS_LEVELS]`:

|Feature&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                    |Units |Description                                                        |
|-------------------------- |-------------- |-------------------------------------------------------------|
|ratiocount`[LEVEL]`         |-|Ratio between the **count** of episodes of a single sleep `[LEVEL]` and the **count** of all episodes of all levels during both `main` and `nap` sleep types. This answers the question: what percentage of all `wake`, `deep`, `light`, and `rem` episodes were `rem`? (e.g., $countepisode[remstages][all] / countepisode[all][all]$)
|ratioduration`[LEVEL]`      |-|Ratio between the **duration** of episodes of a single sleep `[LEVEL]` and the **duration** of all episodes of all levels during both `main` and `nap` sleep types. This answers the question: what percentage of all `wake`, `deep`, `light`, and `rem` time was `rem`? (e.g., $sumduration[remstages][all] / sumduration[all][all]$)


Features description for `[FITBIT_SLEEP_INTRADAY][PROVIDERS][RAPIDS]` RATIOS `[ACROSS_TYPES]`:

|Feature&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                    |Units          |Description                                                  |
|-------------------------- |-------------- |-------------------------------------------------------------|
|ratiocountmain             |-              |Ratio between the **count** of all `main` episodes (independently of the levels inside) divided by the **count** of all `main` and `nap` episodes. This answers the question: what percentage of all sleep episodes (`main` and `nap`) were `main`? We do not provide the ratio for `nap` because is complementary. ($countepisode[all][main] / countepisode[all][all]$)
|ratiodurationmain          |-              |Ratio between the **duration** of all `main` episodes (independently of the levels inside) divided by the **duration** of all `main` and `nap` episodes. This answers the question: what percentage of all sleep time (`main` and `nap`) was `main`? We do not provide the ratio for `nap` because is complementary. ($sumduration[all][main] / sumduration[all][all]$)


Features description for `[FITBIT_SLEEP_INTRADAY][PROVIDERS][RAPIDS]` RATIOS `[WITHIN_LEVELS]`:

|Feature&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                           |Units          |Description                                                  |
|--------------------------------- |-------------- |-------------------------------------------------------------|
|ratiocountmainwithin`[LEVEL]`    |-              |Ratio between the **count** of episodes of a single sleep `[LEVEL]` during `main` sleep divided by the **count** of episodes of a single sleep `[LEVEL]` during `main` **and** `nap`. This answers the question: are `rem` episodes more frequent during `main` than `nap` sleep? We do not provide the ratio for `nap` because is complementary. ($countepisode[remstages][main] / countepisode[remstages][all]$)
|ratiodurationmainwithin`[LEVEL]` |-              |Ratio between the **duration** of episodes of a single sleep `[LEVEL]` during `main` sleep divided by the **duration** of episodes of a single sleep `[LEVEL]` during `main` **and** `nap`. This answers the question: is `rem` time more frequent during `main` than `nap` sleep? We do not provide the ratio for `nap` because is complementary. ($countepisode[remstages][main] / countepisode[remstages][all]$)


Features description for `[FITBIT_SLEEP_INTRADAY][PROVIDERS][RAPIDS]` RATIOS `[WITHIN_TYPES]`:

|Feature&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|Units|Description|
| - |- | - |
|ratiocount`[LEVEL]`within`[TYPE]`    |-|Ratio between the **count** of episodes of a single sleep `[LEVEL]` and the **count** of all episodes of all levels during either `main` or `nap` sleep types. This answers the question: what percentage of all `wake`, `deep`, `light`, and `rem` episodes were `rem` during `main`/`nap` sleep time? (e.g., $countepisode[remstages][main] / countepisode[all][main]$)
|ratioduration`[LEVEL]`within`[TYPE]` |-|Ratio between the **duration** of episodes of a single sleep `[LEVEL]` and the **duration** of all episodes of all levels during either `main` or `nap` sleep types. This answers the question: what percentage of all `wake`, `deep`, `light`, and `rem` time was `rem` during `main`/`nap` sleep time? (e.g., $sumduration[remstages][main] / sumduration[all][main]$)



!!! note "Assumptions/Observations"
    1. [This diagram](../../img/sleep_intraday_rapids.png) will help you understand how sleep episodes are chunked and grouped within time segments for the RAPIDS provider.
    1. Features listed in `[LEVELS_AND_TYPES]` are computed for any levels and types listed in `[SLEEP_LEVELS]` or `[SLEEP_TYPES]`. For example if `STAGES` only contains `[rem, light]` you will not get `countepisode[wake|deep][TYPE]` or sum, max, min, avg, median, or std `duration`. Levels or types in these lists do not influence `RATIOS` or `ROUTINE` features.
    2. Any `[LEVEL]` grouping is done within the elements of each class `CLASSIC`, `STAGES`, and `UNIFIED`. That is, we never combine `CLASSIC` or `STAGES` types to compute features.
    3. The categories for `all` levels (when `INCLUDE_ALL_GROUPS` is `True`) and `all` `SLEEP_TYPES` are not considered for `RATIOS` features as they are always 1.
    3. These features can be computed in time segments of any length, but only the 1-minute sleep chunks within each segment instance will be used.



## PRICE provider

!!! hint "Understanding PRICE features"
    [This diagram](../../img/sleep_intraday_price.png) will help you understand how sleep episodes are chunked and grouped within time segments and `LNE-LNE` intervals for the PRICE provider.

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
|`[SLEEP_LEVELS]`                             | Fitbit’s sleep API Version 1 only provides `CLASSIC` records. However, Version 1.2 provides 2 types of records: `CLASSIC` and `STAGES`. `STAGES` is only available in devices with a heart rate sensor and even those devices will fail to report it if the battery is low or the device is not tight enough. While `CLASSIC` contains 3 sleep levels (`awake`, `restless`, and `asleep`), `STAGES` contains 4 sleep levels (`wake`, `deep`, `light`, `rem`). To make it consistent, RAPIDS groups them into 2 `UNIFIED` sleep levels: `awake` (`CLASSIC`: `awake` and `restless`; `STAGES`: `wake`) and `asleep` (`CLASSIC`: `asleep`; `STAGES`: `deep`, `light`, and `rem`). In this section, there is a boolean flag named `INCLUDE_ALL_GROUPS` that if set to TRUE, computes avgdurationallmain`[DAY_TYPE]` features grouping all levels together in a single `all` category.
|`[DAY_TYPE]`                                 | The features of this provider can be computed using daily averages/standard deviations that were extracted on `WEEKEND` days only, `WEEK` days only, or `ALL` days|
|`[LAST_NIGHT_END]`                    | Only `main` sleep episodes that start within the `LNE-LNE` interval [`LAST_NIGHT_END`, `LAST_NIGHT_END` + 23H 59M 59S] are taken into account to compute the features described below. `[LAST_NIGHT_END]` is a number ranging from 0 (midnight) to 1439 (23:59). |


Features description for `[FITBIT_SLEEP_INTRADAY][PROVIDERS][PRICE]`:

|Feature&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                                  |Units          |Description                                                  |
|------------------------------------- |----------------- |-------------------------------------------------------------|
|avgduration`[LEVEL]`main`[DAY_TYPE]`             |minutes           | Average duration of daily sleep chunks of a `LEVEL`. Use the `DAY_TYPE` flag to include daily durations from weekend days only, weekdays, or both. Use `[LEVEL]` to group all levels in a single `all` category.
|avgratioduration`[LEVEL]`withinmain`[DAY_TYPE]`  |-                 | Average of the daily ratio between the duration of sleep chunks of a `LEVEL` and total duration of all `main` sleep episodes in a day. When `INCLUDE_ALL_GROUPS` is `True` the `all` `LEVEL` is ignored since this feature is always 1. Use the `DAY_TYPE` flag to include start times from weekend days only, weekdays, or both.
|avgstarttimeofepisodemain`[DAY_TYPE]` |minutes           | Average of all start times of the first `main` sleep episode within each `LNE-LNE` interval in a time segment. Use the `DAY_TYPE` flag to include start times from `LNE-LNE` intervals that start on weekend days only, weekdays, or both.
|avgendtimeofepisodemain`[DAY_TYPE]`   |minutes           | Average of all end times of the last `main` sleep episode within each `LNE-LNE` interval in a time segment. Use the `DAY_TYPE` flag to include end times from `LNE-LNE` intervals that start on weekend days only, weekdays, or both.
|avgmidpointofepisodemain`[DAY_TYPE]`  |minutes           | Average of all the differences between `avgendtime...` and `avgstarttime..` in a time segment. Use the `DAY_TYPE` flag to include end times from `LNE-LNE` intervals that start on weekend days only, weekdays, or both.
|stdstarttimeofepisodemain`[DAY_TYPE]` |minutes           | Standard deviation of all start times of the first `main` sleep episode within each `LNE-LNE` interval in a time segment. Use the `DAY_TYPE` flag to include start times from `LNE-LNE` intervals that start on weekend days only, weekdays, or both.
|stdendtimeofepisodemain`[DAY_TYPE]`   |minutes           | Standard deviation of all end times of the last `main` sleep episode within each `LNE-LNE` interval in a time segment. Use the `DAY_TYPE` flag to include end times from `LNE-LNE` intervals that start on weekend days only, weekdays, or both.
|stdmidpointofepisodemain`[DAY_TYPE]`  |minutes           | Standard deviation of all the differences between `avgendtime...` and `avgstarttime..` in a time segment. Use the `DAY_TYPE` flag to include end times from `LNE-LNE` intervals that start on weekend days only, weekdays, or both.
|socialjetlag                          |minutes           | Difference in minutes between the avgmidpointofepisodemain of weekends and weekdays that belong to each time segment instance. If your time segment does not contain at least one week day and one weekend day this feature will be NA. 
|rmssdmeanstarttimeofepisodemain       |minutes           | Square root of the **mean** squared successive difference (RMSSD) between today's and yesterday's `starttimeofepisodemain` values across the entire participant's sleep data grouped per time segment instance. It represents the mean of how someone's `starttimeofepisodemain` (bedtime) changed from night to night.
|rmssdmeanendtimeofepisodemain         |minutes           | Square root of the **mean** squared successive difference (RMSSD) between today's and yesterday's `endtimeofepisodemain` values across the entire participant's sleep data grouped per time segment instance. It represents the mean of how someone's `endtimeofepisodemain` (wake time) changed from night to night.
|rmssdmeanmidpointofepisodemain        |minutes           | Square root of the **mean** squared successive difference (RMSSD) between today's and yesterday's `midpointofepisodemain` values across the entire participant's sleep data grouped per time segment instance. It represents the mean of how someone's `midpointofepisodemain` (mid time between bedtime and wake time) changed from night to night.
|rmssdmedianstarttimeofepisodemain     |minutes           | Square root of the **median** squared successive difference (RMSSD) between today's and yesterday's `starttimeofepisodemain` values across the entire participant's sleep data grouped per time segment instance. It represents the median of how someone's `starttimeofepisodemain` (bedtime) changed from night to night.
|rmssdmedianendtimeofepisodemain       |minutes           | Square root of the **median** squared successive difference (RMSSD) between today's and yesterday's `endtimeofepisodemain` values across the entire participant's sleep data grouped per time segment instance. It represents the median of how someone's `endtimeofepisodemain` (wake time) changed from night to night.
|rmssdmedianmidpointofepisodemain      |minutes           | Square root of the **median** squared successive difference (RMSSD) between today's and yesterday's `midpointofepisodemain` values across the entire participant's sleep data grouped per time segment instance. It represents the median of how someone's `midpointofepisodemain` (average mid time between bedtime and wake time) changed from night to night.



!!! note "Assumptions/Observations"
    1. [This diagram](../../img/sleep_intraday_price.png) will help you understand how sleep episodes are chunked and grouped within time segments and `LNE-LNE` intervals for the PRICE provider.
    1. We recommend you use periodic segments that start in the morning so RAPIDS can chunk and group sleep episodes overnight. Shifted segments (as any other segments) are labelled based on their start and end date times.
    5. `avgstarttime...` and `avgendtime...` are roughly equivalent to an average bed and awake time only if you are using shifted segments.
    1. The features of this provider are only available on time segments that are longer than 24 hours because they are based on descriptive statistics computed across daily values.
    2. Even though Fitbit provides 2 types of sleep episodes (`main` and `nap`), only `main` sleep episodes are considered.
    4. The reference point for all times is 00:00 of the first day in the LNE-LNE interval.
    5. Sleep episodes are formed by 1-minute chunks that we group overnight starting from today’s LNE and ending on tomorrow’s LNE or the end of that segment (whatever is first). 
    5. The features `avgstarttime...` and `avgendtime...` are the average of the first and last sleep episode across every LNE-LNE interval within a segment (`avgmidtime...` is the mid point between start and end). Therefore, only segments longer than 24hrs will be averaged across more than one LNE-LNE interval.
    5. `socialjetlag` is only available on segment instances equal or longer than 48hrs that contain at least one weekday day and one weekend day, for example seven-day (weekly) segments.
