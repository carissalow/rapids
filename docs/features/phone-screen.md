# Phone Screen

Sensor parameters description for `[PHONE_SCREEN]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Data stream [container](../../datastreams/data-streams-introduction/) (database table, CSV file, etc.) where the screen data is stored

## RAPIDS provider

!!! info "Available time segments and platforms"
    - Available for all time segments
    - Available for Android and iOS

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/phone_screen_raw.csv
    - data/raw/{pid}/phone_screen_with_datetime.csv
    - data/interim/{pid}/phone_screen_episodes.csv
    - data/interim/{pid}/phone_screen_episodes_resampled.csv
    - data/interim/{pid}/phone_screen_episodes_resampled_with_datetime.csv
    - data/interim/{pid}/phone_screen_features/phone_screen_{language}_{provider_key}.csv
    - data/processed/features/{pid}/phone_screen.csv
    ```


Parameters description for `[PHONE_SCREEN][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`| Set to `True` to extract `PHONE_SCREEN` features from the `RAPIDS` provider|
|`[FEATURES]` |         Features to be computed, see table below
|`[REFERENCE_HOUR_FIRST_USE]` |  The reference point from which `firstuseafter` is to be computed, default is midnight
|`[IGNORE_EPISODES_SHORTER_THAN]` |  Ignore episodes that are shorter than this threshold (minutes). Set to 0 to disable this filter.
|`[IGNORE_EPISODES_LONGER_THAN]` |  Ignore episodes that are longer than this threshold (minutes), default is 6 hours. Set to 0 to disable this filter.
|`[EPISODE_TYPES]` |  Currently we only support `unlock` episodes (from when the phone is unlocked until the screen is off)


Features description for `[PHONE_SCREEN][PROVIDERS][RAPIDS]`:

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
|sumduration               |minutes           |Total duration of all unlock episodes.
|maxduration               |minutes           |Longest duration of any unlock episode.
|minduration               |minutes           |Shortest duration of any unlock episode.
|avgduration               |minutes           |Average duration of all unlock episodes.
|stdduration               |minutes           |Standard deviation duration of all unlock episodes.
|countepisode              |episodes          |Number of all unlock episodes
|firstuseafter             |minutes           |Minutes until the first unlock episode.

<!-- |episodepersensedminutes   |episodes/minute   |The ratio between the total number of episodes in an epoch divided by the total time (minutes) the phone was sensing data. -->

!!! note "Assumptions/Observations"
    1. In Android, `lock` events can happen right after an `off` event, after a few seconds of an `off` event, or never happen depending on the phone\'s settings, therefore, an `unlock` episode is defined as the time between an `unlock` and a `off` event. In iOS, `on` and `off` events do not exist, so an `unlock` episode is defined as the time between an `unlock` and a `lock` event.

    2. Events in iOS are recorded reliably albeit some duplicated `lock` events within milliseconds from each other, so we only keep consecutive unlock/lock pairs. In Android you cand find multiple consecutive `unlock` or `lock` events, so we only keep consecutive unlock/off pairs. In our experiments these cases are less than 10% of the screen events collected and this happens because `ACTION_SCREEN_OFF` and `ACTION_SCREEN_ON` are `sent when the device becomes non-interactive which may have nothing to do with the screen turning off`. In addition to unlock/off episodes, in Android it is possible to measure the time spent on the lock screen before an `unlock` event as well as the total screen time (i.e. `ON` to `OFF`) but these are not implemented at the moment.

    3. We transform iOS screen events to match Android's format, we replace `lock` episodes with `off` episodes (2 with 0) in iOS. However, as mentioned above this is still computing `unlock` to `lock` episodes.
