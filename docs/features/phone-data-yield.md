# Phone Data Yield

This is a combinatorial sensor which means that we use the data from multiple sensors to extract data yield features. Data yield features can be used to remove rows ([time segments](../../setup/configuration/#time-segments)) that do not contain enough data. You should decide what is your "enough" threshold depending on the type of sensors you collected (frequency vs event based, e.g. acceleroemter vs calls), the length of your study, and the rates of missing data that your analysis could handle.

!!! hint "Why is data yield important?"
    Imagine that you want to extract `PHONE_CALL` features on daily segments (`00:00` to `23:59`). Let's say that on day 1 the phone logged 10 calls and 23 hours of data from other sensors and on day 2 the phone logged 10 calls and only 2 hours of data from other sensors. It's more likely that other calls were placed on the 22 hours of data that you didn't log on day 2 than on the 1 hour of data you didn't log on day 1, and so including day 2 in your analysis could bias your results.

Sensor parameters description for `[PHONE_DATA_YIELD]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;          | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[SENSORS]`| One or more phone sensor config keys (e.g. `PHONE_MESSAGE`). The more keys you include the more accurately RAPIDS can approximate the time an smartphone was sensing data. The supported phone sensors you can include in this list are outlined below (**do NOT include Fitbit sensors, ONLY include phone sensors**).

!!! info "Supported phone sensors for `[PHONE_DATA_YIELD][SENSORS]`"
    ```yaml
    PHONE_ACCELEROMETER
    PHONE_ACTIVITY_RECOGNITION
    PHONE_APPLICATIONS_CRASHES
    PHONE_APPLICATIONS_FOREGROUND
    PHONE_APPLICATIONS_NOTIFICATIONS
    PHONE_BATTERY
    PHONE_BLUETOOTH
    PHONE_CALLS
    PHONE_CONVERSATION
    PHONE_KEYBOARD
    PHONE_LIGHT
    PHONE_LOCATIONS
    PHONE_LOG
    PHONE_MESSAGES
    PHONE_SCREEN
    PHONE_WIFI_CONNECTED
    PHONE_WIFI_VISIBLE
    ```

## RAPIDS provider

Before explaining the data yield features, let's define the following relevant concepts:

- A valid minute is any 60 second window when any phone sensor logged at least 1 row of data
- A valid hour is any 60 minute window with at least X valid minutes. The X or threshold is given by `[MINUTE_RATIO_THRESHOLD_FOR_VALID_YIELDED_HOURS]`

The timestamps of all sensors are concatenated and then grouped per time segment. Minute and hour windows are created from the beginning of each time segment instance and these windows are marked as valid based on the definitions above. The duration of each time segment is taken into account to compute the features described below.

!!! info "Available time segments and platforms"
    - Available for all time segments
    - Available for Android and iOS

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/{sensor}_raw.csv # one for every [PHONE_DATA_YIELD][SENSORS]
    - data/interim/{pid}/phone_yielded_timestamps.csv
    - data/interim/{pid}/phone_yielded_timestamps_with_datetime.csv
    - data/interim/{pid}/phone_data_yield_features/phone_data_yield_{language}_{provider_key}.csv
    - data/processed/features/{pid}/phone_data_yield.csv
    ```


Parameters description for `[PHONE_DATA_YIELD][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`| Set to `True` to extract `PHONE_DATA_YIELD` features from the `RAPIDS` provider|
|`[FEATURES]` |  Features to be computed, see table below
|`[MINUTE_RATIO_THRESHOLD_FOR_VALID_YIELDED_HOURS]` | The proportion `[0.0 ,1.0]` of valid minutes in a 60-minute window necessary to flag that window as valid.


Features description for `[PHONE_DATA_YIELD][PROVIDERS][RAPIDS]`:

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
|ratiovalidyieldedminutes   |-          | The ratio between the number of valid minutes and the duration in minutes of a time segment.
|ratiovalidyieldedhours     |-          | The ratio between the number of valid hours and the duration in hours of a time segment. If the time segment is shorter than 1 hour this feature will always be 1.


!!! note "Assumptions/Observations"
    1. We recommend using `ratiovalidyieldedminutes` on time segments that are shorter than two or three hours and `ratiovalidyieldedhours` for longer segments. This is because relying on yielded minutes only can be misleading when a big chunk of those missing minutes are clustered together. 
    
        For example, let's assume we are working with a 24-hour time segment that is missing 12 hours of data. Two extreme cases can occur: 

        <ol type="A">
        <li>the 12 missing hours are from the beginning of the segment or </li>
        <li>30 minutes could be missing from every hour (24 * 30 minutes = 12 hours).</li>
        </ol>
        
        `ratiovalidyieldedminutes` would be 0.5 for both `a` and `b` (hinting the missing circumstances are similar). However, `ratiovalidyieldedhours` would be 0.5 for `a` and 1.0 for `b` if `[MINUTE_RATIO_THRESHOLD_FOR_VALID_YIELDED_HOURS]` is between [0.0 and 0.49] (hinting that the missing circumstances might be more favorable for `b`. In other words, sensed data for `b` is more evenly spread compared to `a`.
