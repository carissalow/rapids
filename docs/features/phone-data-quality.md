# Phone Data Quality

## Phone Valid Sensed Bins

A valid bin is any period of `BIN_SIZE` minutes starting from midnight with at least 1 row from any phone sensor. `PHONE_VALID_SENSED_BINS` are used to compute `PHONE_VALID_SENSED_DAYS`, to resample fused location data and to compute some screen features.

!!! hint
    `PHONE_VALID_SENSED_DAYS` are an approximation to the time the phone was sensing data so add as many sensors as you have to `[PHONE_VALID_SENSED_BINS][PHONE_SENSORS]`

Parameters description for `PHONE_VALID_SENSED_BINS`:

| Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Description|
|-----|------------|
| `[COMPUTE]`| Set to `True` to compute |
| `[BIN_SIZE]` | Size of each bin in minutes | 
| `[PHONE_SENSORS]` | One or more sensor config keys (e.g. `PHONE_MESSAGE`) to be used to flag a bin as valid or not (whether or not a bin contains at least one row from any sensor)|

!!! info "Possible values for `[PHONE_VALID_SENSED_BINS][PHONE_SENSORS]`"
    ```yaml
    PHONE_MESSAGES
    PHONE_CALLS
    PHONE_LOCATIONS
    PHONE_BLUETOOTH
    PHONE_ACTIVITY_RECOGNITION
    PHONE_BATTERY
    PHONE_SCREEN
    PHONE_LIGHT
    PHONE_ACCELEROMETER
    PHONE_APPLICATIONS_FOREGROUND
    PHONE_WIFI_VISIBLE
    PHONE_WIFI_CONNECTED
    PHONE_CONVERSATION
    ```


## Phone Valid Sensed Days

On any given day, a phone could have sensed data only for a few minutes or for 24 hours. Features should considered more reliable the more hours the phone was logging data, for example, 10 calls logged on a day when only 1 hour of data was recorded is a less reliable feature compared to 10 calls on a day when 23 hours of data were recorded.

Therefore, we define a valid hour as those that contain a minimum number of valid bins (see above). We mark an hour as valid when contains at least `MIN_VALID_BINS_PER_HOUR` (out of 60min/`BIN_SIZE` bins). In turn, we mark a day as valid if it has at least `MIN_VALID_HOURS_PER_DAY`. You can use `PHONE_VALID_SENSED_DAYS` to manually discard days when not enough data was collected after your features are computed.

Parameters description for `PHONE_VALID_SENSED_DAYS`:

| Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Description|
|-----|------------|
| `[COMPUTE]`| Set to `True` to compute |
| `[MIN_VALID_BINS_PER_HOUR]` | An array of integer values, 6 by default. Minimum number of valid bins to mark an hour as valid|
| `[MIN_VALID_HOURS_PER_DAY]` | An array of integer values, 16 by default. Minimum number of valid hours to mark a day as valid | 

