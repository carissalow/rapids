# Phone Bluetooth

Sensor parameters description for `[PHONE_BLUETOOTH]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[TABLE]`| Database table where the bluetooth data is stored

## RAPIDS provider

!!! info "Available day segments and platforms"
    - Available for all day segments
    - Available for Android only

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/phone_bluetooth_raw.csv
    - data/raw/{pid}/phone_bluetooth_with_datetime.csv
    - data/interim/{pid}/phone_bluetooth_features/phone_bluetooth_{language}_{provider_key}.csv
    - data/processed/features/{pid}/phone_bluetooth.csv"
    ```


Parameters description for `[PHONE_BLUETOOTH][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`| Set to `True` to extract `PHONE_BLUETOOTH` features from the `RAPIDS` provider|
|`[FEATURES]` |         Features to be computed, see table below


Features description for `[PHONE_BLUETOOTH][PROVIDERS][RAPIDS]`:

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
| countscans                 | devices | Number of scanned devices during a `day_segment`, a device can be detected multiple times over time and these appearances are counted separately |
| uniquedevices              | devices | Number of unique devices during a `day_segment` as identified by their hardware (`bt_address`) address                                                          |
| countscansmostuniquedevice | scans   | Number of scans of the most scanned device during a `day_segment` across the whole monitoring period                                             |

!!! note "Assumptions/Observations"
    NA
