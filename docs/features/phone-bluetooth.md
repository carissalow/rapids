# Phone Bluetooth

Sensor parameters description for `[PHONE_BLUETOOTH]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[TABLE]`| Database table where the bluetooth data is stored

## RAPIDS provider

!!! info "Available time segments and platforms"
    - Available for all time segments
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
| countscans                 | devices | Number of scanned devices during a `time_segment`, a device can be detected multiple times over time and these appearances are counted separately |
| uniquedevices              | devices | Number of unique devices during a `time_segment` as identified by their hardware (`bt_address`) address                                                          |
| countscansmostuniquedevice | scans   | Number of scans of the most scanned device during a `time_segment` across the whole monitoring period                                             |

!!! note "Assumptions/Observations"
    NA

## DORYAB provider

!!! info "Available time segments and platforms"
    - Available for all time segments
    - Available for Android only

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/phone_bluetooth_raw.csv
    - data/raw/{pid}/phone_bluetooth_with_datetime.csv
    - data/interim/{pid}/phone_bluetooth_features/phone_bluetooth_{language}_{provider_key}.csv
    - data/processed/features/{pid}/phone_bluetooth.csv"
    ```


Parameters description for `[PHONE_BLUETOOTH][PROVIDERS][DORYAB]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`| Set to `True` to extract `PHONE_BLUETOOTH` features from the `DORYAB` provider|
|`[FEATURES]` |         Features to be computed, see table below. These features are computed for three device categories: `all` devices, `own` devices and `other` devices.


Features description for `[PHONE_BLUETOOTH][PROVIDERS][DORYAB]`:

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
| countscans                 | scans | Number of scans (rows) from the devices sensed during a time segment instance. The more scans a bluetooth device has the longer it remained within range of the participant's phone |
| uniquedevices              | devices | Number of unique bluetooth devices sensed during a time segment instance as identified by their hardware addresses (`bt_address`) |
| countscansmostuniquedevice | scans   | Number of scans of the most sensed device within each time segment instance|
| countscansleastuniquedevice| scans| Number of scans of the least sensed device within each time segment instance |
| meanscans | scans| Mean of the scans of every sensed device within each time segment instance|
| stdscans | scans| Standard deviation of the scans of every sensed device within each time segment instance|

!!! note "Assumptions/Observations"
    - This provider is adapted from the work by [Doryab et al](../../citation#doryab-bluetooth). Devices are clasified as belonging to the participant (`own`) or to other people (`others`) using k-means based on the number of times and the number of days each device was detected across each participant's dataset.
    - If ownership cannot be computed because all devices were detected on only one day, they are all considered as `other`. Thus `all` and `other` features will be equal.
    - These features are computed for devices detected within each time segment instance. For example, let's say that we logged the following devices on three different time segment instances (days) for `p01`:
    ```csv
    local_date                            bt_address
    2016-11-29  55C836F5-487E-405F-8E28-21DBD40FA4FF
    2016-11-29  55C836F5-487E-405F-8E28-21DBD40FA4FF
    2016-11-29  55C836F5-487E-405F-8E28-21DBD40FA4FF
    2016-11-29  48872A52-68DE-420D-98DA-73339A1C4685
    2016-11-29  48872A52-68DE-420D-98DA-73339A1C4685
    2016-11-30  55C836F5-487E-405F-8E28-21DBD40FA4FF
    2016-11-30  55C836F5-487E-405F-8E28-21DBD40FA4FF
    2016-11-30  48872A52-68DE-420D-98DA-73339A1C4685
    2017-05-07  5C5A9C41-2F68-4CEB-96D0-77DE3729B729
    2017-05-07  25262DC7-780C-4AD5-AD3A-D9776AEF7FC1
    2017-05-07  5B1E6981-2E50-4D9A-99D8-67AED430C5A8
    2017-05-07  6C444841-FE64-4375-BC3F-FA410CDC0AC7
    2017-05-07  5B1E6981-2E50-4D9A-99D8-67AED430C5A8
    2017-05-07  4DC7A22D-9F1F-4DEF-8576-086910AABCB5
    ```
    - For each device we compute `days_scanned` (the number of days on which each device was detected), `scans` (the number of times each device was detected), `scans_per_day` that's equal to `scans/days_scanned`, and whether a devices is labelled as `own` or `other` (note the last device is labelled as a `own` device because it was detected 6 times over two time segment instances):
    ```csv
    bt_address                            days_scanned  scans  scans_per_day own_device
    25262DC7-780C-4AD5-AD3A-D9776AEF7FC1             1      1            1.0          0
    4DC7A22D-9F1F-4DEF-8576-086910AABCB5             1      1            1.0          0
    5C5A9C41-2F68-4CEB-96D0-77DE3729B729             1      1            1.0          0
    6C444841-FE64-4375-BC3F-FA410CDC0AC7             1      1            1.0          0
    5B1E6981-2E50-4D9A-99D8-67AED430C5A8             1      2            2.0          0
    48872A52-68DE-420D-98DA-73339A1C4685             2      3            1.5          0
    55C836F5-487E-405F-8E28-21DBD40FA4FF             2      5            2.5          1
    ```
    - These are the metrics for each time instance (day) for `own` and `other` devices (we ignore `all` for brevity). The only `own` device (`55C836F5-487E-405F-8E28-21DBD40FA4FF`) was detected on the first two days, 3 and 2 times respectively, the `other` devices where detected on all three days. On the last day (`2017-05-07`) there were 6 scans from 5 unique devices, the most frequent device for that day was `5B1E6981-2E50-4D9A-99D8-67AED430C5A8` with 2 scans, and the mean number of scans among all devices was 1.2 (`[1 + 1 + 1 + 1 + 2] / 5`)
    ```csv
    local_segment countscansown uniquedevicesown countscansmostuniquedeviceown countscansleastuniquedeviceown meanscansown stdscansown countscansothers uniquedevicesothers countscansmostuniquedeviceothers countscansleastuniquedeviceothers meanscansothers stdscansothers
    2016-11-29 3.0 1.0 3.0 3.0 3.0 NaN 2 1 2 2 2.0 NaN
    2016-11-30 2.0 1.0 2.0 2.0 2.0 NaN 1 1 1 1 1.0 NaN
    2017-05-07 NaN NaN NaN NaN NaN NaN 6 5 2 1 1.2 0.447214
    ```
