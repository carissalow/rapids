# Phone WiFi Visible

Sensor parameters description for `[PHONE_WIFI_VISIBLE]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Data stream [container](../../datastreams/data-streams-introduction/) (database table, CSV file, etc.) where the wifi (visible) data is stored

## RAPIDS provider

!!! info "Available time segments and platforms"
    - Available for all time segments
    - Available for Android only

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/phone_wifi_visible_raw.csv
    - data/raw/{pid}/phone_wifi_visible_with_datetime.csv
    - data/interim/{pid}/phone_wifi_visible_features/phone_wifi_visible_{language}_{provider_key}.csv
    - data/processed/features/{pid}/phone_wifi_visible.csv
    ```


Parameters description for `[PHONE_WIFI_VISIBLE][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`| Set to `True` to extract `PHONE_WIFI_VISIBLE` features from the `RAPIDS` provider|
|`[FEATURES]` |         Features to be computed, see table below


Features description for `[PHONE_WIFI_VISIBLE][PROVIDERS][RAPIDS]`:

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
| countscans                 | devices | Number of scanned WiFi access points visible during a time_segment, an access point can be detected multiple times over time and these appearances are counted separately |
| uniquedevices              | devices | Number of unique access point during a time_segment as identified by their hardware address                                                                       |
| countscansmostuniquedevice | scans   | Number of scans of the most scanned access point during a time_segment across the whole monitoring period                                                         |

!!! note "Assumptions/Observations"
    1. A visible WiFI access point is one that a phone sensed around itself but that it was not connected to. Due to API restrictions, this sensor is not available on iOS.
    2. By default AWARE stores this data in the `wifi` table.
