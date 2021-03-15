# Empatica Electrodermal Activity

Sensor parameters description for `[EMPATICA_ELECTRODERMAL_ACTIVITY]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Name of the CSV file containing electrodermal activity data that is compressed inside an Empatica zip file. Since these zip files are created [automatically](https://support.empatica.com/hc/en-us/articles/201608896-Data-export-and-formatting-from-E4-connect-) by Empatica, there is no need to change the value of this attribute.

## DBDP provider

!!! info "Available time segments and platforms"
    - Available for all time segments

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/empatica_electrodermal_activity_raw.csv
    - data/raw/{pid}/empatica_electrodermal_activity_with_datetime.csv
    - data/interim/{pid}/empatica_electrodermal_activity_features/empatica_electrodermal activity_{language}_{provider_key}.csv
    - data/processed/features/{pid}/empatica_electrodermal_activity.csv
    ```


Parameters description for `[EMPATICA_ELECTRODERMAL_ACTIVITY][PROVIDERS][DBDP]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`  | Set to `True` to extract `EMPATICA_ELECTRODERMAL_ACTIVITY` features from the `DBDP` provider|
|`[FEATURES]` |         Features to be computed from electrodermal activity intraday data, see table below          |


Features description for `[EMPATICA_ELECTRODERMAL ACTIVITY][PROVIDERS][DBDP]`:

|Feature                    |Units          |Description|
|-------------------------- |-------------- |---------------------------|
|maxeda                      |microsiemens     |The maximum electrical conductance during a time segment.
|mineda                      |microsiemens     |The minimum electrical conductance during a time segment.
|avgeda                      |microsiemens     |The average electrical conductance during a time segment.
|medianeda                   |microsiemens     |The median of electrical conductance during a time segment.
|modeeda                     |microsiemens     |The mode of electrical conductance during a time segment.
|stdeda                      |microsiemens     |The standard deviation of electrical conductance during a time segment.
|diffmaxmodeeda              |microsiemens     |The difference between the maximum and mode electrical conductance during a time segment.
|diffminmodeeda              |microsiemens     |The difference between the mode and minimum electrical conductance during a time segment.
|entropyeda                  |nats           |Shannon’s entropy measurement based on electrical conductance during a time segment.

!!! note "Assumptions/Observations"
    None