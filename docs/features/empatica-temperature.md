# Empatica Temperature

Sensor parameters description for `[EMPATICA_TEMPERATURE]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Name of the CSV file containing temperature data that is compressed inside an Empatica zip file. Since these zip files are created [automatically](https://support.empatica.com/hc/en-us/articles/201608896-Data-export-and-formatting-from-E4-connect-) by Empatica, there is no need to change the value of this attribute.

## DBDP provider

!!! info "Available time segments and platforms"
    - Available for all time segments

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/empatica_temperature_raw.csv
    - data/raw/{pid}/empatica_temperature_with_datetime.csv
    - data/interim/{pid}/empatica_temperature_features/empatica_temperature_{language}_{provider_key}.csv
    - data/processed/features/{pid}/empatica_temperature.csv
    ```


Parameters description for `[EMPATICA_TEMPERATURE][PROVIDERS][DBDP]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`  | Set to `True` to extract `EMPATICA_TEMPERATURE` features from the `DBDP` provider|
|`[FEATURES]` |         Features to be computed from temperature intraday data, see table below          |


Features description for `[EMPATICA_TEMPERATURE][PROVIDERS][DBDP]`:

|Feature                    |Units          |Description|
|-------------------------- |-------------- |---------------------------|
|maxtemp                      |degrees C     |The maximum temperature during a time segment.
|mintemp                      |degrees C     |The minimum temperature during a time segment.
|avgtemp                      |degrees C     |The average temperature during a time segment.
|mediantemp                   |degrees C     |The median of temperature during a time segment.
|modetemp                     |degrees C     |The mode of temperature during a time segment.
|stdtemp                      |degrees C     |The standard deviation of temperature during a time segment.
|diffmaxmodetemp              |degrees C     |The difference between the maximum and mode temperature during a time segment.
|diffminmodetemp              |degrees C     |The difference between the mode and minimum temperature during a time segment.
|entropytemp                  |nats           |Shannon’s entropy measurement based on temperature during a time segment.

!!! note "Assumptions/Observations"
    None