# Empatica Heart Rate

Sensor parameters description for `[EMPATICA_HEARTRATE]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Name of the CSV file containing heart rate data that is compressed inside an Empatica zip file. Since these zip files are created [automatically](https://support.empatica.com/hc/en-us/articles/201608896-Data-export-and-formatting-from-E4-connect-) by Empatica, there is no need to change the value of this attribute.

## DBDP provider

!!! info "Available time segments and platforms"
    - Available for all time segments

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/empatica_heartrate_raw.csv
    - data/raw/{pid}/empatica_heartrate_with_datetime.csv
    - data/interim/{pid}/empatica_heartrate_features/empatica_heartrate_{language}_{provider_key}.csv
    - data/processed/features/{pid}/empatica_heartrate.csv
    ```


Parameters description for `[EMPATICA_HEARTRATE][PROVIDERS][DBDP]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`  | Set to `True` to extract `EMPATICA_HEARTRATE` features from the `DBDP` provider|
|`[FEATURES]` |         Features to be computed from heart rate intraday data, see table below          |


Features description for `[EMPATICA_HEARTRATE][PROVIDERS][DBDP]`:

|Feature                    |Units          |Description|
|-------------------------- |-------------- |---------------------------|
|maxhr                      |beats     |The maximum heart rate during a time segment.
|minhr                      |beats     |The minimum heart rate during a time segment.
|avghr                      |beats     |The average heart rate during a time segment.
|medianhr                   |beats     |The median of heart rate during a time segment.
|modehr                     |beats     |The mode of heart rate during a time segment.
|stdhr                      |beats     |The standard deviation of heart rate during a time segment.
|diffmaxmodehr              |beats     |The difference between the maximum and mode heart rate during a time segment.
|diffminmodehr              |beats     |The difference between the mode and minimum heart rate during a time segment.
|entropyhr                  |nats           |Shannon’s entropy measurement based on heart rate during a time segment.

!!! note "Assumptions/Observations"
    We extract the previous features based on the average heart rate values computed in [10-second windows](https://support.empatica.com/hc/en-us/articles/360029469772-E4-data-HR-csv-explanation).