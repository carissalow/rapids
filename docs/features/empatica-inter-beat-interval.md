# Empatica Inter Beat Interval

Sensor parameters description for `[EMPATICA_INTER_BEAT_INTERVAL]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Name of the CSV file containing inter beat interval data that is compressed inside an Empatica zip file. Since these zip files are created [automatically](https://support.empatica.com/hc/en-us/articles/201608896-Data-export-and-formatting-from-E4-connect-) by Empatica, there is no need to change the value of this attribute.

## DBDP provider

!!! info "Available time segments and platforms"
    - Available for all time segments

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/empatica_inter_beat_interval_raw.csv
    - data/raw/{pid}/empatica_inter_beat_interval_with_datetime.csv
    - data/interim/{pid}/empatica_inter_beat_interval_features/empatica_inter_beat_interval_{language}_{provider_key}.csv
    - data/processed/features/{pid}/empatica_inter_beat_interval.csv
    ```


Parameters description for `[EMPATICA_INTER_BEAT_INTERVAL][PROVIDERS][DBDP]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`  | Set to `True` to extract `EMPATICA_INTER_BEAT_INTERVAL` features from the `DBDP` provider|
|`[FEATURES]` |         Features to be computed from inter beat interval intraday data, see table below          |


Features description for `[EMPATICA_INTER_BEAT_INTERVAL][PROVIDERS][DBDP]`:

|Feature                    |Units          |Description|
|-------------------------- |-------------- |---------------------------|
|maxibi                      |seconds     |The maximum inter beat interval during a time segment.
|minibi                      |seconds     |The minimum inter beat interval during a time segment.
|avgibi                      |seconds     |The average inter beat interval during a time segment.
|medianibi                   |seconds     |The median of inter beat interval during a time segment.
|modeibi                     |seconds     |The mode of inter beat interval during a time segment.
|stdibi                      |seconds     |The standard deviation of inter beat interval during a time segment.
|diffmaxmodeibi              |seconds     |The difference between the maximum and mode inter beat interval during a time segment.
|diffminmodeibi              |seconds     |The difference between the mode and minimum inter beat interval during a time segment.
|entropyibi                  |nats           |Shannon’s entropy measurement based on inter beat interval during a time segment.

!!! note "Assumptions/Observations"
    For more information about IBI read [this](https://support.empatica.com/hc/en-us/articles/360030058011-E4-data-IBI-expected-signal).