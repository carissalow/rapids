# Empatica Blood Volume Pulse

Sensor parameters description for `[EMPATICA_BLOOD_VOLUME_PULSE]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Name of the CSV file containing blood volume pulse data that is compressed inside an Empatica zip file. Since these zip files are created [automatically](https://support.empatica.com/hc/en-us/articles/201608896-Data-export-and-formatting-from-E4-connect-) by Empatica, there is no need to change the value of this attribute.

## DBDP provider

!!! info "Available time segments and platforms"
    - Available for all time segments

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/empatica_blood_volume_pulse_raw.csv 
    - data/raw/{pid}/empatica_blood_volume_pulse_with_datetime.csv
    - data/interim/{pid}/empatica_blood_volume_pulse_features/empatica_blood_volume_pulse_{language}_{provider_key}.csv
    - data/processed/features/{pid}/empatica_blood_volume_pulse.csv
    ```


Parameters description for `[EMPATICA_BLOOD_VOLUME_PULSE][PROVIDERS][DBDP]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`  | Set to `True` to extract `EMPATICA_BLOOD_VOLUME_PULSE` features from the `DBDP` provider|
|`[FEATURES]` |         Features to be computed from blood volume pulse intraday data, see table below          |


Features description for `[EMPATICA_BLOOD_VOLUME_PULSE][PROVIDERS][DBDP]`:

|Feature                    |Units          |Description|
|-------------------------- |-------------- |---------------------------|
|maxbvp                      |-     |The maximum blood volume pulse during a time segment.
|minbvp                      |-     |The minimum blood volume pulse during a time segment.
|avgbvp                      |-     |The average blood volume pulse during a time segment.
|medianbvp                   |-     |The median of blood volume pulse during a time segment.
|modebvp                     |-     |The mode of blood volume pulse during a time segment.
|stdbvp                      |-     |The standard deviation of blood volume pulse during a time segment.
|diffmaxmodebvp              |-     |The difference between the maximum and mode blood volume pulse during a time segment.
|diffminmodebvp              |-     |The difference between the mode and minimum blood volume pulse during a time segment.
|entropybvp                  |nats           |Shannon’s entropy measurement based on blood volume pulse during a time segment.

!!! note "Assumptions/Observations"
    For more information about BVP read [this](https://support.empatica.com/hc/en-us/articles/360029719792-E4-data-BVP-expected-signal).