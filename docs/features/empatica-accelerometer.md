# Empatica Accelerometer

Sensor parameters description for `[EMPATICA_ACCELEROMETER]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Name of the CSV file containing accelerometer data that is compressed inside an Empatica zip file. Since these zip files are created [automatically](https://support.empatica.com/hc/en-us/articles/201608896-Data-export-and-formatting-from-E4-connect-) by Empatica, there is no need to change the value of this attribute.

## DBDP provider

!!! info "Available time segments and platforms"
    - Available for all time segments

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/empatica_accelerometer_raw.csv
    - data/raw/{pid}/empatica_accelerometer_with_datetime.csv
    - data/interim/{pid}/empatica_accelerometer_features/empatica_accelerometer_{language}_{provider_key}.csv
    - data/processed/features/{pid}/empatica_accelerometer.csv
    ```


Parameters description for `[EMPATICA_ACCELEROMETER][PROVIDERS][DBDP]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`| Set to `True` to extract `EMPATICA_ACCELEROMETER` features from the `DBDP` provider|
|`[FEATURES]` |         Features to be computed, see table below


Features description for `[EMPATICA_ACCELEROMETER][PROVIDERS][RAPDBDPIDS]`:

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
|maxmagnitude      |m/s^2^    |The maximum magnitude of acceleration ($\|acceleration\| = \sqrt{x^2 + y^2 + z^2}$).
|minmagnitude      |m/s^2^    |The minimum magnitude of acceleration.
|avgmagnitude      |m/s^2^    |The average magnitude of acceleration.
|medianmagnitude   |m/s^2^    |The median magnitude of acceleration.
|stdmagnitude      |m/s^2^    |The standard deviation of acceleration.

!!! note "Assumptions/Observations"
    1. Analyzing accelerometer data is a memory intensive task. If RAPIDS crashes is likely because the accelerometer dataset for a participant is too big to fit in memory. We are considering different alternatives to overcome this problem, if this is something you need, get in touch and we can discuss how to implement it.
