# Phone Light

Sensor parameters description for `[PHONE_LIGHT]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Data stream [container](../../datastreams/data-streams-introduction/) (database table, CSV file, etc.) where the light data is stored

## RAPIDS provider

!!! info "Available time segments and platforms"
    - Available for all time segments
    - Available for Android only

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/phone_light_raw.csv
    - data/raw/{pid}/phone_light_with_datetime.csv
    - data/interim/{pid}/phone_light_features/phone_light_{language}_{provider_key}.csv
    - data/processed/features/{pid}/phone_light.csv
    ```


Parameters description for `[PHONE_LIGHT][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`| Set to `True` to extract `PHONE_LIGHT` features from the `RAPIDS` provider|
|`[FEATURES]` |         Features to be computed, see table below


Features description for `[PHONE_LIGHT][PROVIDERS][RAPIDS]`:

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
|count       |rows    | Number light sensor rows recorded.
|maxlux      |lux     | The maximum ambient luminance.
|minlux      |lux     | The minimum ambient luminance.
|avglux      |lux     | The average ambient luminance.
|medianlux   |lux     | The median ambient luminance.
|stdlux      |lux     | The standard deviation of ambient luminance.

!!! note "Assumptions/Observations"
    NA
