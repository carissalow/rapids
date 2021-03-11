# `empatica_zip`
This [data stream](../../datastreams/data-streams-introduction) handles Empatica sensor data downloaded as zip files using the [E4 Connect](https://support.empatica.com/hc/en-us/articles/201608896-Data-export-and-formatting-from-E4-connect-). 

## Container

You need to create a subfolder for every participant named after their `device id` inside the folder specified by `[EMPATICA_DATA_STREAMS][empatica_zipfiles][FOLDER]`. You can add one or more Empatica zip files to any subfolder. 

The script to connect and download data from this container is at:
```bash
src/data/streams/empatica_zip/container.R
```

## Format


The `format.yaml` maps and transforms columns in your raw data stream to the [mandatory columns RAPIDS needs for Empatica sensors](../mandatory-empatica-format). This file is at:

```bash
src/data/streams/empatica_zip/format.yaml
```

All columns are mutated from the raw data in the zip files so you don't need to modify any column mappings.

??? info "EMPATICA_ACCELEROMETER"

    
    **RAPIDS_COLUMN_MAPPINGS**

    | RAPIDS column   | Stream column   |
    |-----------------|-----------------|
    | TIMESTAMP | timestamp|
    | DEVICE_ID | device_id|
    | DOUBLE_VALUES_0 | double_values_0|
    | DOUBLE_VALUES_1 | double_values_1|
    | DOUBLE_VALUES_2 | double_values_2|

    **MUTATION**

    - **COLUMN_MAPPINGS** (None)
    - **SCRIPTS** (None)

??? info "EMPATICA_HEARTRATE"

    
    **RAPIDS_COLUMN_MAPPINGS**

    | RAPIDS column   | Stream column   |
    |-----------------|-----------------|
    |TIMESTAMP | timestamp|
    |DEVICE_ID | device_id|
    |HEARTRATE | heartrate|

    **MUTATION**

    - **COLUMN_MAPPINGS** (None)
    - **SCRIPTS** (None)

??? info "EMPATICA_TEMPERATURE"

    
    **RAPIDS_COLUMN_MAPPINGS**

    | RAPIDS column   | Stream column   |
    |-----------------|-----------------|
    |TIMESTAMP | timestamp|
    |DEVICE_ID | device_id|
    |TEMPERATURE | temperature|

    **MUTATION**

    - **COLUMN_MAPPINGS** (None)
    - **SCRIPTS** (None)

??? info "EMPATICA_ELECTRODERMAL_ACTIVITY"

    
    **RAPIDS_COLUMN_MAPPINGS**

    | RAPIDS column   | Stream column   |
    |-----------------|-----------------|
    |TIMESTAMP | timestamp|
    |DEVICE_ID | device_id|
    |ELECTRODERMAL_ACTIVITY | electrodermal_activity|

    **MUTATION**

    - **COLUMN_MAPPINGS** (None)
    - **SCRIPTS** (None)

??? info "EMPATICA_BLOOD_VOLUME_PULSE"

    
    **RAPIDS_COLUMN_MAPPINGS**

    | RAPIDS column   | Stream column   |
    |-----------------|-----------------|
    |TIMESTAMP | timestamp|
    |DEVICE_ID | device_id|
    |BLOOD_VOLUME_PULSE | blood_volume_pulse|

    **MUTATION**

    - **COLUMN_MAPPINGS** (None)
    - **SCRIPTS** (None)

??? info "EMPATICA_INTER_BEAT_INTERVAL"

    
    **RAPIDS_COLUMN_MAPPINGS**

    | RAPIDS column   | Stream column   |
    |-----------------|-----------------|
    |TIMESTAMP | timestamp|
    |DEVICE_ID | device_id|
    |INTER_BEAT_INTERVAL | inter_beat_interval|

    **MUTATION**

    - **COLUMN_MAPPINGS** (None)
    - **SCRIPTS** (None)

??? info "EMPATICA_EMPATICA_TAGS"

    
    **RAPIDS_COLUMN_MAPPINGS**

    | RAPIDS column   | Stream column   |
    |-----------------|-----------------|
    |TIMESTAMP | timestamp|
    |DEVICE_ID | device_id|
    |TAGS | tags|

    **MUTATION**

    - **COLUMN_MAPPINGS** (None)
    - **SCRIPTS** (None)