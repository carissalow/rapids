# `aware-mysql`

This [data stream](../../datastreams/data-streams-introduction) handles iOS and Android sensor data collected with the [AWARE Framework](https://awareframework.com/) and stored in a MySQL database.

## Container
A MySQL database with a table per sensor, each containing the data for all participants. This is the default database created by the old PHP AWARE server (as opposed to the new JavaScript Micro server).

The script to connect and download data from this container is at:
```bash
src/data/streams/aware_mysql/container.R
```

## Format
If you collected sensor data with the vanilla (original) AWARE mobile clients you shouldn't need to modify this format (described below). 

Remember that a format maps and transforms columns in your raw data stream to the [mandatory columns RAPIDS needs](../mandatory-phone-format).

The yaml file that describes the format of this data stream is at:
```bash
src/data/streams/aware_mysql/format.yaml
```

!!! hint
    The mappings in this stream (RAPIDS/Stream) are the same names because AWARE data was the first stream RAPIDS supported, meaning that it considers AWARE column names the default.

??? info "PHONE_ACCELEROMETER"

    === "ANDROID"
    
        **COLUMN_MAPPINGS**

        | RAPIDS column   | Stream column   |
        |-----------------|-----------------|
        | TIMESTAMP       | timestamp       |
        | DEVICE_ID       | device_id       |
        | DOUBLE_VALUES_0 | double_values_0 |
        | DOUBLE_VALUES_1 | double_values_1 |
        | DOUBLE_VALUES_2 | double_values_2 |

        **MUTATION_SCRIPTS**

        None

    === "IOS"
    
        Same as ANDROID
