# Mandatory Phone Format

This is a description of the format RAPIDS needs to process data for the following PHONE sensors.

??? info "PHONE_ACCELEROMETER"

    === "ANDROID"

        | RAPIDS column   | Description   |
        |-----------------|-----------------|
        | TIMESTAMP       | A UNIX timestamp (13 digits) when a row of data was logged   |
        | DEVICE_ID       | A string that uniquely identifies a device       |
        | DOUBLE_VALUES_0 | x axis of acceleration |
        | DOUBLE_VALUES_1 | y axis of acceleration |
        | DOUBLE_VALUES_2 | z axis of acceleration |

    === "IOS"
        Same as ANDROID
    