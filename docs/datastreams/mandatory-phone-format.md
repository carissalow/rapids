# Mandatory Phone Format

This is a description of the format RAPIDS needs to process data for the following PHONE sensors.

??? info "PHONE_ACCELEROMETER"

    | RAPIDS column   | Description                                                  |
    |-----------------|--------------------------------------------------------------|
    | TIMESTAMP       | An UNIX timestamp (13 digits) when a row of data was logged  |
    | DEVICE_ID       | A string that uniquely identifies a device                   |
    | DOUBLE_VALUES_0 | x axis of acceleration                                       |
    | DOUBLE_VALUES_1 | y axis of acceleration                                       |
    | DOUBLE_VALUES_2 | z axis of acceleration                                       |


??? info "PHONE_ACTIVITY_RECOGNITION"

    | RAPIDS column   | Description                                                               |
    |-----------------|---------------------------------------------------------------------------|
    | TIMESTAMP       | An UNIX timestamp (13 digits) when a row of data was logged               |
    | DEVICE_ID       | A string that uniquely identifies a device                                |
    | ACTIVITY_TYPE   | An integer (ranged from 0 to 8) that denotes current activity type        |
    | ACTIVITY_NAME   | An string that denotes current activity name: `in_vehicle`, `on_bicycle`, `on_foot`, `still`, `unknown`, `tilting`, `walking` or `running`   |
    | CONFIDENCE      | An integer (ranged from 0 to 100) that denotes the prediction accuracy    |


??? info "PHONE_CONVERSATION"

    | RAPIDS column        | Description                                                                          |
    |----------------------|--------------------------------------------------------------------------------------|
    | TIMESTAMP            | An UNIX timestamp (13 digits) when a row of data was logged                          |
    | DEVICE_ID            | A string that uniquely identifies a device                                           |
    | DOUBLE_ENERGY        | A number that denotes the amplitude of an audio sample (L2-norm of the audio frame)     |
    | INFERENCE            | An integer (ranged from 0 to 3) that denotes the type of an audio sample: 0 = silence, 1 = noise, 2 = voice, 3 = unknown      |
    | DOUBLE_CONVO_START   | UNIX timestamp (13 digits) of the beginning of a conversation                        |
    | DOUBLE_CONVO_END     | UNIX timestamp (13 digits) of the end of a conversation                              |
