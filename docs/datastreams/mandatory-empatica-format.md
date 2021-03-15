# Mandatory Empatica Format

This is a description of the format RAPIDS needs to process data for the following Empatica sensors.

??? info "EMPATICA_ACCELEROMETER"

    | RAPIDS column   | Description                                                  |
    |-----------------|--------------------------------------------------------------|
    | TIMESTAMP       | An UNIX timestamp (13 digits) when a row of data was logged  |
    | DEVICE_ID       | A string that uniquely identifies a device                   |
    | DOUBLE_VALUES_0 | x axis of acceleration                                       |
    | DOUBLE_VALUES_1 | y axis of acceleration                                       |
    | DOUBLE_VALUES_2 | z axis of acceleration                                       |

??? info "EMPATICA_HEARTRATE"

    | RAPIDS column   | Description   |
    |-----------------|-----------------|
    | TIMESTAMP       |  An UNIX timestamp (13 digits) when a row of data was logged (automatically created by RAPIDS) |
    | DEVICE_ID       |  A string that uniquely identifies a device |
    | HEARTRATE |  Intraday heartrate |

??? info "EMPATICA_TEMPERATURE"

    | RAPIDS column   | Description   |
    |-----------------|-----------------|
    | TIMESTAMP       |  An UNIX timestamp (13 digits) when a row of data was logged (automatically created by RAPIDS) |
    | DEVICE_ID       |  A string that uniquely identifies a device |
    | TEMPERATURE |  temperature |

??? info "EMPATICA_ELECTRODERMAL_ACTIVITY"

    | RAPIDS column   | Description   |
    |-----------------|-----------------|
    | TIMESTAMP       |  An UNIX timestamp (13 digits) when a row of data was logged (automatically created by RAPIDS) |
    | DEVICE_ID       |  A string that uniquely identifies a device |
    | ELECTRODERMAL_ACTIVITY |  electrical conductance |

??? info "EMPATICA_BLOOD_VOLUME_PULSE"

    | RAPIDS column   | Description   |
    |-----------------|-----------------|
    | TIMESTAMP       |  An UNIX timestamp (13 digits) when a row of data was logged (automatically created by RAPIDS) |
    | DEVICE_ID       |  A string that uniquely identifies a device |
    | BLOOD_VOLUME_PULSE |  blood volume pulse |

??? info "EMPATICA_INTER_BEAT_INTERVAL"

    | RAPIDS column   | Description   |
    |-----------------|-----------------|
    | TIMESTAMP       |  An UNIX timestamp (13 digits) when a row of data was logged (automatically created by RAPIDS) |
    | DEVICE_ID       |  A string that uniquely identifies a device |
    | INTER_BEAT_INTERVAL |  inter beat interval |

??? info "EMPATICA_TAGS"

    | RAPIDS column   | Description   |
    |-----------------|-----------------|
    | TIMESTAMP       |  An UNIX timestamp (13 digits) when a row of data was logged (automatically created by RAPIDS) |
    | DEVICE_ID       |  A string that uniquely identifies a device |
    | TAGS |  tags |
