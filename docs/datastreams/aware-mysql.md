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

Stream columns named `FLAG_TO_MUTATE` means they are extracted based on the `MUTATION` section. You can refer to [OS complex mapping](../../datastreams/add-new-data-streams/#os-complex-mapping) for detailed information.

!!! hint
    The mappings in this stream (RAPIDS/Stream) are the same names because AWARE data was the first stream RAPIDS supported, meaning that it considers AWARE column names the default.

??? info "PHONE_ACCELEROMETER"

    === "ANDROID"
    
        **RAPIDS_COLUMN_MAPPINGS**

        | RAPIDS column   | Stream column   |
        |-----------------|-----------------|
        | TIMESTAMP       | timestamp       |
        | DEVICE_ID       | device_id       |
        | DOUBLE_VALUES_0 | double_values_0 |
        | DOUBLE_VALUES_1 | double_values_1 |
        | DOUBLE_VALUES_2 | double_values_2 |

        **MUTATION**

        - **COLUMN_MAPPINGS** (None)
        - **SCRIPTS** (None)

    === "IOS"
    
        Same as ANDROID

??? info "PHONE_ACTIVITY_RECOGNITION"

    === "ANDROID"
    
        **RAPIDS_COLUMN_MAPPINGS**

        | RAPIDS column   | Stream column   |
        |-----------------|-----------------|
        | TIMESTAMP       | timestamp       |
        | DEVICE_ID       | device_id       |
        | ACTIVITY_NAME   | activity_name   |
        | ACTIVITY_TYPE   | activity_type   |
        | CONFIDENCE      | confidence      |

        **MUTATION**

        - **COLUMN_MAPPINGS** (None)
        - **SCRIPTS** (None)

    === "IOS"

        **RAPIDS_COLUMN_MAPPINGS**

        | RAPIDS column   | Stream column   |
        |-----------------|-----------------|
        | TIMESTAMP       | timestamp       |
        | DEVICE_ID       | device_id       |
        | ACTIVITY_NAME   | FLAG_TO_MUTATE  |
        | ACTIVITY_TYPE   | FLAG_TO_MUTATE  |
        | CONFIDENCE      | FLAG_TO_MUTATE  |

        **MUTATION**

        - **COLUMN_MAPPINGS**

        | Script column   | Stream column   |
        |-----------------|-----------------|
        | ACTIVITIES      | activities      |
        | CONFIDENCE      | confidence      |

        - **SCRIPTS**
        
        ```bash
        src/data/streams/mutations/phone/aware/activity_recogniton_ios_unification.R
        ```


        !!! note
            For RAPIDS columns of `ACTIVITY_NAME` and `ACTIVITY_TYPE`:

            - if stream's `activities` field is automotive, set `ACTIVITY_NAME` = in_vehicle and `ACTIVITY_TYPE` = 0
            - if stream's `activities` field is cycling, set `ACTIVITY_NAME` = on_bicycle and `ACTIVITY_TYPE` = 1
            - if stream's `activities` field is walking, set `ACTIVITY_NAME` = walking and `ACTIVITY_TYPE` = 7
            - if stream's `activities` field is running, set `ACTIVITY_NAME` = running and `ACTIVITY_TYPE` = 8
            - if stream's `activities` field is stationary, set `ACTIVITY_NAME` = still and `ACTIVITY_TYPE` = 3
            - if stream's `activities` field is unknown, set `ACTIVITY_NAME` = unknown and  `ACTIVITY_TYPE` = 4

            For RAPIDS `CONFIDENCE` column:
            
            - if stream's `confidence` field is 0, set `CONFIDENCE` = 0
            - if stream's `confidence` field is 1, set `CONFIDENCE` = 50
            - if stream's `confidence` field is 2, set `CONFIDENCE` = 100



??? info "PHONE_APPLICATIONS_FOREGROUND"

    === "ANDROID"
    
         **RAPIDS_COLUMN_MAPPINGS**

        | RAPIDS column      | Stream column      |
        |--------------------|--------------------|
        | TIMESTAMP          | timestamp          |
        | DEVICE_ID          | device_id          |
        | PACKAGE_NAME       | package_name       |
        | APPLICATION_NAME   | application_name   |
        | IS_SYSTEM_APP      | is_system_app      |

        **MUTATION**

        - **COLUMN_MAPPINGS** (None)
        - **SCRIPTS** (None)
    
    === "IOS"

        This sensor is not supported by iOS devices.

??? info "PHONE_BATTERY"

    === "ANDROID"
    
        **RAPIDS_COLUMN_MAPPINGS**

        | RAPIDS column        | Stream column       |
        |----------------------|---------------------|
        | TIMESTAMP            | timestamp           |
        | DEVICE_ID            | device_id           |
        | BATTERY_STATUS       | battery_status      |
        | BATTERY_LEVEL        | battery_level       |
        | BATTERY_SCALE        | battery_scale       |

        **MUTATION**

        - **COLUMN_MAPPINGS** (None)
        - **SCRIPTS** (None)

    === "IOS Client V1"

        **RAPIDS_COLUMN_MAPPINGS**

        | RAPIDS column        | Stream column       |
        |----------------------|---------------------|
        | TIMESTAMP            | timestamp           |
        | DEVICE_ID            | device_id           |
        | BATTERY_STATUS       | FLAG_TO_MUTATE      |
        | BATTERY_LEVEL        | battery_level       |
        | BATTERY_SCALE        | battery_scale       |

        **MUTATION**

        - **COLUMN_MAPPINGS**

        | Script column        | Stream column       |
        |----------------------|---------------------|
        | BATTERY_STATUS       | battery_status      |
        
        - **SCRIPTS**

        ```bash
        src/data/streams/mutations/phone/aware/battery_ios_unification.R
        ```

        !!! note
            For RAPIDS `BATTERY_STATUS` column:

            - if stream's `battery_status` field is 3, set `BATTERY_STATUS` = 5 (full status)
            - if stream's `battery_status` field is 1, set `BATTERY_STATUS` = 3 (discharge)
    
    === "IOS Client V2"

        Same as ANDROID


??? info "PHONE_CONVERSATION"

    === "ANDROID"
    
        **RAPIDS_COLUMN_MAPPINGS**

        | RAPIDS column        | Stream column       |
        |----------------------|---------------------|
        | TIMESTAMP            | timestamp           |
        | DEVICE_ID            | device_id           |
        | DOUBLE_ENERGY        | double_energy       |
        | INFERENCE            | inference           |
        | DOUBLE_CONVO_START   | double_convo_start  |
        | DOUBLE_CONVO_END     | double_convo_end    |

        **MUTATION**

        - **COLUMN_MAPPINGS** (None)
        - **SCRIPTS** (None)

    === "IOS"

        **RAPIDS_COLUMN_MAPPINGS**

        | RAPIDS column        | Stream column       |
        |----------------------|---------------------|
        | TIMESTAMP            | timestamp           |
        | DEVICE_ID            | device_id           |
        | DOUBLE_ENERGY        | double_energy       |
        | INFERENCE            | inference           |
        | DOUBLE_CONVO_START   | FLAG_TO_MUTATE      |
        | DOUBLE_CONVO_END     | FLAG_TO_MUTATE      |        

        **MUTATION**

        - **COLUMN_MAPPINGS**

        | Script column        | Stream column       |
        |----------------------|---------------------|
        | DOUBLE_CONVO_START   | double_convo_start  |
        | DOUBLE_CONVO_END     | double_convo_end    |

        - **SCRIPTS**
        
        ```bash
        src/data/streams/mutations/phone/aware/conversation_ios_timestamp.R
        ```

        !!! note
            For RAPIDS columns of `DOUBLE_CONVO_START` and `DOUBLE_CONVO_END`:

            - if stream's `double_convo_start` field is smaller than 9999999999, it is in seconds instead of milliseconds. Set `DOUBLE_CONVO_START` = 1000 * `double_convo_start`.
            - if stream's `double_convo_end` field is smaller than 9999999999, it is in seconds instead of milliseconds. Set `DOUBLE_CONVO_END` = 1000 * `double_convo_end`.
                



