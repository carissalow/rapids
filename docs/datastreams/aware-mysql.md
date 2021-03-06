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

Stream columns named `FLAG_TO_MUTATE` means they are extracted from the `FLAG_AS_EXTRA` RAPIDS column. You can refer to [OS complex mapping](../../datastreams/add-new-data-streams/#os-complex-mapping) section for detailed information.

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

??? info "PHONE_ACTIVITY_RECOGNITION"

    === "ANDROID"
    
        **COLUMN_MAPPINGS**

        | RAPIDS column   | Stream column   |
        |-----------------|-----------------|
        | TIMESTAMP       | timestamp       |
        | DEVICE_ID       | device_id       |
        | ACTIVITY_TYPE   | activity_type   |
        | ACTIVITY_NAME   | activity_name   |
        | CONFIDENCE      | confidence      |

        **MUTATION_SCRIPTS**

        None

    === "IOS"

        **COLUMN_MAPPINGS**

        | RAPIDS column   | Stream column   |
        |-----------------|-----------------|
        | TIMESTAMP       | timestamp       |
        | DEVICE_ID       | device_id       |
        | ACTIVITY_TYPE   | FLAG_TO_MUTATE  |
        | ACTIVITY_NAME   | FLAG_TO_MUTATE  |
        | CONFIDENCE      | confidence      |
        | FLAG_AS_EXTRA   | activities      |

        **MUTATION_SCRIPTS**
        
        ```bash
        src/data/streams/mutations/phone/aware/activity_recogniton_ios_unification.R
        ```

??? info "PHONE_CONVERSATION"

    === "ANDROID"
    
        **COLUMN_MAPPINGS**

        | RAPIDS column        | Stream column       |
        |----------------------|---------------------|
        | TIMESTAMP            | timestamp           |
        | DEVICE_ID            | device_id           |
        | DOUBLE_ENERGY        | double_energy       |
        | INFERENCE            | inference           |
        | DOUBLE_CONVO_START   | double_convo_start  |
        | DOUBLE_CONVO_END     | double_convo_end    |

        **MUTATION_SCRIPTS**

        None

    === "IOS"

        **COLUMN_MAPPINGS**

        Same as ANDROID

        **MUTATION_SCRIPTS**
        
        ```bash
        src/data/streams/mutations/phone/aware/conversation_ios_timestamp.R
        ```
