If you collected sensor data with the vanilla (original) AWARE mobile clients, you shouldn't need to modify this format (described below). 

Remember that a format maps and transforms columns in your raw data stream to the [mandatory columns RAPIDS needs](../mandatory-phone-format).

The yaml file that describes the format of this data stream is at:
```bash
src/data/streams/aware_csv/format.yaml
```

For some sensors, we need to transform iOS data into Android format; you can refer to [OS complex mapping](../../datastreams/add-new-data-streams/#os-complex-mapping) for learn how this works.

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


??? info "PHONE_APPLICATIONS_CRASHES"

    === "ANDROID"
    
         **RAPIDS_COLUMN_MAPPINGS**

        | RAPIDS column      | Stream column      |
        |--------------------|--------------------|
        | TIMESTAMP          | timestamp          |
        | DEVICE_ID          | device_id          |
        | PACKAGE_NAME       | package_name       |
        | APPLICATION_NAME   | application_name   |
        | APPLICATION_VERSION| application_version|
        | ERROR_SHORT        | error_short        |
        | ERROR_LONG         | error_long         |
        | ERROR_CONDITION    | error_condition    |
        | IS_SYSTEM_APP      | is_system_app      |

        **MUTATION**

        - **COLUMN_MAPPINGS** (None)
        - **SCRIPTS** (None)
    
    === "IOS"

        This sensor is not supported by iOS devices.


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

??? info "PHONE_APPLICATIONS_NOTIFICATIONS"

    === "ANDROID"
    
         **RAPIDS_COLUMN_MAPPINGS**

        | RAPIDS column      | Stream column      |
        |--------------------|--------------------|
        | TIMESTAMP          | timestamp          |
        | DEVICE_ID          | device_id          |
        | PACKAGE_NAME       | package_name       |
        | APPLICATION_NAME   | application_name   |
        | TEXT               | text               |
        | SOUND              | sound              |
        | VIBRATE            | vibrate            |
        | DEFAULTS           | defaults           |
        | FLAGS              | flags              |

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


??? info "PHONE_BLUETOOTH"

    === "ANDROID"
    
        **RAPIDS_COLUMN_MAPPINGS**

        | RAPIDS column        | Stream column       |
        |----------------------|---------------------|
        | TIMESTAMP            | timestamp           |
        | DEVICE_ID            | device_id           |
        | BT_ADDRESS           | bt_address          |
        | BT_NAME              | bt_name             |
        | BT_RSSI              | bt_rssi             |

        **MUTATION**

        - **COLUMN_MAPPINGS** (None)
        - **SCRIPTS** (None)

    === "IOS"

        Only old iOS versions supported this sensor (same mapping as Android).


??? info "PHONE_CALLS"

    === "ANDROID"
    
        **RAPIDS_COLUMN_MAPPINGS**

        | RAPIDS column        | Stream column       |
        |----------------------|---------------------|
        | TIMESTAMP            | timestamp           |
        | DEVICE_ID            | device_id           |
        | CALL_TYPE            | call_type           |
        | CALL_DURATION        | call_duration       |
        | TRACE                | trace               |

        **MUTATION**

        - **COLUMN_MAPPINGS** (None)
        - **SCRIPTS** (None)

    === "IOS"

        **RAPIDS_COLUMN_MAPPINGS**

        | RAPIDS column        | Stream column       |
        |----------------------|---------------------|
        | TIMESTAMP            | timestamp           |
        | DEVICE_ID            | device_id           |
        | CALL_TYPE            | FLAG_TO_MUTATE      |
        | CALL_DURATION        | call_duration       |
        | TRACE                | trace               |

        **MUTATION**

        - **COLUMN_MAPPINGS**

        | Script column        | Stream column       |
        |----------------------|---------------------|
        | CALL_TYPE            | call_type           |


        - **SCRIPTS**

        ```bash
        src/data/streams/mutations/phone/aware/calls_ios_unification.R
        ```

        !!! note

            We transform iOS call logs into Android's format. iOS stores call status: 1=incoming, 2=connected, 3=dialing, 4=disconnected, as opposed to Android's events: 1=incoming, 2=outgoing, 3=missed. 

            We follow this algorithm to convert iOS call data (there are some inaccuracies in the way we handle sequences, see new rules below):

            - Search for the disconnected (4) status as it is common to all calls
            - Group all events that preceded every status 4
            - We convert every 1,2,4 (or 2,1,4) sequence to an incoming call
            - We convert every 3,2,4 (or 2,3,4) sequence to an outgoing call
            - We convert every 1,4 or 3,4 sequence to a missed call (either incoming or outgoing)
            - We set the duration of the call to be the sum of every status (dialing/ringing to hangup) as opposed to the duration of the last status (pick up to hang up)

            **Tested with an Android (OnePlus 7T) and an iPhone XR**

            |Call type | Android (duration) | iOS (duration) | New Rule|
            |---------|----------|--------|------|
            |Outgoing missed ended by me | 2 (0) | 3,4 (0,X) | 3,4 is converted to 2 with duration 0|
            |Outgoing missed ended by them|2(0)|3,2,4 (0,X,X2)| 3,2,4 is converted to 2 with duration X2*|
            |Incoming missed ended by me|NA**|1,4 (0,X)|1,4 is converted to 3 with duration 0|
            |Incoming missed ended by them|3(0)|1,4 (0,X)|1,4 is converted to 3 with duration 0|
            |Outgoing answered|2(X excluding dialing time)|3,2,4 (0,X,X2)|3,2,4 is converted to 2 with duration X2|
            |Incoming answered|1(X excluding dialing time)|1,2,4 (0,X,X2)|1,2,4 is converted to 1 with duration X2|

            .* There is no way to differentiate an outgoing missed call ended by them from an outgoing answered call because the phone goes directly to voice mail and it counts as call time (essentially the voice mail answered).

            .** Android does not record incoming missed calls ended by the participant, just those ended by the person calling or ignored by the participant.


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


??? info "PHONE_KEYBOARD"

    === "ANDROID"
    
        **RAPIDS_COLUMN_MAPPINGS**

        | RAPIDS column        | Stream column       |
        |----------------------|---------------------|
        | TIMESTAMP            | timestamp           |
        | DEVICE_ID            | device_id           |
        | PACKAGE_NAME         | package_name        |
        | BEFORE_TEXT          | before_text         |
        | CURRENT_TEXT         | current_text        |
        | IS_PASSWORD          | is_password         |

        **MUTATION**

        - **COLUMN_MAPPINGS** (None)
        - **SCRIPTS** (None)

    === "IOS"

        This sensor is not supported by iOS devices.


??? info "PHONE_LIGHT"

    === "ANDROID"
    
        **RAPIDS_COLUMN_MAPPINGS**

        | RAPIDS column        | Stream column       |
        |----------------------|---------------------|
        | TIMESTAMP            | timestamp           |
        | DEVICE_ID            | device_id           |
        | DOUBLE_LIGHT_LUX     | double_light_lux    |
        | ACCURACY             | accuracy            |

        **MUTATION**

        - **COLUMN_MAPPINGS** (None)
        - **SCRIPTS** (None)

    === "IOS"

        This sensor is not supported by iOS devices.


??? info "PHONE_LOCATIONS"

    === "ANDROID"
    
        **RAPIDS_COLUMN_MAPPINGS**

        | RAPIDS column        | Stream column       |
        |----------------------|---------------------|
        | TIMESTAMP            | timestamp           |
        | DEVICE_ID            | device_id           |
        | DOUBLE_LATITUDE      | double_latitude     |
        | DOUBLE_LONGITUDE     | double_longitude    |
        | DOUBLE_BEARING       | double_bearing      |
        | DOUBLE_SPEED         | double_speed        |
        | DOUBLE_ALTITUDE      | double_altitude     |
        | PROVIDER             | provider            |
        | ACCURACY             | accuracy            |

        **MUTATION**

        - **COLUMN_MAPPINGS** (None)
        - **SCRIPTS** (None)

    === "IOS"
    
        Same as ANDROID


??? info "PHONE_LOG"

    === "ANDROID"
    
        **RAPIDS_COLUMN_MAPPINGS**

        | RAPIDS column        | Stream column       |
        |----------------------|---------------------|
        | TIMESTAMP            | timestamp           |
        | DEVICE_ID            | device_id           |
        | LOG_MESSAGE          | log_message         |

        **MUTATION**

        - **COLUMN_MAPPINGS** (None)
        - **SCRIPTS** (None)

    === "IOS"
    
        Same as ANDROID


??? info "PHONE_MESSAGES"

    === "ANDROID"
    
        **RAPIDS_COLUMN_MAPPINGS**

        | RAPIDS column        | Stream column       |
        |----------------------|---------------------|
        | TIMESTAMP            | timestamp           |
        | DEVICE_ID            | device_id           |
        | MESSAGE_TYPE         | message_type        |
        | TRACE                | trace               |

        **MUTATION**

        - **COLUMN_MAPPINGS** (None)
        - **SCRIPTS** (None)

    === "IOS"

        This sensor is not supported by iOS devices.


??? info "PHONE_SCREEN"

    === "ANDROID"
    
        **RAPIDS_COLUMN_MAPPINGS**

        | RAPIDS column        | Stream column       |
        |----------------------|---------------------|
        | TIMESTAMP            | timestamp           |
        | DEVICE_ID            | device_id           |
        | SCREEN_STATUS        | screen_status       |

        **MUTATION**

        - **COLUMN_MAPPINGS** (None)
        - **SCRIPTS** (None)

    === "IOS"

        **RAPIDS_COLUMN_MAPPINGS**

        | RAPIDS column        | Stream column       |
        |----------------------|---------------------|
        | TIMESTAMP            | timestamp           |
        | DEVICE_ID            | device_id           |
        | SCREEN_STATUS        | FLAG_TO_MUTATE      |

        **MUTATION**

        - **COLUMN_MAPPINGS**

        | Script column        | Stream column       |
        |----------------------|---------------------|
        | SCREEN_STATUS        | screen_status       |

        - **SCRIPTS**
        
        ```bash
        src/data/streams/mutations/phone/aware/screen_ios_unification.R
        ```

        !!! note
            For `SCREEN_STATUS` RAPIDS column:
            
            - if stream's `screen_status` field is 2 (lock episode), set `SCREEN_STATUS` = 0 (off episode).


??? info "PHONE_WIFI_CONNECTED"

    === "ANDROID"
    
        **RAPIDS_COLUMN_MAPPINGS**

        | RAPIDS column        | Stream column       |
        |----------------------|---------------------|
        | TIMESTAMP            | timestamp           |
        | DEVICE_ID            | device_id           |
        | MAC_ADDRESS          | mac_address         |
        | SSID                 | ssid                |
        | BSSID                | bssid               |

        **MUTATION**

        - **COLUMN_MAPPINGS** (None)
        - **SCRIPTS** (None)

    === "IOS"
    
        Same as ANDROID


??? info "PHONE_WIFI_VISIBLE"

    === "ANDROID"
    
        **RAPIDS_COLUMN_MAPPINGS**

        | RAPIDS column        | Stream column       |
        |----------------------|---------------------|
        | TIMESTAMP            | timestamp           |
        | DEVICE_ID            | device_id           |
        | SSID                 | ssid                |
        | BSSID                | bssid               |
        | SECURITY             | security            |
        | FREQUENCY            | frequency           |
        | RSSI                 | rssi                |

        **MUTATION**

        - **COLUMN_MAPPINGS** (None)
        - **SCRIPTS** (None)

    === "IOS"

        Only old iOS versions supported this sensor (same mapping as Android).

