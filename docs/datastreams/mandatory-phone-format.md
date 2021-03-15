# Mandatory Phone Format

This is a description of the format RAPIDS needs to process data for the following PHONE sensors.

See examples in the CSV files inside [rapids_example_csv.zip](https://osf.io/wbg23/)

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
    | ACTIVITY_NAME   | An string that denotes current activity name: `in_vehicle`, `on_bicycle`, `on_foot`, `still`, `unknown`, `tilting`, `walking` or `running`   |
    | ACTIVITY_TYPE   | An integer (ranged from 0 to 8) that denotes current activity type        |
    | CONFIDENCE      | An integer (ranged from 0 to 100) that denotes the prediction accuracy    |


??? info "PHONE_APPLICATIONS_CRASHES"

    | RAPIDS column      | Description                                                               |
    |--------------------|---------------------------------------------------------------------------|
    | TIMESTAMP          | An UNIX timestamp (13 digits) when a row of data was logged               |
    | DEVICE_ID          | A string that uniquely identifies a device                                |
    | PACKAGE_NAME       | Application’s package name                                                |
    | APPLICATION_NAME   | Application’s localized name                                              |
    | APPLICATION_VERSION| Application’s version code                                                |
    | ERROR_SHORT        | Short description of the error                                            |
    | ERROR_LONG         | More verbose version of the error description                             |
    | ERROR_CONDITION    | 1 = code error; 2 = non-responsive (ANR error)                            |
    | IS_SYSTEM_APP      | Device’s pre-installed application                                        |


??? info "PHONE_APPLICATIONS_FOREGROUND"

    | RAPIDS column      | Description                                                               |
    |--------------------|---------------------------------------------------------------------------|
    | TIMESTAMP          | An UNIX timestamp (13 digits) when a row of data was logged               |
    | DEVICE_ID          | A string that uniquely identifies a device                                |
    | PACKAGE_NAME       | Application’s package name                                                |
    | APPLICATION_NAME   | Application’s localized name                                              |
    | IS_SYSTEM_APP      | Device’s pre-installed application                                        |


??? info "PHONE_APPLICATIONS_NOTIFICATIONS"

    | RAPIDS column      | Description                                                               |
    |--------------------|---------------------------------------------------------------------------|
    | TIMESTAMP          | An UNIX timestamp (13 digits) when a row of data was logged               |
    | DEVICE_ID          | A string that uniquely identifies a device                                |
    | PACKAGE_NAME       | Application’s package name                                                |
    | APPLICATION_NAME   | Application’s localized name                                              |
    | TEXT               | Notification’s header text, not the content                               |
    | SOUND              | Notification’s sound source (if applicable)                               |
    | VIBRATE            | Notification’s vibration pattern (if applicable)                          |
    | DEFAULTS           | If notification was delivered according to device’s default settings      |
    | FLAGS              | An integer that denotes [Android notification flag](https://developer.android.com/reference/android/app/Notification.html)  |


??? info "PHONE_BATTERY"

    | RAPIDS column        | Description                                                                                                            |
    |----------------------|------------------------------------------------------------------------------------------------------------------------|
    | TIMESTAMP            | An UNIX timestamp (13 digits) when a row of data was logged                                                            |
    | DEVICE_ID            | A string that uniquely identifies a device                                                                             |
    | BATTERY_STATUS       | An integer that denotes battery status: 0 or 1 = unknown, 2 = charging, 3 = discharging, 4 = not charging, 5 = full    |
    | BATTERY_LEVEL        | An integer that denotes battery level, between 0 and `BATTERY_SCALE`                                                   |
    | BATTERY_SCALE        | An integer that denotes the maximum battery level                                                                      |


??? info "PHONE_BLUETOOTH"

    | RAPIDS column      | Description                                                               |
    |--------------------|---------------------------------------------------------------------------|
    | TIMESTAMP          | An UNIX timestamp (13 digits) when a row of data was logged               |
    | DEVICE_ID          | A string that uniquely identifies a device                                |
    | BT_ADDRESS         | MAC address of the device’s Bluetooth sensor                              |
    | BT_NAME            | User assigned name of the device’s Bluetooth sensor                       |
    | BT_RSSI            | The RSSI dB to the scanned device                                         |


??? info "PHONE_CALLS"

    | RAPIDS column      | Description                                                               |
    |--------------------|---------------------------------------------------------------------------|
    | TIMESTAMP          | An UNIX timestamp (13 digits) when a row of data was logged               |
    | DEVICE_ID          | A string that uniquely identifies a device                                |
    | CALL_TYPE          | An integer that denotes call type: 1 = incoming, 2 = outgoing, 3 = missed |
    | CALL_DURATION      | Length of the call session                                                |
    | TRACE              | SHA-1 one-way source/target of the call                                   |


??? info "PHONE_CONVERSATION"

    | RAPIDS column        | Description                                                                          |
    |----------------------|--------------------------------------------------------------------------------------|
    | TIMESTAMP            | An UNIX timestamp (13 digits) when a row of data was logged                          |
    | DEVICE_ID            | A string that uniquely identifies a device                                           |
    | DOUBLE_ENERGY        | A number that denotes the amplitude of an audio sample (L2-norm of the audio frame)     |
    | INFERENCE            | An integer (ranged from 0 to 3) that denotes the type of an audio sample: 0 = silence, 1 = noise, 2 = voice, 3 = unknown      |
    | DOUBLE_CONVO_START   | UNIX timestamp (13 digits) of the beginning of a conversation                        |
    | DOUBLE_CONVO_END     | UNIX timestamp (13 digits) of the end of a conversation                              |


??? info "PHONE_KEYBOARD"

    | RAPIDS column      | Description                                                               |
    |--------------------|---------------------------------------------------------------------------|
    | TIMESTAMP          | An UNIX timestamp (13 digits) when a row of data was logged               |
    | DEVICE_ID          | A string that uniquely identifies a device                                |
    | PACKAGE_NAME       | The application’s package name of keyboard interaction                    |
    | BEFORE_TEXT        | The previous keyboard input (empty if password)                           |
    | CURRENT_TEXT       | The current keyboard input (empty if password)                            |
    | IS_PASSWORD        | An integer: 0 = not password; 1 = password                                |


??? info "PHONE_LIGHT"

    | RAPIDS column      | Description                                                                                                          |
    |--------------------|----------------------------------------------------------------------------------------------------------------------|
    | TIMESTAMP          | An UNIX timestamp (13 digits) when a row of data was logged                                                          |
    | DEVICE_ID          | A string that uniquely identifies a device                                                                           |
    | DOUBLE_LIGHT_LUX   | The ambient luminance in lux units                                                                                   |
    | ACCURACY           | An integer that denotes the sensor's accuracy level: 3 = maximum accuracy, 2 = medium accuracy, 1 = low accuracy     |


??? info "PHONE_LOCATIONS"

    | RAPIDS column      | Description                                                               |
    |--------------------|---------------------------------------------------------------------------|
    | TIMESTAMP          | An UNIX timestamp (13 digits) when a row of data was logged               |
    | DEVICE_ID          | A string that uniquely identifies a device                                |
    | DOUBLE_LATITUDE    | The location’s latitude, in degrees                                       |
    | DOUBLE_LONGITUDE   | The location’s longitude, in degrees                                      |
    | DOUBLE_BEARING     | The location’s bearing, in degrees                                        |
    | DOUBLE_SPEED       | The speed if available, in meters/second over ground                      |
    | DOUBLE_ALTITUDE    | The altitude if available, in meters above sea level                      |
    | PROVIDER           | A string that denotes the provider: `gps`, `fused` or `network`           |
    | ACCURACY           | The estimated location accuracy                                           |


??? info "PHONE_LOG"

    | RAPIDS column      | Description                                                               |
    |--------------------|---------------------------------------------------------------------------|
    | TIMESTAMP          | An UNIX timestamp (13 digits) when a row of data was logged               |
    | DEVICE_ID          | A string that uniquely identifies a device                                |
    | LOG_MESSAGE        | A string that denotes log message                                         |


??? info "PHONE_MESSAGES"

    | RAPIDS column      | Description                                                               |
    |--------------------|---------------------------------------------------------------------------|
    | TIMESTAMP          | An UNIX timestamp (13 digits) when a row of data was logged               |
    | DEVICE_ID          | A string that uniquely identifies a device                                |
    | MESSAGE_TYPE       | An integer that denotes message type: 1 = received, 2 = sent              |
    | TRACE              | SHA-1 one-way source/target of the message                                |


??? info "PHONE_SCREEN"

    | RAPIDS column      | Description                                                                       |
    |--------------------|-----------------------------------------------------------------------------------|
    | TIMESTAMP          | An UNIX timestamp (13 digits) when a row of data was logged                       |
    | DEVICE_ID          | A string that uniquely identifies a device                                        |
    | SCREEN_STATUS      | An integer that denotes screen status: 0 = off, 1 = on, 2 = locked, 3 = unlocked  |


??? info "PHONE_WIFI_CONNECTED"

    | RAPIDS column      | Description                                                                       |
    |--------------------|-----------------------------------------------------------------------------------|
    | TIMESTAMP          | An UNIX timestamp (13 digits) when a row of data was logged                       |
    | DEVICE_ID          | A string that uniquely identifies a device                                        |
    | MAC_ADDRESS        | Device’s MAC address                                                              |
    | SSID               | Currently connected access point network name                                     |
    | BSSID              | Currently connected access point MAC address                                      |


??? info "PHONE_WIFI_VISIBLE"

    | RAPIDS column      | Description                                                                       |
    |--------------------|-----------------------------------------------------------------------------------|
    | TIMESTAMP          | An UNIX timestamp (13 digits) when a row of data was logged                       |
    | DEVICE_ID          | A string that uniquely identifies a device                                        |
    | SSID               | Detected access point network name                                                |
    | BSSID              | Detected access point MAC address                                                 |
    | SECURITY           | Active security protocols                                                         |
    | FREQUENCY          | Wi-Fi band frequency (e.g., 2427, 5180), in Hz                                    |
    | RSSI               | RSSI dB to the scanned device                                                     |

