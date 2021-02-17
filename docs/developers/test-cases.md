# Test Cases

Along with the continued development and the addition of new sensors and features to the RAPIDS pipeline, tests for the currently available sensors and features are being implemented. Since this is a Work In Progress this page will be updated with the list of sensors and features for which testing is available. For each of the sensors listed a description of the data used for testing (test cases) are outline. Currently for all intent and testing purposes the `tests/data/raw/test01/` contains all the test data files for testing android data formats and `tests/data/raw/test02/` contains all the test data files for testing iOS data formats. It follows that the expected (verified output) are contained in the `tests/data/processed/test01/` and `tests/data/processed/test02/` for Android and iOS respectively. `tests/data/raw/test03/` and `tests/data/raw/test04/` contain data files for testing empty raw data files for android and iOS respectively.

The following is a list of the sensors that testing is currently available.


| Sensor                        | Provider | Periodic | Frequency | Event |
|-------------------------------|----------|----------|-----------|-------|
| Phone Accelerometer           | Panda    | N        | N         | N     |
| Phone Accelerometer           | RAPIDS   | N        | N         | N     |
| Phone Activity Recognition    | RAPIDS   | N        | N         | N     |
| Phone Applications Foreground | RAPIDS   | N        | N         | N     |
| Phone Battery                 | RAPIDS   | Y        | N         | N     |
| Phone Bluetooth               | Doryab   | N        | N         | N     |
| Phone Bluetooth               | RAPIDS   | Y        | Y         | Y     |
| Phone Calls                   | RAPIDS   | Y        | Y         | N     |
| Phone Conversation            | RAPIDS   | Y        | Y         | N     |
| Phone Data Yield              | RAPIDS   | N        | N         | N     |
| Phone Light                   | RAPIDS   | Y        | Y         | N     |
| Phone Locations               | Doryab   | N        | N         | N     |
| Phone Locations               | Barnett  | N        | N         | N     |
| Phone Messages                | RAPIDS   | Y        | Y         | N     |
| Phone Screen                  | RAPIDS   | N        | N         | N     |
| Phone WiFi Connected          | RAPIDS   | Y        | Y         | N     |
| Phone WiFi Visible            | RAPIDS   | Y        | Y         | N     |
| Fitbit Data Yield             | RAPIDS   | N        | N         | N     |
| Fitbit Heart Rate Summary     | RAPIDS   | N        | N         | N     |
| Fitbit Heart Rate Intraday    | RAPIDS   | N        | N         | N     |
| Fitbit Sleep Summary          | RAPIDS   | N        | N         | N     |
| Fitbit Steps Summary          | RAPIDS   | N        | N         | N     |
| Fitbit Steps Intraday         | RAPIDS   | N        | N         | N     |


## Messages (SMS)

-   The raw message data file contains data for 2 separate days.
-   The data for the first day contains records 5 records for every
    `epoch`.
-   The second day\'s data contains 6 records for each of only 2
    `epoch` (currently `morning` and `evening`)
-   The raw message data contains records for both `message_types`
    (i.e. `recieved` and `sent`) in both days in all epochs. The
    number records with each `message_types` per epoch is randomly
    distributed There is at least one records with each
    `message_types` per epoch.
-   There is one raw message data file each, as described above, for
    testing both iOS and Android data.
-   There is also an additional empty data file for both android and
    iOS for testing empty data files

## Calls

Due to the difference in the format of the raw call data for iOS and Android the following is the expected results the `calls_with_datetime_unified.csv`. This would give a better idea of the use cases being tested since the `calls_with_datetime_unified.csv` would make both the iOS and Android data comparable.

-   The call data would contain data for 2 days.
-   The data for the first day contains 6 records for every `epoch`.
-   The second day\'s data contains 6 records for each of only 2
    `epoch` (currently `morning` and `evening`)
-   The call data contains records for all `call_types` (i.e.
    `incoming`, `outgoing` and `missed`) in both days in all epochs.
    The number records with each of the `call_types` per epoch is
    randomly distributed. There is at least one records with each
    `call_types` per epoch.
-   There is one call data file each, as described above, for testing
    both iOS and Android data.
-   There is also an additional empty data file for both android and
    iOS for testing empty data files

## Screen

Due to the difference in the format of the raw screen data for iOS and Android the following is the expected results the `screen_deltas.csv`. This would give a better idea of the use cases being tested since the `screen_eltas.csv` would make both the iOS and Android data comparable These files are used to calculate the features for the screen sensor

-   The screen delta data file contains data for 1 day.
-   The screen delta data contains 1 record to represent an `unlock`
    episode that falls within an `epoch` for every `epoch`.
-   The screen delta data contains 1 record to represent an `unlock`
    episode that falls across the boundary of 2 epochs. Namely the
    `unlock` episode starts in one epoch and ends in the next, thus
    there is a record for `unlock` episodes that fall across `night`
    to `morning`, `morning` to `afternoon` and finally `afternoon` to
    `night`
-   The testing is done for `unlock` episode\_type.
-   There is one screen data file each for testing both iOS and
    Android data formats.
-   There is also an additional empty data file for both android and
    iOS for testing empty data files

## Battery

Due to the difference in the format of the raw battery data for iOS and Android as well as versions of iOS the following is the expected results the `battery_deltas.csv`. This would give a better idea of the use cases being tested since the `battery_deltas.csv` would make both the iOS and Android data comparable. These files are used to calculate the features for the battery sensor.

-   The battery delta data file contains data for 1 day.
-   The battery delta data contains 1 record each for a `charging` and
    `discharging` episode that falls within an `epoch` for every
    `epoch`. Thus, for the `daily` epoch there would be multiple
    `charging` and `discharging` episodes
-   Since either a `charging` episode or a `discharging` episode and
    not both can occur across epochs, in order to test episodes that
    occur across epochs alternating episodes of `charging` and
    `discharging` episodes that fall across `night` to `morning`,
    `morning` to `afternoon` and finally `afternoon` to `night` are
    present in the battery delta data. This starts with a
    `discharging` episode that begins in `night` and end in `morning`.
-   There is one battery data file each, for testing both iOS and
    Android data formats.
-   There is also an additional empty data file for both android and
    iOS for testing empty data files

## Bluetooth

-   The raw Bluetooth data file contains data for 1 day.
-   The raw Bluetooth data contains at least 2 records for each
    `epoch`. Each `epoch` has a record with a `timestamp` for the
    beginning boundary for that `epoch` and a record with a
    `timestamp` for the ending boundary for that `epoch`. (e.g. For
    the `morning` epoch there is a record with a `timestamp` for
    `6:00AM` and another record with a `timestamp` for `11:59:59AM`.
    These are to test edge cases)
-   An option of 5 Bluetooth devices are randomly distributed
    throughout the data records.
-   There is one raw Bluetooth data file each, for testing both iOS
    and Android data formats.
-   There is also an additional empty data file for both android and
    iOS for testing empty data files.

## WIFI

-   There are 2 data files (`wifi_raw.csv` and `sensor_wifi_raw.csv`)
    for each fake participant for each phone platform. 
-   The raw WIFI data files contain data for 1 day.
-   The `sensor_wifi_raw.csv` data contains at least 2 records for
    each `epoch`. Each `epoch` has a record with a `timestamp` for the
    beginning boundary for that `epoch` and a record with a
    `timestamp` for the ending boundary for that `epoch`. (e.g. For
    the `morning` epoch there is a record with a `timestamp` for
    `6:00AM` and another record with a `timestamp` for `11:59:59AM`.
    These are to test edge cases)
-   The `wifi_raw.csv` data contains 3 records with random timestamps
    for each `epoch` to represent visible broadcasting WIFI network.
    This file is empty for the iOS phone testing data.
-   An option of 10 access point devices is randomly distributed
    throughout the data records. 5 each for `sensor_wifi_raw.csv` and
    `wifi_raw.csv`.
-   There data files for testing both iOS and Android data formats.
-   There are also additional empty data files for both android and
    iOS for testing empty data files.

## Light

-   The raw light data file contains data for 1 day.
-   The raw light data contains 3 or 4 rows of data for each `epoch`
    except `night`. The single row of data for `night` is for testing
    features for single values inputs. (Example testing the standard
    deviation of one input value)
-   Since light is only available for Android there is only one file
    that contains data for Android. All other files (i.e. for iPhone)
    are empty data files.

## Application Foreground

-   The raw application foreground data file contains data for 1 day.
-   The raw application foreground data contains 7 - 9 rows of data
    for each `epoch`. The records for each `epoch` contains apps that
    are randomly selected from a list of apps that are from the
    `MULTIPLE_CATEGORIES` and `SINGLE_CATEGORIES` (See
    [testing\_config.yaml]()). There are also records in each epoch
    that have apps randomly selected from a list of apps that are from
    the `EXCLUDED_CATEGORIES` and `EXCLUDED_APPS`. This is to test
    that these apps are actually being excluded from the calculations
    of features. There are also records to test `SINGLE_APPS`
    calculations.
-   Since application foreground is only available for Android there
    is only one file that contains data for Android. All other files
    (i.e. for iPhone) are empty data files.

## Activity Recognition

-   The raw Activity Recognition data file contains data for 1 day.
-   The raw Activity Recognition data each `epoch` period contains
    rows that records 2 - 5 different `activity_types`. The is such
    that durations of activities can be tested. Additionally, there
    are records that mimic the duration of an activity over the time
    boundary of neighboring epochs. (For example, there a set of
    records that mimic the participant `in_vehicle` from `afternoon`
    into `evening`)
-   There is one file each with raw Activity Recognition data for
    testing both iOS and Android data formats.
    (plugin\_google\_activity\_recognition\_raw.csv for android and
    plugin\_ios\_activity\_recognition\_raw.csv for iOS)
-   There is also an additional empty data file for both android and
    iOS for testing empty data files.

## Conversation

-   The raw conversation data file contains data for 2 day.
-   The raw conversation data contains records with a sample of both
    `datatypes` (i.e. `voice/noise` = `0`, and `conversation` = `2` )
    as well as rows with for samples of each of the `inference` values
    (i.e. `silence` = `0`, `noise` = `1`, `voice` = `2`, and `unknown`
    = `3`) for each `epoch`. The different `datatype` and `inference`
    records are randomly distributed throughout the `epoch`.
-   Additionally there are 2 - 5 records for conversations (`datatype`
    = 2, and `inference` = -1) in each `epoch` and for each `epoch`
    except night, there is a conversation record that has a
    `double_convo_start` `timestamp` that is from the previous
    `epoch`. This is to test the calculations of features across
    `epochs`.
-   There is a raw conversation data file for both android and iOS
    platforms (`plugin_studentlife_audio_android_raw.csv` and
    `plugin_studentlife_audio_raw.csv` respectively).
-   Finally, there are also additional empty data files for both
    android and iOS for testing empty data files
