# Test Cases

Along with the continued development and the addition of new sensors and features to the RAPIDS pipeline, tests for the currently available sensors and features are being implemented. Since this is a Work In Progress this page will be updated with the list of sensors and features for which testing is available. For each of the sensors listed a description of the data used for testing (test cases) are outline. Currently for all intent and testing purposes the `tests/data/raw/test01/` contains all the test data files for testing android data formats and `tests/data/raw/test02/` contains all the test data files for testing iOS data formats. It follows that the expected (verified output) are contained in the `tests/data/processed/test01/` and `tests/data/processed/test02/` for Android and iOS respectively. `tests/data/raw/test03/` and `tests/data/raw/test04/` contain data files for testing empty raw data files for android and iOS respectively.

The following is a list of the sensors that testing is currently available.


| Sensor                        | Provider | Periodic | Frequency | Event |
|-------------------------------|----------|----------|-----------|-------|
| Phone Accelerometer           | Panda    | N        | N         | N     |
| Phone Accelerometer           | RAPIDS   | N        | N         | N     |
| Phone Activity Recognition    | RAPIDS   | N        | N         | N     |
| Phone Applications Foreground | RAPIDS   | N        | N         | N     |
| Phone Battery                 | RAPIDS   | Y        | Y         | N     |
| Phone Bluetooth               | Doryab   | N        | N         | N     |
| Phone Bluetooth               | RAPIDS   | Y        | Y         | Y     |
| Phone Calls                   | RAPIDS   | Y        | Y         | Y     |
| Phone Conversation            | RAPIDS   | Y        | Y         | N     |
| Phone Data Yield              | RAPIDS   | N        | N         | N     |
| Phone Light                   | RAPIDS   | Y        | Y         | N     |
| Phone Locations               | Doryab   | N        | N         | N     |
| Phone Locations               | Barnett  | N        | N         | N     |
| Phone Messages                | RAPIDS   | Y        | Y         | N     |
| Phone Screen                  | RAPIDS   | Y        | Y         | Y     |
| Phone WiFi Connected          | RAPIDS   | Y        | Y         | Y     |
| Phone WiFi Visible            | RAPIDS   | Y        | Y         | Y     |
| Fitbit Calories Intraday      | RAPIDS   | Y        | Y         | Y     |
| Fitbit Data Yield             | RAPIDS   | N        | N         | N     |
| Fitbit Heart Rate Summary     | RAPIDS   | N        | N         | N     |
| Fitbit Heart Rate Intraday    | RAPIDS   | N        | N         | N     |
| Fitbit Sleep Summary          | RAPIDS   | N        | N         | N     |
| Fitbit Sleep Intraday         | RAPIDS   | Y        | Y         | Y     |
| Fitbit Sleep Intraday         | PRICE    | Y        | Y         | Y     |
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

Due to the difference in the format of the raw data for iOS and Android the following is the expected results 
the `phone_calls.csv`. 

Description:

- One missed episode, one outgoing episode and one incoming episode on Friday night, morning, afternoon and evening
- There is at least one episode of each type of phone calls on each day
- One incoming episode crossing two 30-mins segments
- One outgoing episode crossing two 30-mins segments 
- One missed episode before, during and after the `event`
- There is one incoming episode before, during or after the `event`
- There is one outcoming episode before, during or after the `event`
- There is one missed episode before, during or after the `event`

Data format:

| Device | Missed | Outgoing | Incoming |
|-|-|-|-|
|iOS| 3 | 2 | 1 |
|Android| 1,4 or 3,4 | 3,2,4 | 1,2,4 |

Note:
When generating test data, all traces for iOS device need to be unique otherwise the episode with duplicate trace will be dropped 

Checklist:

|time segment| single tz | multi tz|platform|
|-|-|-|-|
|30min|OK|OK|Android, iOS|
|morning|OK|OK|Android, iOS|
|daily|OK|OK|Android, iOS|
|threeday|OK|OK|Android, iOS|
|weekend|OK|OK|Android, iOS|
|beforeMarchEvent|OK|OK|Android, iOS|
|beforeNovemberEvent|OK|OK|Android, iOS|

## Screen

Due to the difference in the format of the raw screen data for iOS and Android the following is the expected results the `phone_screen.csv`. 

Description:

- The screen data file contains data for 4 days.
- The screen data contains 1 record to represent an `unlock`
    episode that falls within an `epoch` for every `epoch`.
- The screen data contains 1 record to represent an `unlock`
    episode that falls across the boundary of 2 epochs. Namely the
    `unlock` episode starts in one epoch and ends in the next, thus
    there is a record for `unlock` episodes that fall across `night`
    to `morning`, `morning` to `afternoon` and finally `afternoon` to
    `night`
- One episode that crossing two `30-min` segments

Data format:

| Device | unlock |
|-|-|
| Android | 3, 0|
| iOS | 3, 2|


Checklist

|time segment| single tz | multi tz|platform|
|-|-|-|-|
|30min|OK|OK|Android, iOS|
|morning|OK|OK|Android, iOS|
|daily|OK|OK|Android, iOS|
|threeday|OK|OK|Android, iOS|
|weekend|OK|OK|Android, iOS|
|beforeMarchEvent|OK|OK|Android, iOS|
|beforeNovemberEvent|OK|OK|Android, iOS|

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

There are two wifi features (`phone wifi connected` and `phone wifi visible`). The raw test data are seperatly stored in the `phone_wifi_connected_raw.csv` and `phone_wifi_visible_raw.csv`.

Description: 

- One episode for each `epoch` (`night`, `morining`, `afternoon` and `evening`)
- Two two episodes in the same time segment (`daily` and `30-min`)
- Two episodes around the transition of `epochs` (e.g. one at the end of `night` and one at the beginning of `morning`) 
- One episode before and after the time switch on Sunday

phone wifi connected:

Checklist

|time segment| single tz | multi tz|platform|
|-|-|-|-|
|30min|OK|OK|android, iOS|
|morning|OK|OK|android, iOS|
|daily|OK|OK|android, iOS|
|threeday|OK|OK|android, iOS|
|weekend|OK|OK|android, iOS|
|beforeMarchEvent|OK|OK|android, iOS|
|beforeNovemberEvent|OK|OK|android, iOS|

phone wifi visible

Checklist

|time segment| single tz | multi tz|platform|
|-|-|-|-|
|30min|OK|OK|android|
|morning|OK|OK|android|
|daily|OK|OK|android|
|threeday|OK|OK|android|
|weekend|OK|OK|android|
|beforeMarchEvent|OK|OK|android|
|beforeNovemberEvent|OK|OK|android|

## Light

-   The raw light data file contains data for 1 day.
-   The raw light data contains 3 or 4 rows of data for each `epoch`
    except `night`. The single row of data for `night` is for testing
    features for single values inputs. (Example testing the standard
    deviation of one input value)
-   Since light is only available for Android there is only one file
    that contains data for Android. All other files (i.e. for iPhone)
    are empty data files.

## Locations

Description

- The participant's home location is (latitude=1, longitude=1).
- From Sat 10:56:00 to Sat 11:04:00, the center of the cluster is (latitude=-100, longitude=-100).
- From Sun 03:30:00 to Sun 03:47:00, the center of the cluster is (latitude=1, longitude=1). Home location is extracted from this period.
- From Sun 11:30:00 to Sun 11:38:00, the center of the cluster is (latitude=100, longitude=100).

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

## Keyboard

- The raw keyboard data file contains data for 4 days.
- The raw keyboard data contains records with difference in `timestamp` ranging from
  milliseconds to seconds.  
  
- With difference in timestamps between consecutive records more than 5 seconds helps us to create separate 
  sessions within the usage of the same app. This helps to verify the case where sessions have to be different. 

- The raw keyboard data contains records where the difference in text is less 
  than 5 seconds which makes it into 1 session but because of difference of app
  new session starts. This edge case determines the behaviour within particular app
  and also within 5 seconds.

- The raw keyboard data also contains the records where length of `current_text` varies between consecutive rows. This helps us to tests on the cases where input text is entered by auto-suggested
  or auto-correct operations.

- One three-minute episode with a 1-minute row on Sun 08:59:54.65 and 09:00:00,another on Sun 12:01:02 that are considering a single episode in multi-timezone event segments to showcase how
 inferring time zone data for Keyboard from phone data can produce inaccurate results around the tz change. This happens because the device was on LA time until 11:59 and switched to NY time at 12pm, in terms of actual time 09 am LA and 12 pm NY represent the same moment in time so 09:00 LA and 12:01 NY are consecutive minutes.

## Fitbit Calories Intraday

Description

- A five-minute sedentary episode on Fri 11:00:00
- A one-minute sedentary episode on Sun 02:00:00. It exists in November but not in February in STZ
- A five-minute sedentary episode on Fri 11:58:00. It is split within two 30-min segments and the morning
- A three-minute lightly active episode on Fri 11:10:00, a one-minute at 11:18:00 and a one-minute 11:24:00. These check for start and end times of first/last/longest episode
- A three-minute fairly active episode on Fri 11:40:00, a one-minute at 11:48:00 and a one-minute 11:54:00. These check for start and end times of first/last/longest episode
- A three-minute very active episode on Fri 12:10:00, a one-minute at 12:18:00 and a one-minute 12:24:00. These check for start and end times of first/last/longest episode
- A eight-minute MVPA episode with intertwined fairly and very active rows on Fri 12:30:00
- The above episodes contain six higmet (>= 3 MET) episodes and nine lowmet episodes.
- One two-minute sedentary episode with a 1-minute row on Sun 09:00:00 and another on Sun 12:01:01 that are considering a single episode in multi-timezone event segments to showcase how inferring time zone data for Fitbit from phone data can produce inaccurate results around the tz change. This happens because the device was on LA time until 11:59 and switched to NY time at 12pm, in terms of actual time 09 am LA and 12 pm NY represent the same moment in time so 09:00 LA and 12:01 NY are consecutive minutes.
- A three-minute sedentary episode on Sat 08:59 that will be ignored for multi-timezone event segments.
- A three-minute sedentary episode on Sat 12:59 of which the first minute will be ignored for multi-timezone event segments since the test segment starts at 13:00
- A three-minute sedentary episode on Sat 16:00
- A four-minute sedentary episode on Sun 10:01 that will be ignored for Novembers's multi-timezone event segments since the test segment ends at 10am on that weekend.
- A three-minute very active episode on Sat 16:03. This episode and the one at 16:00 are counted as one for lowmet episodes

Checklist

|time segment| single tz | multi tz|platform|
|-|-|-|-|
|30min|OK|OK|fitbit|
|morning|OK|OK|fitbit|
|daily|OK|OK|fitbit|
|threeday|OK|OK|fitbit|
|weekend|OK|OK|fitbit|
|beforeMarchEvent|OK|OK|fitbit|
|beforeNovemberEvent|OK|OK|fitbit|

## Fitbit Sleep Summary

Description

- A main sleep episode that starts on Fri 20:00:00 and ends on Sat 02:00:00. This episode starts after 11am (Last Night End) which will be considered as today's (Fri) data.
- A nap that starts on Sat 04:00:00 and ends on Sat 06:00:00. This episode starts before 11am (Last Night End) which will be considered as yesterday's (Fri) data.
- A nap that starts on Sat 13:00:00 and ends on Sat 15:00:00. This episode starts after 11am (Last Night End) which will be considered as today's (Sat) data.
- A main sleep that starts on Sun 01:00:00 and ends on Sun 12:00:00. This episode starts before 11am (Last Night End) which will be considered as yesterday's (Sat) data.
- A main sleep that starts on Sun 23:00:00 and ends on Mon 07:00:00. This episode starts after 11am (Last Night End) which will be considered as today's (Sun) data.
- Any segment shorter than one day will be ignored for sleep RAPIDS features.

Checklist

|time segment| single tz | multi tz|platform|
|-|-|-|-|
|30min|OK|OK|fitbit|
|morning|OK|OK|fitbit|
|daily|OK|OK|fitbit|
|threeday|OK|OK|fitbit|
|weekend|OK|OK|fitbit|
|beforeMarchEvent|OK|OK|fitbit|
|beforeNovemberEvent|OK|OK|fitbit|

## Fitbit Sleep Intraday

Description

- A five-minute main sleep episode with asleep-classic level on Fri 11:00:00.
- An eight-hour main sleep episode on Fri 17:00:00. It is split into 2 parts for daily segment: a seven-hour sleep episode on Fri 17:00:00 and an one-hour sleep episode on Sat 00:00:00.
- A two-hour nap on Sat 01:00:00 that will be ignored for main sleep features.
- An one-hour nap on Sat 13:00:00 that will be ignored for main sleep features.
- An eight-hour main sleep episode on Sat 22:00:00. This episode ends on Sun 08:00:00 (NY) for March and Sun 06:00:00 (NY) for Novembers due to daylight savings. It will be considered for `beforeMarchEvent` segment and ignored for `beforeNovemberEvent` segment.
- A nine-hour main sleep episode on Sun 11:00:00. Start time will be assigned as NY time zone and converted to 14:00:00.
- A seven-hour main sleep episode on Mon 06:00:00. This episode will be split into two parts: a five-hour sleep episode on Mon 06:00:00 and a two-hour sleep episode on Mon 11:00:00. The first part will be discarded as it is before 11am (Last Night End)
- Any segment shorter than one day will be ignored for sleep PRICE features.

Checklist

|time segment| single tz | multi tz|platform|
|-|-|-|-|
|30min|OK|OK|fitbit|
|morning|OK|OK|fitbit|
|daily|OK|OK|fitbit|
|threeday|OK|OK|fitbit|
|weekend|OK|OK|fitbit|
|beforeMarchEvent|OK|OK|fitbit|
|beforeNovemberEvent|OK|OK|fitbit|
