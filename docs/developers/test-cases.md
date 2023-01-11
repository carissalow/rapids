# Test Cases

Along with the continued development and the addition of new sensors and features to the RAPIDS pipeline, tests for the currently available sensors and features are being implemented. Since this is a Work In Progress this page will be updated with the list of sensors and features for which testing is available. For each of the sensors listed a description of the data used for testing (test cases) are outline. Currently for all intent and testing purposes the `tests/data/raw/test01/` contains all the test data files for testing android data formats and `tests/data/raw/test02/` contains all the test data files for testing iOS data formats. It follows that the expected (verified output) are contained in the `tests/data/processed/test01/` and `tests/data/processed/test02/` for Android and iOS respectively. `tests/data/raw/test03/` and `tests/data/raw/test04/` contain data files for testing empty raw data files for android and iOS respectively.

The following is a list of the sensors that testing is currently available.


| Sensor                        | Provider | Periodic | Frequency | Event |
|-------------------------------|----------|----------|-----------|-------|
| Phone Accelerometer           | Panda    | Y        | Y         | Y     |
| Phone Accelerometer           | RAPIDS   | Y        | Y         | Y     |
| Phone Activity Recognition    | RAPIDS   | Y        | Y         | Y     |
| Phone Applications Foreground | RAPIDS   | Y        | Y         | Y     |
| Phone Battery                 | RAPIDS   | Y        | Y         | Y     |
| Phone Bluetooth               | Doryab   | Y        | Y         | Y     |
| Phone Bluetooth               | RAPIDS   | Y        | Y         | Y     |
| Phone Calls                   | RAPIDS   | Y        | Y         | Y     |
| Phone Conversation            | RAPIDS   | Y        | Y         | Y     |
| Phone Data Yield              | RAPIDS   | Y        | Y         | Y     |
| Phone Light                   | RAPIDS   | Y        | Y         | Y     |
| Phone Locations               | Doryab   | Y        | Y         | Y     |
| Phone Locations               | Barnett  | N        | N         | N     |
| Phone Messages                | RAPIDS   | Y        | Y         | Y     |
| Phone Screen                  | RAPIDS   | Y        | Y         | Y     |
| Phone WiFi Connected          | RAPIDS   | Y        | Y         | Y     |
| Phone WiFi Visible            | RAPIDS   | Y        | Y         | Y     |
| Fitbit Calories Intraday      | RAPIDS   | Y        | Y         | Y     |
| Fitbit Data Yield             | RAPIDS   | Y        | Y         | Y     |
| Fitbit Heart Rate Summary     | RAPIDS   | Y        | Y         | Y     |
| Fitbit Heart Rate Intraday    | RAPIDS   | Y        | Y         | Y     |
| Fitbit Sleep Summary          | RAPIDS   | Y        | Y         | Y     |
| Fitbit Sleep Intraday         | RAPIDS   | Y        | Y         | Y     |
| Fitbit Sleep Intraday         | PRICE    | Y        | Y         | Y     |
| Fitbit Steps Summary          | RAPIDS   | Y        | Y         | Y     |
| Fitbit Steps Intraday         | RAPIDS   | Y        | Y         | Y     |

## Accelerometer

Description

- The raw accelerometer data file, `phone_accelerometer_raw.csv`, contains data for 4 separate days
- One episode for each daily segment (night, morning, afternoon and evening)
- Two episodes locate in the same 30-min segment (`Fri 00:15:00` and `Fri 00:21:21`)
- Two episodes locate in the same daily segment (`Fri 00:15:00` and `Fri 18:12:00`)
- One episode before the time switch (`Sun 00:02:00`) and one episode after the time switch (`Sun 04:18:00`)
- Multiple episodes within one min which cause variance in magnitude (`Fri 00:10:25`, `Fri 00:10:27` and `Fri 00:10:46`)

Checklist

|time segment| single tz | multi tz|platform|
|-|-|-|-|
|30min|OK|OK|android, ios|
|morning|OK|OK|android, ios|
|daily|OK|OK|android, ios|
|threeday|OK|OK|android, ios|
|weekend|OK|OK|android, ios|
|beforeMarchEvent|OK|OK|android, ios|
|beforeNovemberEvent|OK|OK|android, ios|

## Messages (SMS)

Description

- The raw message data file, `phone_messages_raw.csv`, contains data for 4 separate days
- One episode for each daily segment (night, morning, afternoon and evening)
- Two `sent` episodes locate in the same 30-min segment (`Fri 16:08:03.000` and `Fri 16:19:35.000`)
- Two `received` episodes locate in the same 30-min segment (`Sat 06:45:05.000` and `Fri 06:45:05.000`)
- Two episodes locate in the same daily segment (`Fri 11:57:56.385` and `Sat 10:54:10.000`)
- One episode before the time switch (`Sun 00:48:01.000`) and one episode after the time switch (`Sun 06:21:01.000`)

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

## Calls

Due to the difference in the format of the raw data for iOS and Android the following is the expected results 
the `phone_calls.csv`. 

Description

- One missed episode, one outgoing episode and one incoming episode on Friday night, morning, afternoon and evening
- There is at least one episode of each type of phone calls on each day
- One incoming episode crossing two 30-mins segments
- One outgoing episode crossing two 30-mins segments
- One missed episode before, during and after the `event`
- There is one incoming episode before, during or after the `event`
- There is one outcoming episode before, during or after the `event`
- There is one missed episode before, during or after the `event`

Data format

| Device | Missed | Outgoing | Incoming |
|-|-|-|-|
|android| 3 | 2 | 1 |
|ios| 1,4 or 3,4 | 3,2,4 | 1,2,4 |

Note
When generating test data, all traces for iOS device need to be unique otherwise the episode with duplicate trace will be dropped 

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

## Screen

Due to the difference in the format of the raw screen data for iOS and Android the following is the expected results the `phone_screen.csv`. 

Description

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

Data format

| Device | unlock |
|-|-|
| Android | 3, 0|
| iOS | 3, 2|


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

## Battery

Description

- The 4-day raw data is contained in `phone_battery_raw.csv`
- One discharge episode acrossing two 30-min time segements (`Fri 05:57:30.123` to `Fri 06:04:32.456`)
- One charging episode acrossing two 30-min time segments (`Fri 11:55:58.416` to `Fri 12:08:07.876`)
- One discharge episode and one charging episode locate within the same 30-min time segement (`Fri 21:30:00` to `Fri 22:00:00`)
- One episode before the time switch (`Sun 00:24:00.000`) and one episode after the time switch (`Sun 21:58:00`)
- Two episodes locate in the same daily segment
  
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

## Bluetooth

Description 

- The 4-day raw data is contained in `phone_bluetooth_raw.csv`
- One episode for each daily segment (`night`, `morning`, `afternoon` and `evening`)
- Two episodes locate in the same 30-min segment (`Fri 23:38:45.789` and `Fri 23:59:59.465`)
- Two episodes locate in the same daily segment (`Fri 00:00:00.798` and `Fri 00:49:04.132`)
- One episode before the time switch (`Sun 00:24:00.000`) and one episode after the time switch (`Sun 17:32:00.000`)

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

## WIFI

There are two wifi features (`phone wifi connected` and `phone wifi visible`). The raw test data are seperatly stored in the `phone_wifi_connected_raw.csv` and `phone_wifi_visible_raw.csv`.

Description 

- One episode for each `epoch` (`night`, `morining`, `afternoon` and `evening`)
- Two two episodes in the same time segment (`daily` and `30-min`)
- Two episodes around the transition of `epochs` (e.g. one at the end of `night` and one at the beginning of `morning`) 
- One episode before and after the time switch on Sunday

phone wifi connected

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

Description

- The 4-day raw light data is contained in `phone_light_raw.csv`
- One episode for each daily segment (`night`, `morning`, `afternoon` and `evening`)
- Two episodes locate in the same 30-min segment (`Fri 00:07:27.000` and `Fri 00:12:00.000`)
- Two episodes locate in the same daily segment (`Fri 01:00:00` and `Fri 03:59:59.654`)
- One episode before the time switch (`Sun 00:08:00.000`) and one episode after the time switch (`Sun 05:36:00.000`)

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

## Locations

Description

- The participant's home location is (latitude=1, longitude=1).
- From Sat 10:56:00 to Sat 11:04:00, the center of the cluster is (latitude=-100, longitude=-100).
- From Sun 03:30:00 to Sun 03:47:00, the center of the cluster is (latitude=1, longitude=1). Home location is extracted from this period.
- From Sun 11:30:00 to Sun 11:38:00, the center of the cluster is (latitude=100, longitude=100).

## Application Foreground

- The 4-day raw application data is contained in `phone_applications_foreground_raw.csv`
- One episode for each daily segment (night, morning, afternoon and evening)
- Two episodes locate in the same 30-min segment (`Fri 10:12:56.385` and `Fri 10:18:48.895`)
- Two episodes locate in the same daily segment (`Fri 11:57:56.385` and `Fri 12:02:56.385`)
- One episode before the time switch (`Sun 00:07:48.001`) and one episode after the time switch (`Sun 05:10:30.001`)
- Two custom category (`Dating`) episode, one at `Fri 06:05:10.385`, another one at ` Fri 11:53:00.385`
- One episode for an application containing a special character (`Pokémon GO`) at `Mon 21:29:46.001`

Checklist:

|time segment| single tz | multi tz|platform|
|-|-|-|-|
|30min|OK|OK|android|
|morning|OK|OK|android|
|daily|OK|OK|android|
|threeday|OK|OK|android|
|weekend|OK|OK|android|
|beforeMarchEvent|OK|OK|android|
|beforeNovemberEvent|OK|OK|android|

## Activity Recognition

Description

- The 4-day raw activity data is contained in `plugin_google_activity_recognition_raw.csv` and `plugin_ios_activity_recognition_raw.csv`.
- Two episodes locate in the same 30-min segment (`Fri 04:01:54` and `Fri 04:13:52`)
- One episode for each daily segment (`night`, `morning`, `afternoon` and `evening`)
- Two episodes locate in the same daily segment (`Fri  05:03:09` and `Fri 05:50:36`)
- Two episodes with the time difference less than `5 mins` threshold (`Fri 07:14:21` and `Fri 07:18:50`)
- One episode before the time switch (`Sun 00:46:00`) and one episode after the time switch (`Sun 03:42:00`)

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

## Conversation

The 4-day raw conversation data is contained in `phone_conversation_raw.csv`. The different `inference` records are 
randomly distributed throughout the `epoch`. 

Description

- One episode for each daily segment (`night`, `morning`, `afternoon` and `evening`) on each day
- Two episodes near the transition of the daily segment, one starts at the end of the afternoon, `Fri 17:10:00` and another one starts at the beginning of the evening, `Fri 18:01:00`
- One episode across two segments, `daily` and `30-mins`, (from `Fri 05:55:00` to `Fri 06:00:41`)
- Two episodes locate in the same daily segment (`Sat 12:45:36` and `Sat 16:48:22`)
- One episode before the time switch, `Sun 00:15:06`, and one episode after the time switch, `Sun 06:01:00`

Data format

| inference | type |
| - | - |
| 0 | silence |
| 1 | noise | 
| 2 | voice |
| 3 | unknown | 

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
## Application Episodes

-   The feature requires raw application foreground data file and raw phone screen data file
-   The raw data files contains data for 4 day.
-   The raw conversation data contains records with difference in `timestamp` ranging from milliseconds to minutes.
-   An app episode starts when an app is launched and ends when another app is launched, marking the episode end of the first one,
or when the screen locks. Thus, we are taking into account the screen unlock episodes.
-   There are multiple apps usage within each screen unlock episode to verify creation of different app episodes in each 
screen unlock session. In the screen unlock episode starting from Fri 05:56:51, Fri 10:00:24, Sat 17:48:01, Sun 22:02:00, and Mon 21:05:00 we have multiple apps, both system and non-system apps, to check this.
-   The 22 minute chunk starting from Fri 10:03:56 checks app episodes for system apps only.
-   The screen unlock episode starting from Mon 21:05:00 and Sat 17:48:01 checks if the screen lock marks the end of episode for that particular app which was launched a few milliseconds to 8 mins before the screen lock.
-   Finally, since application foreground is only for Android devices, this feature is also for Android devices only. All other files are empty data files


## Data Yield

Description

- Two sensors were picked for testing, `phone_screen` and `phone_light`. `phone_screen` is event based and `phone_light` is sampling at regular frequency
- A 31-min episode (from `Fri 01:00:00` to `Fri 01:30:00`) in phone_light data, which is considered as a `validyieldedhours`


Checklist

|time segment| single tz | multi tz|platform|
|-|-|-|-|
|30min|OK|OK|android, ios|
|morning|OK|OK|android, ios|
|daily|OK|OK|android, ios|
|threeday|OK|OK|android, ios|
|weekend|OK|OK|android, ios|
|beforeMarchEvent|OK|OK|android, ios|
|beforeNovemberEvent|OK|OK|android, ios|


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


## Fitbit Heartrate intraday 

Description:

- The 4-day raw heartrate data is contained in `fitbit_heartrate_intraday_raw.csv`
- One episode for each daily segment (`night`, `morning`, `afternoon` and `evening`)
- Two episodes locate in the same 30-min segment (`Fri 00:49:00` and `Fri 00:52:00`)
- Two different types of heartrate zone episodes locate in the same 30-min segment (`Fri 05:49:00 outofrange` and `Fri 05:57:00 fatburn`)
- Two episodes locate in the same daily segment (`Fri 12:02:00` and `Fri 19:38:00`)
- One episode before the time switch, `Sun 00:08:00`, and one episode after the time switch, `Sun 07:28:00`


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


## Fitbit Heartrate Summary

Description

- The 4-day raw heartrate summary data is contained in `fitbit_heartrate_summary_raw.csv`.
- As heartrate summary is periodic, it only generates results in periodic feature, there will be no result in frequency and event. 


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

## Fitbit Step Intraday

Description

- The 4-day raw step summary data is contained in `fitbit_steps_intraday_raw.csv`
- One episode for each daily segment (`night`, `morning`, `afternoon` and `evening`) on each day
- Two episodes within the same 30-min segment (`Fri 05:58:00` and `Fri 05:59:00`)
- A one-min episode at `2020-03-07 09:00:00` that will be converted to New York time `2020-03-07 12:00:00`
- One episode before the time switch, `Sun 00:19:00`, and one episode after the time switch, `Sun 09:01:00`
- Episodes cross two 30-min segments (`Fri 11:59:00` and `Fri 12:00:00`)

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


## Fitbit Step Summary

Description

- The 4-day calculated step summary data is contained in `fitbit_steps_summary_raw.csv`.
- The 4-day calculated step volatility summary data is contained in`fitbit_steps_summary_raw.csv`.
- Step summary including max, min, median, mean and standard deviation value.
- Volatility summary including max, min, median, mean, standard deviation and annulized volatility value.
- As step summary is periodic, it only generates results in periodic feature, there will be no result in frequency and event. 


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

## Fitbit Data Yield

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