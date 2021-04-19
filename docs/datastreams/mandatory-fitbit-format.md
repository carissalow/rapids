# Mandatory Fitbit Format

This is a description of the format RAPIDS needs to process data for the following Fitbit\ sensors.

??? info "FITBIT_HEARTRATE_SUMMARY"

    | RAPIDS column   | Description   |
    |-----------------|-----------------|
    | TIMESTAMP       |  An UNIX timestamp (13 digits) when a row of data was logged (automatically created by RAPIDS) |
    | LOCAL_DATE_TIME       |  Date time string with format `yyyy-mm-dd hh:mm:ss` |
    | DEVICE_ID       |  A string that uniquely identifies a device |
    | HEARTRATE_DAILY_RESTINGHR |  Daily resting heartrate |
    | HEARTRATE_DAILY_CALORIESOUTOFRANGE |  Calories spent while heartrate was oustide a heartrate [zone](https://help.fitbit.com/articles/en_US/Help_article/1565.htm#) |
    | HEARTRATE_DAILY_CALORIESFATBURN |  Calories spent while heartrate was inside the fat burn [zone](https://help.fitbit.com/articles/en_US/Help_article/1565.htm#) |
    | HEARTRATE_DAILY_CALORIESCARDIO |  Calories spent while heartrate was inside the cardio [zone](https://help.fitbit.com/articles/en_US/Help_article/1565.htm#) |
    | HEARTRATE_DAILY_CALORIESPEAK |  Calories spent while heartrate was inside the peak [zone](https://help.fitbit.com/articles/en_US/Help_article/1565.htm#) |

??? info "FITBIT_HEARTRATE_INTRADAY"

    | RAPIDS column   | Description   |
    |-----------------|-----------------|
    | TIMESTAMP       |  An UNIX timestamp (13 digits) when a row of data was logged (automatically created by RAPIDS) |
    | LOCAL_DATE_TIME       |  Date time string with format `yyyy-mm-dd hh:mm:ss` |
    | DEVICE_ID       |  A string that uniquely identifies a device |
    | HEARTRATE |  Intraday heartrate |
    | HEARTRATE_ZONE |  Heartrate [zone](https://help.fitbit.com/articles/en_US/Help_article/1565.htm#) that HEARTRATE belongs to. It is based on the heartrate zone ranges of each device |

??? info "FITBIT_SLEEP_SUMMARY"

    | RAPIDS column   | Description   |
    |-----------------|-----------------|
    | TIMESTAMP       |  An UNIX timestamp (13 digits) when a row of data was logged (automatically created by RAPIDS) |
    | LOCAL_DATE_TIME       |  Date time string with format `yyyy-mm-dd 00:00:00`, the date is the same as the start date of a daily sleep episode if its time is after SLEEP_SUMMARY_LAST_NIGHT_END, otherwise it is the day before the start date of that sleep episode |
    | LOCAL_START_DATE_TIME       |  Date time string with format `yyyy-mm-dd hh:mm:ss` representing the start of a daily sleep episode |
    | LOCAL_END_DATE_TIME       |  Date time string with format `yyyy-mm-dd hh:mm:ss`  representing the end of a daily sleep episode|
    | DEVICE_ID       |  A string that uniquely identifies a device |
    | EFFICIENCY | Sleep efficiency computed by fitbit as time asleep / (total time in bed - time to fall asleep)|
    | MINUTES_AFTER_WAKEUP | Minutes the participant spent in bed after waking up|
    | MINUTES_ASLEEP | Minutes the participant was asleep |
    | MINUTES_AWAKE | Minutes the participant was awake |
    | MINUTES_TO_FALL_ASLEEP | Minutes the participant spent in bed before falling asleep|
    | MINUTES_IN_BED | Minutes the participant spent in bed across the sleep episode|
    | IS_MAIN_SLEEP | 0 if this episode is a nap, or 1 if it is a main sleep episode|
    | TYPE | stages or classic [sleep data](https://dev.fitbit.com/build/reference/web-api/sleep/)|

??? info "FITBIT_SLEEP_INTRADAY"

    | RAPIDS column   | Description   |
    |-----------------|-----------------|
    | TIMESTAMP       |  An UNIX timestamp (13 digits) when a row of data was logged (automatically created by RAPIDS)|
    | LOCAL_DATE_TIME       |  Date time string with format `yyyy-mm-dd hh:mm:ss`, this either is a copy of LOCAL_START_DATE_TIME or LOCAL_END_DATE_TIME depending on which column is used to assign an episode to a specific day|
    | DEVICE_ID       |  A string that uniquely identifies a device |
    | TYPE_EPISODE_ID | An id for each unique main or nap episode. Main and nap episodes have different levels, each row in this table is one of such levels, so multiple rows can have the same TYPE_EPISODE_ID|
    | DURATION | Duration of the episode level in minutes|
    | IS_MAIN_SLEEP | 0 if this episode level belongs to a nap, or 1 if it belongs to a main sleep episode|
    | TYPE | type of level: stages or classic [sleep data](https://dev.fitbit.com/build/reference/web-api/sleep/)|
    | LEVEL | For stages levels one of `wake`, `deep`, `light`, or `rem`. For classic levels one of `awake`, `restless`, and `asleep`|

??? info "FITBIT_STEPS_SUMMARY"

    | RAPIDS column   | Description   |
    |-----------------|-----------------|
    | TIMESTAMP       |  An UNIX timestamp (13 digits) when a row of data was logged (automatically created by RAPIDS) |
    | LOCAL_DATE_TIME       |  Date time string with format `yyyy-mm-dd hh:mm:ss` |
    | DEVICE_ID       |  A string that uniquely identifies a device |
    | STEPS |  Daily step count |

??? info "FITBIT_STEPS_INTRADAY"

    | RAPIDS column   | Description   |
    |-----------------|-----------------|
    | TIMESTAMP       |  An UNIX timestamp (13 digits) when a row of data was logged (automatically created by RAPIDS) |
    | LOCAL_DATE_TIME       |  Date time string with format `yyyy-mm-dd hh:mm:ss` |
    | DEVICE_ID       |  A string that uniquely identifies a device |
    | STEPS |  Intraday step count (usually every minute)|