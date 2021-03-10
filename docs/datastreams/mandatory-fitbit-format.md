# Mandatory Fitbit Format

This is a description of the format RAPIDS needs to process data for the following Fitbit\ sensors.

??? info "FITBIT_HEARTRATE_SUMMARY"

    | RAPIDS column   | Description   |
    |-----------------|-----------------|
    | TIMESTAMP       |  An UNIX timestamp (13 digits) when a row of data was logged |
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
    | TIMESTAMP       |  An UNIX timestamp (13 digits) when a row of data was logged |
    | LOCAL_DATE_TIME       |  Date time string with format `yyyy-mm-dd hh:mm:ss` |
    | DEVICE_ID       |  A string that uniquely identifies a device |
    | HEARTRATE |  Intraday heartrate |
    | HEARTRATE_ZONE |  Heartrate [zone](https://help.fitbit.com/articles/en_US/Help_article/1565.htm#) that HEARTRATE belongs to. It is based on the heartrate zone ranges of each device |

??? info "FITBIT_STEPS_SUMMARY"

    | RAPIDS column   | Description   |
    |-----------------|-----------------|
    | TIMESTAMP       |  An UNIX timestamp (13 digits) when a row of data was logged |
    | LOCAL_DATE_TIME       |  Date time string with format `yyyy-mm-dd hh:mm:ss` |
    | DEVICE_ID       |  A string that uniquely identifies a device |
    | STEPS |  Daily step count |

??? info "FITBIT_STEPS_INTRADAY"

    | RAPIDS column   | Description   |
    |-----------------|-----------------|
    | TIMESTAMP       |  An UNIX timestamp (13 digits) when a row of data was logged |
    | LOCAL_DATE_TIME       |  Date time string with format `yyyy-mm-dd hh:mm:ss` |
    | DEVICE_ID       |  A string that uniquely identifies a device |
    | STEPS |  Intraday step count (usually every minute)|