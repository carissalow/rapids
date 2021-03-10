# Mandatory Fitbit Format

This is a description of the format RAPIDS needs to process data for the following Fitbit\ sensors.

??? info "FITBIT_HEARTRATE_SUMMARY"

    | RAPIDS column   | Description   |
    |-----------------|-----------------|
    | LOCAL_DATE_TIME       |  Date time string with format `yyyy-mm-dd hh:mm:ss` |
    | DEVICE_ID       |  A string that uniquely identifies a device |
    | HEARTRATE_DAILY_RESTINGHR |  TODO |
    | HEARTRATE_DAILY_CALORIESOUTOFRANGE |  TODO |
    | HEARTRATE_DAILY_CALORIESFATBURN |  TODO |
    | HEARTRATE_DAILY_CALORIESCARDIO |  TODO |
    | HEARTRATE_DAILY_CALORIESPEAK |  TODO |

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