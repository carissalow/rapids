# `fitbitparsed_mysql`
This [data stream](../../datastreams/data-streams-introduction) handles Fitbit sensor data downloaded using the [Fitbit Web API](https://dev.fitbit.com/build/reference/web-api/), **parsed**, and stored in a MySQL database. Please note that RAPIDS cannot query the API directly, you need to use other available tools or implement your own. Once you have your sensor data in a MySQL database, RAPIDS can process it.

!!! info "What is the difference between JSON and plain data streams"
    Most people will only need `fitbitjson_mysql` because they downloaded and stored their data directly from Fitbit's API. However, if for some reason you don't have access to that JSON data and instead only have the parsed data (columns and rows) you can use this data stream.

## Container
A MySQL database with a table per sensor, each containing the data for all participants.

The script to connect and download data from this container is at:
```bash
src/data/streams/fitbitjson_mysql/container.R
```

## Format

The `format.yaml` maps and transforms columns in your raw data stream to the [mandatory columns RAPIDS needs for Fitbit sensors](../mandatory-fitbit-format). This file is at:

```bash
src/data/streams/fitbitparsed_mysql/format.yaml
```

If you want RAPIDS to process Fitbit sensor data using this stream, you will need to replace any `COLUMN_MAPPINGS` inside **each sensor** section in  `format.yaml` to match your raw data column names. 

All columns are mandatory, however, all except `device_id` and `local_date_time` can be empty if you don't have that data. Just have in mind that some features will be empty if some of these columns are empty.

??? info "FITBIT_HEARTRATE_SUMMARY section"

    
    **COLUMN_MAPPINGS**

    | RAPIDS column   | Stream column   |
    |-----------------|-----------------|
    | LOCAL_DATE_TIME       | local_date_time |
    | DEVICE_ID       | device_id |
    | HEARTRATE_DAILY_RESTINGHR | heartrate_daily_restinghr |
    | HEARTRATE_DAILY_CALORIESOUTOFRANGE | heartrate_daily_caloriesoutofrange |
    | HEARTRATE_DAILY_CALORIESFATBURN | heartrate_daily_caloriesfatburn |
    | HEARTRATE_DAILY_CALORIESCARDIO | heartrate_daily_caloriescardio |
    | HEARTRATE_DAILY_CALORIESPEAK | heartrate_daily_caloriespeak |

    **MUTATION_SCRIPTS**

    TODO list our parsing script

    ??? "Example of the raw data RAPIDS expects for this data stream"

        |device_id                              |local_date_time   |heartrate_daily_restinghr |heartrate_daily_caloriesoutofrange  |heartrate_daily_caloriesfatburn  |heartrate_daily_caloriescardio  |heartrate_daily_caloriespeak   |
        |-------------------------------------- |----------------- |------- |-------------- |------------- |------------ |-------|
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-07        |72      |1200.6102      |760.3020      |15.2048      |0      |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-08        |70      |1100.1120      |660.0012      |23.7088      |0      |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-09        |69      |750.3615       |734.1516      |131.8579     |0      |
