# `fitbitjson_mysql`
This [data stream](../../datastreams/data-streams-introduction) handles Fitbit sensor data downloaded using the [Fitbit Web API](https://dev.fitbit.com/build/reference/web-api/) and stored in a MySQL database. Please note that RAPIDS cannot query the API directly, you need to use other available tools or implement your own. Once you have your sensor data in a MySQL database, RAPIDS can process it.

## Container
A MySQL database with a table per sensor, each containing the data for all participants.

The script to connect and download data from this container is at:
```bash
src/data/streams/fitbitjson_mysql/container.R
```

## Format

The `format.yaml` maps and transforms columns in your raw data stream to the [mandatory columns RAPIDS needs for Fitbit sensors](../mandatory-fitbit-format). This file is at:

```bash
src/data/streams/fitbitjson_mysql/format.yaml
```

If you want RAPIDS to process Fitbit sensor data using this stream, you will need to replace the following `COLUMN_MAPPINGS` inside **each sensor** section in `format.yaml` to match your raw data column names:

| Column   | Description   |
|-----------------|-----------------|
| device_id       | A string that uniquely identifies a device |
| fitbit_data       | A string column that contains the JSON objects downloaded from Fitbit's API |



??? info "FITBIT_HEARTRATE_SUMMARY section"

    **COLUMN_MAPPINGS**

    | RAPIDS column   | Stream column   |
    |-----------------|-----------------|
    | LOCAL_DATE_TIME       | FLAG_TO_MUTATE |
    | DEVICE_ID       | device_id |
    | HEARTRATE_DAILY_RESTINGHR | FLAG_TO_MUTATE |
    | HEARTRATE_DAILY_CALORIESOUTOFRANGE | FLAG_TO_MUTATE |
    | HEARTRATE_DAILY_CALORIESFATBURN | FLAG_TO_MUTATE |
    | HEARTRATE_DAILY_CALORIESCARDIO | FLAG_TO_MUTATE |
    | HEARTRATE_DAILY_CALORIESPEAK | FLAG_TO_MUTATE |
    | FLAG_AS_EXTRA: | fitbit_data |


    **MUTATION_SCRIPTS**

    TODO list our parsing script

    ??? "Example of the raw data RAPIDS expects for this data stream"

        |device_id                                |fitbit_data                                               |
        |---------------------------------------- |--------------------------------------------------------- |
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"activities-heart":[{"dateTime":"2020-10-07","value":{"customHeartRateZones":[],"heartRateZones":[{"caloriesOut":1200.6102,"max":88,"min":31,"minutes":1058,"name":"Out of Range"},{"caloriesOut":760.3020,"max":120,"min":86,"minutes":366,"name":"Fat Burn"},{"caloriesOut":15.2048,"max":146,"min":120,"minutes":2,"name":"Cardio"},{"caloriesOut":0,"max":221,"min":148,"minutes":0,"name":"Peak"}],"restingHeartRate":72}}],"activities-heart-intraday":{"dataset":[{"time":"00:00:00","value":68},{"time":"00:01:00","value":67},{"time":"00:02:00","value":67},...],"datasetInterval":1,"datasetType":"minute"}}
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"activities-heart":[{"dateTime":"2020-10-08","value":{"customHeartRateZones":[],"heartRateZones":[{"caloriesOut":1100.1120,"max":89,"min":30,"minutes":921,"name":"Out of Range"},{"caloriesOut":660.0012,"max":118,"min":82,"minutes":361,"name":"Fat Burn"},{"caloriesOut":23.7088,"max":142,"min":108,"minutes":3,"name":"Cardio"},{"caloriesOut":0,"max":221,"min":148,"minutes":0,"name":"Peak"}],"restingHeartRate":70}}],"activities-heart-intraday":{"dataset":[{"time":"00:00:00","value":77},{"time":"00:01:00","value":75},{"time":"00:02:00","value":73},...],"datasetInterval":1,"datasetType":"minute"}}
        |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"activities-heart":[{"dateTime":"2020-10-09","value":{"customHeartRateZones":[],"heartRateZones":[{"caloriesOut":750.3615,"max":77,"min":30,"minutes":851,"name":"Out of Range"},{"caloriesOut":734.1516,"max":107,"min":77,"minutes":550,"name":"Fat Burn"},{"caloriesOut":131.8579,"max":130,"min":107,"minutes":29,"name":"Cardio"},{"caloriesOut":0,"max":220,"min":130,"minutes":0,"name":"Peak"}],"restingHeartRate":69}}],"activities-heart-intraday":{"dataset":[{"time":"00:00:00","value":90},{"time":"00:01:00","value":89},{"time":"00:02:00","value":88},...],"datasetInterval":1,"datasetType":"minute"}}