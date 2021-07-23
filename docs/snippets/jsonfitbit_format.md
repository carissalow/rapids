
The `format.yaml` maps and transforms columns in your raw data stream to the [mandatory columns RAPIDS needs for Fitbit sensors](../mandatory-fitbit-format). This file is at:

```bash
src/data/streams/fitbitjson_csv/format.yaml
```

If you want RAPIDS to process Fitbit sensor data using this stream, you will need to map `DEVICE_ID` and `JSON_FITBIT_COLUMN` to your own raw data columns inside **each sensor** section in `format.yaml`.

??? info "FITBIT_HEARTRATE_SUMMARY"

    **RAPIDS_COLUMN_MAPPINGS**

    | RAPIDS column   | Stream column   |
    |-----------------|-----------------|
    | LOCAL_DATE_TIME       | FLAG_TO_MUTATE |
    | DEVICE_ID       | device_id |
    | HEARTRATE_DAILY_RESTINGHR | FLAG_TO_MUTATE |
    | HEARTRATE_DAILY_CALORIESOUTOFRANGE | FLAG_TO_MUTATE |
    | HEARTRATE_DAILY_CALORIESFATBURN | FLAG_TO_MUTATE |
    | HEARTRATE_DAILY_CALORIESCARDIO | FLAG_TO_MUTATE |
    | HEARTRATE_DAILY_CALORIESPEAK | FLAG_TO_MUTATE |


    **MUTATION**

    - **COLUMN_MAPPINGS**

        | Script column   | Stream column   |
        |-----------------|-----------------|
        | JSON_FITBIT_COLUMN      | fitbit_data      |
    
    - **SCRIPTS**
    
        ```bash
        - src/data/streams/mutations/fitbit/parse_heartrate_summary_json.py
        - src/data/streams/mutations/fitbit/add_zero_timestamp.py
        ```

        !!! note
            All columns except `DEVICE_ID` are parsed from `JSON_FITBIT_COLUMN`. `JSON_FITBIT_COLUMN` is a string column containing the JSON objects returned by Fitbit's API. See an example of the raw data RAPIDS expects for this data stream:


            ??? example "Example of the raw data RAPIDS expects for this data stream"

                |device_id                                |fitbit_data                                               |
                |---------------------------------------- |--------------------------------------------------------- |
                |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"activities-heart":[{"dateTime":"2020-10-07","value":{"customHeartRateZones":[],"heartRateZones":[{"caloriesOut":1200.6102,"max":88,"min":31,"minutes":1058,"name":"Out of Range"},{"caloriesOut":760.3020,"max":120,"min":86,"minutes":366,"name":"Fat Burn"},{"caloriesOut":15.2048,"max":146,"min":120,"minutes":2,"name":"Cardio"},{"caloriesOut":0,"max":221,"min":148,"minutes":0,"name":"Peak"}],"restingHeartRate":72}}],"activities-heart-intraday":{"dataset":[{"time":"00:00:00","value":68},{"time":"00:01:00","value":67},{"time":"00:02:00","value":67},...],"datasetInterval":1,"datasetType":"minute"}}
                |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"activities-heart":[{"dateTime":"2020-10-08","value":{"customHeartRateZones":[],"heartRateZones":[{"caloriesOut":1100.1120,"max":89,"min":30,"minutes":921,"name":"Out of Range"},{"caloriesOut":660.0012,"max":118,"min":82,"minutes":361,"name":"Fat Burn"},{"caloriesOut":23.7088,"max":142,"min":108,"minutes":3,"name":"Cardio"},{"caloriesOut":0,"max":221,"min":148,"minutes":0,"name":"Peak"}],"restingHeartRate":70}}],"activities-heart-intraday":{"dataset":[{"time":"00:00:00","value":77},{"time":"00:01:00","value":75},{"time":"00:02:00","value":73},...],"datasetInterval":1,"datasetType":"minute"}}
                |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"activities-heart":[{"dateTime":"2020-10-09","value":{"customHeartRateZones":[],"heartRateZones":[{"caloriesOut":750.3615,"max":77,"min":30,"minutes":851,"name":"Out of Range"},{"caloriesOut":734.1516,"max":107,"min":77,"minutes":550,"name":"Fat Burn"},{"caloriesOut":131.8579,"max":130,"min":107,"minutes":29,"name":"Cardio"},{"caloriesOut":0,"max":220,"min":130,"minutes":0,"name":"Peak"}],"restingHeartRate":69}}],"activities-heart-intraday":{"dataset":[{"time":"00:00:00","value":90},{"time":"00:01:00","value":89},{"time":"00:02:00","value":88},...],"datasetInterval":1,"datasetType":"minute"}}

??? info "FITBIT_HEARTRATE_INTRADAY"

    **RAPIDS_COLUMN_MAPPINGS**

    | RAPIDS column   | Stream column   |
    |-----------------|-----------------|
    | LOCAL_DATE_TIME       | FLAG_TO_MUTATE |
    | DEVICE_ID       | device_id |
    | HEARTRATE | FLAG_TO_MUTATE |
    | HEARTRATE_ZONE | FLAG_TO_MUTATE |


    **MUTATION**

    - **COLUMN_MAPPINGS**

        | Script column   | Stream column   |
        |-----------------|-----------------|
        | JSON_FITBIT_COLUMN      | fitbit_data      |
    
    - **SCRIPTS**
    
        ```bash
        - src/data/streams/mutations/fitbit/parse_heartrate_intraday_json.py
        - src/data/streams/mutations/fitbit/add_zero_timestamp.py
        ```

        !!! note
            All columns except `DEVICE_ID` are parsed from `JSON_FITBIT_COLUMN`. `JSON_FITBIT_COLUMN` is a string column containing the JSON objects returned by Fitbit's API. See an example of the raw data RAPIDS expects for this data stream:


            ??? example "Example of the raw data RAPIDS expects for this data stream"

                |device_id                                |fitbit_data                                               |
                |---------------------------------------- |--------------------------------------------------------- |
                |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"activities-heart":[{"dateTime":"2020-10-07","value":{"customHeartRateZones":[],"heartRateZones":[{"caloriesOut":1200.6102,"max":88,"min":31,"minutes":1058,"name":"Out of Range"},{"caloriesOut":760.3020,"max":120,"min":86,"minutes":366,"name":"Fat Burn"},{"caloriesOut":15.2048,"max":146,"min":120,"minutes":2,"name":"Cardio"},{"caloriesOut":0,"max":221,"min":148,"minutes":0,"name":"Peak"}],"restingHeartRate":72}}],"activities-heart-intraday":{"dataset":[{"time":"00:00:00","value":68},{"time":"00:01:00","value":67},{"time":"00:02:00","value":67},...],"datasetInterval":1,"datasetType":"minute"}}
                |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"activities-heart":[{"dateTime":"2020-10-08","value":{"customHeartRateZones":[],"heartRateZones":[{"caloriesOut":1100.1120,"max":89,"min":30,"minutes":921,"name":"Out of Range"},{"caloriesOut":660.0012,"max":118,"min":82,"minutes":361,"name":"Fat Burn"},{"caloriesOut":23.7088,"max":142,"min":108,"minutes":3,"name":"Cardio"},{"caloriesOut":0,"max":221,"min":148,"minutes":0,"name":"Peak"}],"restingHeartRate":70}}],"activities-heart-intraday":{"dataset":[{"time":"00:00:00","value":77},{"time":"00:01:00","value":75},{"time":"00:02:00","value":73},...],"datasetInterval":1,"datasetType":"minute"}}
                |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"activities-heart":[{"dateTime":"2020-10-09","value":{"customHeartRateZones":[],"heartRateZones":[{"caloriesOut":750.3615,"max":77,"min":30,"minutes":851,"name":"Out of Range"},{"caloriesOut":734.1516,"max":107,"min":77,"minutes":550,"name":"Fat Burn"},{"caloriesOut":131.8579,"max":130,"min":107,"minutes":29,"name":"Cardio"},{"caloriesOut":0,"max":220,"min":130,"minutes":0,"name":"Peak"}],"restingHeartRate":69}}],"activities-heart-intraday":{"dataset":[{"time":"00:00:00","value":90},{"time":"00:01:00","value":89},{"time":"00:02:00","value":88},...],"datasetInterval":1,"datasetType":"minute"}}

??? info "FITBIT_SLEEP_SUMMARY"

    **RAPIDS_COLUMN_MAPPINGS**

    | RAPIDS column   | Stream column   |
    |-----------------|-----------------|
    | TIMESTAMP | FLAG_TO_MUTATE |
    | LOCAL_DATE_TIME       | FLAG_TO_MUTATE |
    | LOCAL_START_DATE_TIME | FLAG_TO_MUTATE |
    | LOCAL_END_DATE_TIME | FLAG_TO_MUTATE |
    | DEVICE_ID | device_id |
    | EFFICIENCY | FLAG_TO_MUTATE |
    | MINUTES_AFTER_WAKEUP | FLAG_TO_MUTATE |
    | MINUTES_ASLEEP | FLAG_TO_MUTATE |
    | MINUTES_AWAKE | FLAG_TO_MUTATE |
    | MINUTES_TO_FALL_ASLEEP | FLAG_TO_MUTATE |
    | MINUTES_IN_BED | FLAG_TO_MUTATE |
    | IS_MAIN_SLEEP | FLAG_TO_MUTATE |
    | TYPE | FLAG_TO_MUTATE |

    **MUTATION**

    - **COLUMN_MAPPINGS**

        | Script column   | Stream column   |
        |-----------------|-----------------|
        | JSON_FITBIT_COLUMN      | fitbit_data      |
    
    - **SCRIPTS**
    
        ```bash
        - src/data/streams/mutations/fitbit/parse_sleep_summary_json.py
        - src/data/streams/mutations/fitbit/add_local_date_time.py
        - src/data/streams/mutations/fitbit/add_zero_timestamp.py
        ```

        !!! note

            Fitbit API has two versions for sleep data, v1 and v1.2. We support both but ignore v1's `count_awake`, `duration_awake`, and `count_awakenings`, `count_restless`, `duration_restless` columns.
            
            All columns except `DEVICE_ID` are parsed from `JSON_FITBIT_COLUMN`. `JSON_FITBIT_COLUMN` is a string column containing the JSON objects returned by Fitbit's API. See an example of the raw data RAPIDS expects for this data stream:

            ??? example "Example of the expected raw data"

                |device_id                                |fitbit_data                                               |
                |---------------------------------------- |--------------------------------------------------------- |
                |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"sleep":[{"dateOfSleep":"2020-10-10","duration":3600000,"efficiency":92,"endTime":"2020-10-10T16:37:00.000","infoCode":2,"isMainSleep":false,"levels":{"data":[{"dateTime":"2020-10-10T15:36:30.000","level":"restless","seconds":60},{"dateTime":"2020-10-10T15:37:30.000","level":"asleep","seconds":660},{"dateTime":"2020-10-10T15:48:30.000","level":"restless","seconds":60},...], "summary":{"asleep":{"count":0,"minutes":56},"awake":{"count":0,"minutes":0},"restless":{"count":3,"minutes":4}}},"logId":26315914306,"minutesAfterWakeup":0,"minutesAsleep":55,"minutesAwake":5,"minutesToFallAsleep":0,"startTime":"2020-10-10T15:36:30.000","timeInBed":60,"type":"classic"},{"dateOfSleep":"2020-10-10","duration":22980000,"efficiency":88,"endTime":"2020-10-10T08:10:00.000","infoCode":0,"isMainSleep":true,"levels":{"data":[{"dateTime":"2020-10-10T01:46:30.000","level":"light","seconds":420},{"dateTime":"2020-10-10T01:53:30.000","level":"deep","seconds":1230},{"dateTime":"2020-10-10T02:14:00.000","level":"light","seconds":360},...], "summary":{"deep":{"count":3,"minutes":92,"thirtyDayAvgMinutes":0},"light":{"count":29,"minutes":193,"thirtyDayAvgMinutes":0},"rem":{"count":4,"minutes":33,"thirtyDayAvgMinutes":0},"wake":{"count":28,"minutes":65,"thirtyDayAvgMinutes":0}}},"logId":26311786557,"minutesAfterWakeup":0,"minutesAsleep":318,"minutesAwake":65,"minutesToFallAsleep":0,"startTime":"2020-10-10T01:46:30.000","timeInBed":383,"type":"stages"}],"summary":{"stages":{"deep":92,"light":193,"rem":33,"wake":65},"totalMinutesAsleep":373,"totalSleepRecords":2,"totalTimeInBed":443}}
                |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"sleep":[{"dateOfSleep":"2020-10-11","duration":41640000,"efficiency":89,"endTime":"2020-10-11T11:47:00.000","infoCode":0,"isMainSleep":true,"levels":{"data":[{"dateTime":"2020-10-11T00:12:30.000","level":"wake","seconds":450},{"dateTime":"2020-10-11T00:20:00.000","level":"light","seconds":870},{"dateTime":"2020-10-11T00:34:30.000","level":"wake","seconds":780},...], "summary":{"deep":{"count":4,"minutes":52,"thirtyDayAvgMinutes":62},"light":{"count":32,"minutes":442,"thirtyDayAvgMinutes":364},"rem":{"count":6,"minutes":68,"thirtyDayAvgMinutes":58},"wake":{"count":29,"minutes":132,"thirtyDayAvgMinutes":94}}},"logId":26589710670,"minutesAfterWakeup":1,"minutesAsleep":562,"minutesAwake":132,"minutesToFallAsleep":0,"startTime":"2020-10-11T00:12:30.000","timeInBed":694,"type":"stages"}],"summary":{"stages":{"deep":52,"light":442,"rem":68,"wake":132},"totalMinutesAsleep":562,"totalSleepRecords":1,"totalTimeInBed":694}}
                |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"sleep":[{"dateOfSleep":"2020-10-12","duration":28980000,"efficiency":93,"endTime":"2020-10-12T09:34:30.000","infoCode":0,"isMainSleep":true,"levels":{"data":[{"dateTime":"2020-10-12T01:31:00.000","level":"wake","seconds":600},{"dateTime":"2020-10-12T01:41:00.000","level":"light","seconds":60},{"dateTime":"2020-10-12T01:42:00.000","level":"deep","seconds":2340},...], "summary":{"deep":{"count":4,"minutes":63,"thirtyDayAvgMinutes":59},"light":{"count":27,"minutes":257,"thirtyDayAvgMinutes":364},"rem":{"count":5,"minutes":94,"thirtyDayAvgMinutes":58},"wake":{"count":24,"minutes":69,"thirtyDayAvgMinutes":95}}},"logId":26589710673,"minutesAfterWakeup":0,"minutesAsleep":415,"minutesAwake":68,"minutesToFallAsleep":0,"startTime":"2020-10-12T01:31:00.000","timeInBed":483,"type":"stages"}],"summary":{"stages":{"deep":63,"light":257,"rem":94,"wake":69},"totalMinutesAsleep":415,"totalSleepRecords":1,"totalTimeInBed":483}}

??? info "FITBIT_SLEEP_INTRADAY"

    **RAPIDS_COLUMN_MAPPINGS**

    | RAPIDS column   | Stream column   |
    |-----------------|-----------------|
    | TIMESTAMP | FLAG_TO_MUTATE |
    | LOCAL_DATE_TIME       | FLAG_TO_MUTATE |
    | DEVICE_ID | device_id |
    | TYPE_EPISODE_ID | FLAG_TO_MUTATE |
    | DURATION | FLAG_TO_MUTATE |
    | IS_MAIN_SLEEP | FLAG_TO_MUTATE |
    | TYPE | FLAG_TO_MUTATE |
    | LEVEL | FLAG_TO_MUTATE |

    **MUTATION**

    - **COLUMN_MAPPINGS**

        | Script column   | Stream column   |
        |-----------------|-----------------|
        | JSON_FITBIT_COLUMN      | fitbit_data      |
    
    - **SCRIPTS**
    
        ```bash
        - src/data/streams/mutations/fitbit/parse_sleep_intraday_json.py
        - src/data/streams/mutations/fitbit/add_zero_timestamp.py
        ```

        !!! note

            Fitbit API has two versions for sleep data, v1 and v1.2, we support both.
            
            All columns except `DEVICE_ID` are parsed from `JSON_FITBIT_COLUMN`. `JSON_FITBIT_COLUMN` is a string column containing the JSON objects returned by Fitbit's API. See an example of the raw data RAPIDS expects for this data stream:

            ??? example "Example of the expected raw data"

                |device_id                                |fitbit_data                                               |
                |---------------------------------------- |--------------------------------------------------------- |
                |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"sleep":[{"dateOfSleep":"2020-10-10","duration":3600000,"efficiency":92,"endTime":"2020-10-10T16:37:00.000","infoCode":2,"isMainSleep":false,"levels":{"data":[{"dateTime":"2020-10-10T15:36:30.000","level":"restless","seconds":60},{"dateTime":"2020-10-10T15:37:30.000","level":"asleep","seconds":660},{"dateTime":"2020-10-10T15:48:30.000","level":"restless","seconds":60},...], "summary":{"asleep":{"count":0,"minutes":56},"awake":{"count":0,"minutes":0},"restless":{"count":3,"minutes":4}}},"logId":26315914306,"minutesAfterWakeup":0,"minutesAsleep":55,"minutesAwake":5,"minutesToFallAsleep":0,"startTime":"2020-10-10T15:36:30.000","timeInBed":60,"type":"classic"},{"dateOfSleep":"2020-10-10","duration":22980000,"efficiency":88,"endTime":"2020-10-10T08:10:00.000","infoCode":0,"isMainSleep":true,"levels":{"data":[{"dateTime":"2020-10-10T01:46:30.000","level":"light","seconds":420},{"dateTime":"2020-10-10T01:53:30.000","level":"deep","seconds":1230},{"dateTime":"2020-10-10T02:14:00.000","level":"light","seconds":360},...], "summary":{"deep":{"count":3,"minutes":92,"thirtyDayAvgMinutes":0},"light":{"count":29,"minutes":193,"thirtyDayAvgMinutes":0},"rem":{"count":4,"minutes":33,"thirtyDayAvgMinutes":0},"wake":{"count":28,"minutes":65,"thirtyDayAvgMinutes":0}}},"logId":26311786557,"minutesAfterWakeup":0,"minutesAsleep":318,"minutesAwake":65,"minutesToFallAsleep":0,"startTime":"2020-10-10T01:46:30.000","timeInBed":383,"type":"stages"}],"summary":{"stages":{"deep":92,"light":193,"rem":33,"wake":65},"totalMinutesAsleep":373,"totalSleepRecords":2,"totalTimeInBed":443}}
                |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"sleep":[{"dateOfSleep":"2020-10-11","duration":41640000,"efficiency":89,"endTime":"2020-10-11T11:47:00.000","infoCode":0,"isMainSleep":true,"levels":{"data":[{"dateTime":"2020-10-11T00:12:30.000","level":"wake","seconds":450},{"dateTime":"2020-10-11T00:20:00.000","level":"light","seconds":870},{"dateTime":"2020-10-11T00:34:30.000","level":"wake","seconds":780},...], "summary":{"deep":{"count":4,"minutes":52,"thirtyDayAvgMinutes":62},"light":{"count":32,"minutes":442,"thirtyDayAvgMinutes":364},"rem":{"count":6,"minutes":68,"thirtyDayAvgMinutes":58},"wake":{"count":29,"minutes":132,"thirtyDayAvgMinutes":94}}},"logId":26589710670,"minutesAfterWakeup":1,"minutesAsleep":562,"minutesAwake":132,"minutesToFallAsleep":0,"startTime":"2020-10-11T00:12:30.000","timeInBed":694,"type":"stages"}],"summary":{"stages":{"deep":52,"light":442,"rem":68,"wake":132},"totalMinutesAsleep":562,"totalSleepRecords":1,"totalTimeInBed":694}}
                |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |{"sleep":[{"dateOfSleep":"2020-10-12","duration":28980000,"efficiency":93,"endTime":"2020-10-12T09:34:30.000","infoCode":0,"isMainSleep":true,"levels":{"data":[{"dateTime":"2020-10-12T01:31:00.000","level":"wake","seconds":600},{"dateTime":"2020-10-12T01:41:00.000","level":"light","seconds":60},{"dateTime":"2020-10-12T01:42:00.000","level":"deep","seconds":2340},...], "summary":{"deep":{"count":4,"minutes":63,"thirtyDayAvgMinutes":59},"light":{"count":27,"minutes":257,"thirtyDayAvgMinutes":364},"rem":{"count":5,"minutes":94,"thirtyDayAvgMinutes":58},"wake":{"count":24,"minutes":69,"thirtyDayAvgMinutes":95}}},"logId":26589710673,"minutesAfterWakeup":0,"minutesAsleep":415,"minutesAwake":68,"minutesToFallAsleep":0,"startTime":"2020-10-12T01:31:00.000","timeInBed":483,"type":"stages"}],"summary":{"stages":{"deep":63,"light":257,"rem":94,"wake":69},"totalMinutesAsleep":415,"totalSleepRecords":1,"totalTimeInBed":483}}

??? info "FITBIT_STEPS_SUMMARY"

    **RAPIDS_COLUMN_MAPPINGS**

    | RAPIDS column   | Stream column   |
    |-----------------|-----------------|
    | TIMESTAMP       | FLAG_TO_MUTATE       |
    | DEVICE_ID       | device_id       |
    | LOCAL_DATE_TIME       | FLAG_TO_MUTATE       |
    | STEPS | FLAG_TO_MUTATE |

    **MUTATION**

    - **COLUMN_MAPPINGS**

        | Script column   | Stream column   |
        |-----------------|-----------------|
        | JSON_FITBIT_COLUMN      | fitbit_data      |
    
    - **SCRIPTS**
    
        ```bash
        - src/data/streams/mutations/fitbit/parse_steps_summary_json.py
        - src/data/streams/mutations/fitbit/add_zero_timestamp.py
        ```

        !!! note
            `TIMESTAMP`, `LOCAL_DATE_TIME`, and `STEPS` are parsed from `JSON_FITBIT_COLUMN`. `JSON_FITBIT_COLUMN` is a string column containing the JSON objects returned by Fitbit's API. See an example of the raw data RAPIDS expects for this data stream:

            ??? example "Example of the expected raw data"

                |device_id                                |fitbit_data                                               |
                |---------------------------------------- |--------------------------------------------------------- |
                |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |"activities-steps":[{"dateTime":"2020-10-07","value":"1775"}],"activities-steps-intraday":{"dataset":[{"time":"00:00:00","value":5},{"time":"00:01:00","value":3},{"time":"00:02:00","value":0},...],"datasetInterval":1,"datasetType":"minute"}}
                |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |"activities-steps":[{"dateTime":"2020-10-08","value":"3201"}],"activities-steps-intraday":{"dataset":[{"time":"00:00:00","value":14},{"time":"00:01:00","value":11},{"time":"00:02:00","value":10},...],"datasetInterval":1,"datasetType":"minute"}}
                |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |"activities-steps":[{"dateTime":"2020-10-09","value":"998"}],"activities-steps-intraday":{"dataset":[{"time":"00:00:00","value":0},{"time":"00:01:00","value":0},{"time":"00:02:00","value":0},...],"datasetInterval":1,"datasetType":"minute"}}
        
??? info "FITBIT_STEPS_INTRADAY"

    **RAPIDS_COLUMN_MAPPINGS**

    | RAPIDS column   | Stream column   |
    |-----------------|-----------------|
    | TIMESTAMP       | FLAG_TO_MUTATE       |
    | DEVICE_ID       | device_id       |
    | LOCAL_DATE_TIME       | FLAG_TO_MUTATE       |
    | STEPS | FLAG_TO_MUTATE |

    **MUTATION**

    - **COLUMN_MAPPINGS**

        | Script column   | Stream column   |
        |-----------------|-----------------|
        | JSON_FITBIT_COLUMN      | fitbit_data      |
    
    - **SCRIPTS**
    
        ```bash
        - src/data/streams/mutations/fitbit/parse_steps_intraday_json.py
        - src/data/streams/mutations/fitbit/add_zero_timestamp.py
        ```

        !!! note
            `TIMESTAMP`, `LOCAL_DATE_TIME`, and `STEPS` are parsed from `JSON_FITBIT_COLUMN`. `JSON_FITBIT_COLUMN` is a string column containing the JSON objects returned by [Fitbit's API](https://dev.fitbit.com/build/reference/web-api/activity/#get-activity-intraday-time-series). See an example of the raw data RAPIDS expects for this data stream:

            ??? example "Example of the expected raw data"

                |device_id                                |fitbit_data                                               |
                |---------------------------------------- |--------------------------------------------------------- |
                |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |"activities-steps":[{"dateTime":"2020-10-07","value":"1775"}],"activities-steps-intraday":{"dataset":[{"time":"00:00:00","value":5},{"time":"00:01:00","value":3},{"time":"00:02:00","value":0},...],"datasetInterval":1,"datasetType":"minute"}}
                |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |"activities-steps":[{"dateTime":"2020-10-08","value":"3201"}],"activities-steps-intraday":{"dataset":[{"time":"00:00:00","value":14},{"time":"00:01:00","value":11},{"time":"00:02:00","value":10},...],"datasetInterval":1,"datasetType":"minute"}}
                |a748ee1a-1d0b-4ae9-9074-279a2b6ba524     |"activities-steps":[{"dateTime":"2020-10-09","value":"998"}],"activities-steps-intraday":{"dataset":[{"time":"00:00:00","value":0},{"time":"00:01:00","value":0},{"time":"00:02:00","value":0},...],"datasetInterval":1,"datasetType":"minute"}}
    
        