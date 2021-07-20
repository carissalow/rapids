
The `format.yaml` maps and transforms columns in your raw data stream to the [mandatory columns RAPIDS needs for Fitbit sensors](../mandatory-fitbit-format). This file is at:

```bash
src/data/streams/fitbitparsed_mysql/format.yaml
```

If you want to use this stream with your data, modify every sensor in `format.yaml` to map all columns except `TIMESTAMP` in `[RAPIDS_COLUMN_MAPPINGS]` to your raw data column names.

All columns are mandatory; however, all except `device_id` and `local_date_time` can be empty if you don't have that data. Just have in mind that some features will be empty if some of these columns are empty.


??? info "FITBIT_HEARTRATE_SUMMARY"

    **RAPIDS_COLUMN_MAPPINGS**

    | RAPIDS column   | Stream column   |
    |-----------------|-----------------|
    | TIMESTAMP| FLAG_TO_MUTATE |
    | LOCAL_DATE_TIME | local_date_time |
    | DEVICE_ID | device_id |
    | HEARTRATE_DAILY_RESTINGHR | heartrate_daily_restinghr |
    | HEARTRATE_DAILY_CALORIESOUTOFRANGE | heartrate_daily_caloriesoutofrange |
    | HEARTRATE_DAILY_CALORIESFATBURN | heartrate_daily_caloriesfatburn |
    | HEARTRATE_DAILY_CALORIESCARDIO | heartrate_daily_caloriescardio |
    | HEARTRATE_DAILY_CALORIESPEAK | heartrate_daily_caloriespeak |


    **MUTATION**

    - **COLUMN_MAPPINGS** (None)

    - **SCRIPTS** 

        ```bash
        src/data/streams/mutations/fitbit/add_zero_timestamp.py
        ```

    !!! note
        `add_zero_timestamp` adds an all-zero column called `timestamp` that will be filled in later in the pipeline by `readable_time.R` converting LOCAL_DATE_TIME to a unix timestamp taking into account single or multiple time zones.

        ??? example "Example of the raw data RAPIDS expects for this data stream"

            |device_id                              |local_date_time   |heartrate_daily_restinghr |heartrate_daily_caloriesoutofrange  |heartrate_daily_caloriesfatburn  |heartrate_daily_caloriescardio  |heartrate_daily_caloriespeak   |
            |-------------------------------------- |----------------- |------- |-------------- |------------- |------------ |-------|
            |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-07 00:00:00  |72      |1200.6102      |760.3020      |15.2048      |0      |
            |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-08 00:00:00  |70      |1100.1120      |660.0012      |23.7088      |0      |
            |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-09 00:00:00  |69      |750.3615       |734.1516      |131.8579     |0      |

??? info "FITBIT_HEARTRATE_INTRADAY"

    **RAPIDS_COLUMN_MAPPINGS**

    | RAPIDS column   | Stream column   |
    |-----------------|-----------------|
    | TIMESTAMP| FLAG_TO_MUTATE |
    | LOCAL_DATE_TIME | local_date_time |
    | DEVICE_ID | device_id |
    | HEARTRATE | heartrate |
    | HEARTRATE_ZONE | heartrate_zone |


    **MUTATION**

    - **COLUMN_MAPPINGS** (None)
    
    - **SCRIPTS** 
    
        ```bash
        src/data/streams/mutations/fitbit/add_zero_timestamp.py
        ```

    !!! note
        `add_zero_timestamp` adds an all-zero column called `timestamp` that will be filled in later in the pipeline by `readable_time.R` converting LOCAL_DATE_TIME to a unix timestamp taking into account single or multiple time zones.

        ??? example "Example of the raw data RAPIDS expects for this data stream"

            |device_id                              |local_date_time        |heartrate |heartrate_zone  |
            |-------------------------------------- |---------------------- |--------- |--------------- |
            |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-07 00:00:00    |68        |outofrange      |
            |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-07 00:01:00    |67        |outofrange      |
            |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-07 00:02:00    |67        |outofrange      |

??? info "FITBIT_SLEEP_SUMMARY"

    **RAPIDS_COLUMN_MAPPINGS**

    | RAPIDS column   | Stream column   |
    |-----------------|-----------------|
    | TIMESTAMP| FLAG_TO_MUTATE |
    | LOCAL_DATE_TIME| FLAG_TO_MUTATE |
    | LOCAL_START_DATE_TIME| local_start_date_time |
    | LOCAL_END_DATE_TIME| local_end_date_time |
    | DEVICE_ID| device_id |
    | EFFICIENCY| efficiency |
    | MINUTES_AFTER_WAKEUP| minutes_after_wakeup |
    | MINUTES_ASLEEP| minutes_asleep |
    | MINUTES_AWAKE| minutes_awake |
    | MINUTES_TO_FALL_ASLEEP| minutes_to_fall_asleep |
    | MINUTES_IN_BED| minutes_in_bed |
    | IS_MAIN_SLEEP| is_main_sleep |
    | TYPE| type |

    **MUTATION**

    - **COLUMN_MAPPINGS** (None)

    - **SCRIPTS** 

        ```bash
        - src/data/streams/mutations/fitbit/add_local_date_time.py
        - src/data/streams/mutations/fitbit/add_zero_timestamp.py
        ```

    !!! note
        `add_zero_timestamp` adds an all-zero column called `timestamp` that will be filled in later in the pipeline by `readable_time.R` converting LOCAL_DATE_TIME to a unix timestamp taking into account single or multiple time zones.

        Fitbit API has two versions for sleep data, v1 and v1.2. We support both but ignore v1's `count_awake`, `duration_awake`, and `count_awakenings`, `count_restless`, `duration_restless` columns.
        
        ??? example "Example of the expected raw data"

            |device_id                              |local_start_date_time  |local_end_date_time    |efficiency  |minutes_after_wakeup  |minutes_asleep  |minutes_awake  |minutes_to_fall_asleep  |minutes_in_bed  |is_main_sleep  |type     |
            |-------------------------------------- |---------------------- |---------------------- |----------- |--------------------- |--------------- |-------------- |----------------------- |--------------- |-------------- |-------- |
            |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-10 15:36:30    |2020-10-10 16:37:00    |92          |0                     |55              |5              |0                       |60              |0              |classic  |
            |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-10 01:46:30    |2020-10-10 08:10:00    |88          |0                     |318             |65             |0                       |383             |1              |stages   |
            |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-11 00:12:30    |2020-10-11 11:47:00    |89          |1                     |562             |132            |0                       |694             |1              |stages   |
            |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-12 01:31:00    |2020-10-12 09:34:30    |93          |0                     |415             |68             |0                       |483             |1              |stages   |


??? info "FITBIT_SLEEP_INTRADAY"

    **RAPIDS_COLUMN_MAPPINGS**

    | RAPIDS column   | Stream column   |
    |-----------------|-----------------|
    | TIMESTAMP | FLAG_TO_MUTATE |
    | LOCAL_DATE_TIME | local_date_time |
    | DEVICE_ID | device_id |
    | TYPE_EPISODE_ID | type_episode_id |
    | DURATION | duration |
    | IS_MAIN_SLEEP | is_main_sleep |
    | TYPE | type |
    | LEVEL | level |

    **MUTATION**

    - **COLUMN_MAPPINGS** (None)

    - **SCRIPTS** 

        ```bash
        src/data/streams/mutations/fitbit/add_zero_timestamp.py
        ```

    !!! note
        `add_zero_timestamp` adds an all-zero column called `timestamp` that will be filled in later in the pipeline by `readable_time.R` converting LOCAL_DATE_TIME to a unix timestamp taking into account single or multiple time zones.

        Fitbit API has two versions for sleep data, v1 and v1.2, we support both.

        ??? example "Example of the expected raw data"

            |device_id                              |type_episode_id  |local_date_time     |duration  |level      |is_main_sleep  |type           |
            |------------------------------------   |---------------- |------------------- |--------- |---------- |-------------- |-------------- |
            |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |0                |2020-10-10 15:36:30 |60        |restless   |0              |classic        |
            |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |0                |2020-10-10 15:37:30 |660       |asleep     |0              |classic        |
            |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |0                |2020-10-10 15:48:30 |60        |restless   |0              |classic        |
            |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |...              |...                 |...       |...        |...            |...            |
            |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |1                |2020-10-10 01:46:30 |420       |light      |1              |stages         |
            |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |1                |2020-10-10 01:53:30 |1230      |deep       |1              |stages         |

??? info "FITBIT_STEPS_SUMMARY"

    **RAPIDS_COLUMN_MAPPINGS**

    | RAPIDS column   | Stream column   |
    |-----------------|-----------------|
    | TIMESTAMP | FLAG_TO_MUTATE |
    | DEVICE_ID | device_id |
    | LOCAL_DATE_TIME | local_date_time |
    | STEPS | steps |

    **MUTATION**

    - **COLUMN_MAPPINGS** (None)

    - **SCRIPTS** 
    
        ```bash
        src/data/streams/mutations/fitbit/add_zero_timestamp.py
        ```

    !!! note
        `add_zero_timestamp` adds an all-zero column called `timestamp` that will be filled in later in the pipeline by `readable_time.R` converting LOCAL_DATE_TIME to a unix timestamp taking into account single or multiple time zones.

        ??? example "Example of the expected raw data"

            |device_id                              |local_date_time        |steps     |
            |-------------------------------------- |---------------------- |--------- |
            |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-07             |1775      |
            |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-08             |3201      |
            |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-09             |998       |

??? info "FITBIT_STEPS_INTRADAY"

    **RAPIDS_COLUMN_MAPPINGS**

    | RAPIDS column   | Stream column   |
    |-----------------|-----------------|
    | TIMESTAMP | FLAG_TO_MUTATE |
    | DEVICE_ID | device_id |
    | LOCAL_DATE_TIME | local_date_time |
    | STEPS | steps |

    **MUTATION**

    - **COLUMN_MAPPINGS** (None)

    - **SCRIPTS** 

        ```bash
        src/data/streams/mutations/fitbit/add_zero_timestamp.py
        ```

    !!! note
        `add_zero_timestamp` adds an all-zero column called `timestamp` that will be filled in later in the pipeline by `readable_time.R` converting LOCAL_DATE_TIME to a unix timestamp taking into account single or multiple time zones.

        ??? example "Example of the expected raw data"

            |device_id                              |local_date_time        |steps     |
            |-------------------------------------- |---------------------- |--------- |
            |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-07 00:00:00    |5         |
            |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-07 00:01:00    |3         |
            |a748ee1a-1d0b-4ae9-9074-279a2b6ba524   |2020-10-07 00:02:00    |0         |
        
