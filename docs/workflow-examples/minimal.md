Minimal Working Example
=======================

This is a quick guide for creating and running a simple pipeline to extract missing, outgoing, and incoming `call` features for `daily` (`00:00:00` to `23:59:59`) and `night` (`00:00:00` to `05:59:59`) epochs of every day of data of one participant monitored on the US East coast with an Android smartphone.

!!! hint
    If you don't have `call` data that you can use to try this example you can restore this [CSV file](../img/calls.csv) as a table in a MySQL database.


1. Install RAPIDS and make sure your `conda` environment is active (see [Installation](../../setup/installation))
2. Make the changes listed below for the corresponding [Configuration](../../setup/configuration) step (we provide an example of what the relevant sections in your `config.yml` will look like after you are done)
    
    ??? info "Required configuration changes"
        1. **Add your [database credentials](../../setup/configuration#database-credentials).** 
            
            Setup your database connection credentials in `.env`, we assume your credentials group in the `.env` file is called `MY_GROUP`.

        2. **Choose the [timezone of your study](../../setup/configuration#timezone-of-your-study).** 
        
            Since this example is processing data collected on the US East cost, `America/New_York` should be the configured timezone, change this according to your data.

        3. **Create your [participants files](../../setup/configuration#participant-files).**
        
            Since we are processing data from a single participant, you only need to create a single participant file called `p01.yaml`. This participant file only has a `PHONE` section because this hypothetical participant was only monitored with an smartphone. You also need to add `p01` to `[PIDS]` in `config.yaml`. The following would be the content of your `p01.yaml` participant file:
            ```yaml
            PHONE:
                DEVICE_IDS: [a748ee1a-1d0b-4ae9-9074-279a2b6ba524] # the participant's AWARE device id
                PLATFORMS: [android] # or ios
                LABEL: MyTestP01 # any string
                START_DATE: 2020-01-01 # this can also be empty
                END_DATE: 2021-01-01 # this can also be empty
            ```
        
        4. **Select what [time segments](../../setup/configuration#time-segments) you want to extract features on.** 
        
            `[TIME_SEGMENTS][TYPE]` should be the default `PERIODIC`. Change `[TIME_SEGMENTS][FILE]` with the path (for example `data/external/timesegments_periodic.csv`) of a file containing the following lines:
             ```csv
             label,start_time,length,repeats_on,repeats_value
             daily,00:00:00,23H 59M 59S,every_day,0
             night,00:00:00,5H 59M 59S,every_day,0
             ```

         5. **Modify your [device data source configuration](../../setup/configuration#device-data-source-configuration)**
            
            In this example we do not need to modify this section because we are using smartphone data collected with AWARE stored on a MySQL database.

         6. **Select what [sensors and features](../../setup/configuration#sensor-and-features-to-process) you want to process.** 
         
            Set `[PHONE_CALLS][PROVIDERS][RAPIDS][COMPUTE]` to `True` in the `config.yaml` file.


    ??? example "Example of the `config.yaml` sections after the changes outlined above"
        Highlighted lines are related to the configuration steps above.
        ``` yaml hl_lines="1 4 7 12 13 38"
        PIDS: [p01]

        TIMEZONE: &timezone
        America/New_York

        DATABASE_GROUP: &database_group
        MY_GROUP

        # ... other irrelevant sections

        TIME_SEGMENTS: &time_segments
            TYPE: PERIODIC
            FILE: "data/external/timesegments_periodic.csv" # make sure the three lines specified above are in the file
            INCLUDE_PAST_PERIODIC_SEGMENTS: FALSE

        # No need to change this if you collected AWARE data on a database and your credentials are grouped under `MY_GROUP` in `.env`
        DEVICE_DATA:
            PHONE:
                SOURCE: 
                    TYPE: DATABASE
                    DATABASE_GROUP: *database_group
                    DEVICE_ID_COLUMN: device_id # column name
                TIMEZONE: 
                    TYPE: SINGLE # SINGLE or MULTIPLE
                    VALUE: *timezone 


        ############## PHONE ###########################################################
        ################################################################################

        # ... other irrelevant sections

        # Communication call features config, TYPES and FEATURES keys need to match
        PHONE_CALLS:
            TABLE: calls # change if your calls table has a different name
            PROVIDERS:
                RAPIDS:
                    COMPUTE: True # set this to True!
                    CALL_TYPES: ...
        ```

3. Run RAPIDS
    ```bash
    ./rapids -j1
    ```
4. The call features for daily and morning time segments will be in 
   ```
   /data/processed/features/p01/phone_calls.csv
   ```


