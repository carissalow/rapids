Minimal Working Example
=======================

This is a quick guide for creating and running a simple pipeline to extract missing, outgoing, and incoming call features for `daily` and `night` epochs of one participant monitored on the US East coast.

1. Install RAPIDS and make sure your `conda` environment is active (see [Installation](../../setup/installation))
2. Make the changes listed below for the corresponding [Initial Configuration](../../setup/configuration) step (we provide an example of what the relevant sections in your `config.yml` will look like after you are done)
    
    !!! info "Things to change on each configuration step"
        1\. Setup your database connection credentials in `.env`. We assume your credentials group is called `MY_GROUP`.

        2\. `America/New_York` should be the default timezone

        3\. Create a participant file `p01.yaml` based on one of your participants and add `p01` to `[PIDS]` in `config.yaml`. The following would be the content of your `p01.yaml` participant file:
            ```yaml
            PHONE:
                DEVICE_IDS: [aaaaaaaa-1111-bbbb-2222-cccccccccccc] # your participant's AWARE device id
                PLATFORMS: [android] # or ios
                LABEL: MyTestP01 # any string
                START_DATE: 2020-01-01 # this can also be empty
                END_DATE: 2021-01-01 # this can also be empty
            ```
        
        4\. `[DAY_SEGMENTS][TYPE]` should be the default `PERIODIC`. Change `[DAY_SEGMENTS][FILE]` with the path of a file containing the following lines:
             ```csv
             label,start_time,length,repeats_on,repeats_value
             daily,00:00:00,23H 59M 59S,every_day,0
             night,00:00:00,5H 59M 59S,every_day,0
             ```

         5\. If you collected data with AWARE you won't need to modify the attributes of `[DEVICE_DATA][PHONE]`

         6\. Set `[PHONE_CALLS][PROVIDERS][RAPIDS][COMPUTE]` to `True`


    !!! example "Example of the `config.yaml` sections after the changes outlined above"
        ```
        PIDS: [p01]

        TIMEZONE: &timezone
        America/New_York

        DATABASE_GROUP: &database_group
        MY_GROUP

        # ... other irrelevant sections

        DAY_SEGMENTS: &day_segments
            TYPE: PERIODIC
            FILE: "data/external/daysegments_periodic.csv" # make sure the three lines specified above are in the file
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
4. The call features for daily and morning day segments will be in 
   ```
   /data/processed/features/p01/phone_calls.csv
   ```


