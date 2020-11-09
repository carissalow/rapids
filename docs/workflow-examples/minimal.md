Minimal Working Example
=======================

This is a quick guide for creating and running a simple pipeline to extract missing, outgoing, and incoming call features for `daily` and `night` epochs of one participant monitored on the US East coast.

1. Install RAPIDS and make sure your `conda` environment is active (see [Installation](../../setup/installation))
2. For the [Initial Configuration](../../setup/configuration) steps do the following and use the example as a guide:
    
    !!! info "Things to change on each configuration step"
        1\. Setup your database connection credentials in `.env`. We assume your credentials group is called `MY_GROUP`.

        2\. `America/New_York` should be the default timezone

        3\. Create a participant file `p01.yaml` based on one of your participants and add `p01` to `[PIDS]` in `config.yaml`
        
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

        ....

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


