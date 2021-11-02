Minimal Working Example
=======================

This is a quick guide for creating and running a simple pipeline to extract missing, outgoing, and incoming `call` features for `24 hr` (`00:00:00` to `23:59:59`) and `night` (`00:00:00` to `05:59:59`) time segments of every day of data of one participant that was monitored on the US East coast with an Android smartphone.

1. Install RAPIDS and make sure your `conda` environment is active (see [Installation](../../setup/installation))
3. Download this [CSV file](../img/calls.csv) and save it as `data/external/aware_csv/calls.csv`
2. Make the changes listed below for the corresponding [Configuration](../../setup/configuration) step (we provide an example of what the relevant sections in your `config.yml` will look like after you are done)
    
    ??? info "Required configuration changes (*click to expand*)"
        1. **Supported [data streams](../../setup/configuration#supported-data-streams).** 
            
            Based on the docs, we decided to use the `aware_csv` data stream because we are processing aware data saved in a CSV file. We will use this label in a later step; there's no need to type it or save it anywhere yet.

        3. **Create your [participants file](../../setup/configuration#participant-files).**
        
            Since we are processing data from a single participant, you only need to create a single participant file called `p01.yaml` in `data/external/participant_files`. This participant file only has a `PHONE` section because this hypothetical participant was only monitored with a smartphone. Note that for a real analysis, you can do this [automatically with a CSV file](../../setup/configuration##automatic-creation-of-participant-files)
            
            1. Add `p01` to `[PIDS]` in `config.yaml`

            1. Create a file in `data/external/participant_files/p01.yaml` with the following content:

                ```yaml
                PHONE:
                    DEVICE_IDS: [a748ee1a-1d0b-4ae9-9074-279a2b6ba524] # the participant's AWARE device id
                    PLATFORMS: [android] # or ios
                    LABEL: MyTestP01 # any string
                    START_DATE: 2020-01-01 # this can also be empty
                    END_DATE: 2021-01-01 # this can also be empty
                ```
        
        4. **Select what [time segments](../../setup/configuration#time-segments) you want to extract features on.** 
        
            1. Set `[TIME_SEGMENTS][FILE]` to `data/external/timesegments_periodic.csv` 

            1. Create a file in `data/external/timesegments_periodic.csv` with the following content
            
                ```csv
                label,start_time,length,repeats_on,repeats_value
                daily,00:00:00,23H 59M 59S,every_day,0
                night,00:00:00,5H 59M 59S,every_day,0
                ```
        
        2. **Choose the [timezone of your study](../../setup/configuration#timezone-of-your-study).** 
        
            We will use the default time zone settings since this example is processing data collected on the US East Coast (`America/New_York`)

            ```yaml
            TIMEZONE: 
                TYPE: SINGLE
                SINGLE:
                    TZCODE: America/New_York
            ```

         5. **Modify your [device data stream configuration](../../setup/configuration#data-stream-configuration)**
            
            1. Set `[PHONE_DATA_STREAMS][USE]` to `aware_csv`. 
            
            2. We will use the default value for `[PHONE_DATA_STREAMS][aware_csv][FOLDER]` since we already stored the test calls CSV file there.

         6. **Select what [sensors and features](../../setup/configuration#sensor-and-features-to-process) you want to process.** 
         
            1. Set `[PHONE_CALLS][CONTAINER]` to `calls.csv` in the `config.yaml` file.

            1. Set `[PHONE_CALLS][PROVIDERS][RAPIDS][COMPUTE]` to `True` in the `config.yaml` file.


    !!! example "Example of the `config.yaml` sections after the changes outlined above"

        This will be your `config.yaml` after following the instructions above. Click on the numbered markers to know more.

        ``` { .yaml .annotate } 
        PIDS: [p01] # (1)
        
        TIMEZONE:
            TYPE: SINGLE # (2)
            SINGLE:
                TZCODE: America/New_York

        # ... other irrelevant sections

        TIME_SEGMENTS: &time_segments
            TYPE: PERIODIC # (3)
            FILE: "data/external/timesegments_periodic.csv" # (4)
            INCLUDE_PAST_PERIODIC_SEGMENTS: FALSE

        PHONE_DATA_STREAMS:
            USE: aware_csv # (5)

            aware_csv:
                FOLDER: data/external/aware_csv # (6)

        # ... other irrelevant sections

        ############## PHONE ###########################################################
        ################################################################################

        # ... other irrelevant sections

        # Communication call features config, TYPES and FEATURES keys need to match
        PHONE_CALLS:
            CONTAINER: calls.csv  # (7) 
            PROVIDERS:
                RAPIDS:
                    COMPUTE: True # (8)
                    CALL_TYPES: ...
        ```

        1. We added `p01` to PIDS after creating the participant file:
            ```bash
            data/external/participant_files/p01.yaml
            ```

            With the following content:
            ```yaml
            PHONE:
                DEVICE_IDS: [a748ee1a-1d0b-4ae9-9074-279a2b6ba524] # the participant's AWARE device id
                PLATFORMS: [android] # or ios
                LABEL: MyTestP01 # any string
                START_DATE: 2020-01-01 # this can also be empty
                END_DATE: 2021-01-01 # this can also be empty
            ```

        2. We use the default `SINGLE` time zone.

        3. We use the default `PERIODIC` time segment `[TYPE]`

        4. We created this time segments file with these lines:

            ```csv
            label,start_time,length,repeats_on,repeats_value
            daily,00:00:00,23H 59M 59S,every_day,0
            night,001:00:00,5H 59M 59S,every_day,0
            ```

        5. We set `[USE]` to `aware_device` to tell RAPIDS to process sensor data collected with the AWARE Framework stored in CSV files.

        6. We used the default `[FOLDER]` for `awre_csv` since we already stored our test `calls.csv` file there

        7. We changed `[CONTAINER]` to `calls.csv` to process our test call data.

        8. We flipped `[COMPUTE]` to `True` to extract call behavioral features using the `RAPIDS` feature provider.

3. Run RAPIDS
    ```bash
    ./rapids -j1
    ```
4. The call features for daily and morning time segments will be in 
   ```
   data/processed/features/all_participants/all_sensor_features.csv
   ```


