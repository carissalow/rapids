# Testing

The following is a simple guide to run RAPIDS' tests. All files necessary for testing are stored in the `./tests/` directory

## Steps for Testing

1. **Add raw data.**
    1. Add the raw data to the corresponding sensor CSV file in `tests/data/external/aware_csv`. Create the CSV if it does not exist.
2. **Link raw data.**
    1. Make sure that you link the new raw data to a participant by using the same `device_id` in the data and in `[DEVICE_IDS]` inside their participant file (`tests/data/external/participant_files/testXX.yaml`). 
    2. Create the participant file if it does not exist, and don't forget to edit `[PIDS]` in the config file of the time segments you are testing (see below). For simplicity, we use a participant's id (`testXX`) as their `device_id`.
3. **Edit the config file.**
    1. Activate the sensor provider you are testing if it isn't already. Set `[SENSOR][PROVIDER][COMPUTE]` to `TRUE` in the `config.yaml` of the time segments you are testing:
    ```yaml
    - tests/settings/frequency_config.yaml # For frequency time segments
    - tests/settings/periodic_config.yaml # For periodic time segments
    # We have not tested events time segments yet
    ```
4. **Run the pipeline and tests.**
    1. You can run all time segments pipelines and their tests
    ```bash
    tests/scripts/run_tests.sh -t all
    ```
    2. You can run only the pipeline of a specific time segment and its tests
    ```bash
    tests/scripts/run_tests.sh -t frequency -a both
    ```
    2. Or, if you are working on your tests and you want to run a pipeline and its tests independently
    ```bash
    tests/scripts/run_tests.sh -t frequency -a run
    tests/scripts/run_tests.sh -t frequency -a test
    ```

## Output example
The following is a snippet of the output you should see after running your test.

```bash
test_sensors_files_exist (test_sensor_features.TestSensorFeatures) ... periodic
ok
test_sensors_features_calculations (test_sensor_features.TestSensorFeatures) ... periodic
ok

test_sensors_files_exist (test_sensor_features.TestSensorFeatures) ... frequency
ok
test_sensors_features_calculations (test_sensor_features.TestSensorFeatures) ... frequency
FAIL
```

The results above show that the for periodic both `test_sensors_files_exist` and `test_sensors_features_calculations` passed while for frequency first test `test_sensors_files_exist` passed while `test_sensors_features_calculations` failed. Additionally, you should get the traceback of the failure (not shown here). For more information on how to implement test scripts and use unittest please see [Unittest Documentation](https://docs.python.org/3.7/library/unittest.html#command-line-interface)

Testing of the RAPIDS sensors and features is a work-in-progress. Please see [Test Cases](../test-cases) for a list of sensors and features that have testing currently available.

## How do we execute the tests?
This bash script `tests/scripts/run_tests.sh` executes one or all pipelines for different time segment types (`frequency`, `periodic`, and `events`) as well as their tests (see below).

This python script `tests/scripts/run_tests.py` runs the tests. It parses the involved participants and active sensor providers in the `config.yaml` file of the time segment type being tested. We test that the output file we expect exists and that its content matches the expected values.

??? example "Example of raw data for PHONE_APPLICATIONS_FOREGROUND testing"
    ```json hl_lines="1 2 4" linenums="1"
    --8<---- "tests/data/external/aware_csv/phone_applications_foreground_raw.csv"
    ```

## What cases do we test?
The sample data includes 7 tests cases. Take phone battery as an example, on this platform, battery status 2 represents `charging` and battery status 4 represents `discharge`. 

??? "1. A daily segment instance with no battery episodes"

    ??? "Example"

        Input time segments:

        | timestamp | device_id | battery_status | battery_level | battery_scale | battery_voltage | battery_temperature | battery_adaptor | battery_health | battery_technology |
        |---|---|---|---|---|---|---|---|---|---|
        | 00:17:38.602 | test02 | 4 | 77 | 100 | 4157 | 23 | 0 | 2 | Li-ion |
        | 03:20:30.415 | test02 | 2 | 77 | 100 | 4170 | 23 | 0 | 2 | Li-ion |

        Output results

        | local_segment | local_segment_label | local_segment_start_datetime | local_segment_end_datetime | phone_battery_rapids_countdischarge | phone_battery_rapids_sumdurationdischarge | phone_battery_rapids_avgconsumptionrate | phone_battery_rapids_maxconsumptionrate | phone_battery_rapids_countcharge | phone_battery_rapids_sumdurationcharge |
        | --- |---|---|---|---|---|---|---|---|---|
        | 00:00:00,00:29:59  | thirtyminutes0000 | 2020-07-01 00:00:00 | 2020-07-01 00:29:59 | 1 | 21.8259833333333 | 0.137450851775292 | 0.137450851775292 | 0 | 0 |
        | 00:03:00,03:29:59  | thirtyminutes0006 | 2020-07-01 03:00:00 | 2020-07-01 03:29:59 | 0 | 0 | 0 | 0 | 1 | 9.49288333333333 |
    
        Since there is no battery episode between 00:00:30 and 03:00:00, no result will be generated for this epoch.

??? "2. A daily segment instance with two battery episodes (one charging, one discharge)"

    ??? "Periodic (daily)"

        Input time segments:

        | timestamp             | device_id | battery_status | battery_level | battery_scale | battery_voltage | battery_temperature | battery_adaptor | battery_health | battery_technology |
        |-----------------------|-----------|----------------|---------------|---------------|-----------------|---------------------|-----------------|----------------|--------------------|
        | 17:59:41.434 | test02    | 4              | 59            | 100           | 4094            | 23                  | 0               | 2              | Li-ion             |
        | 18:04:14.321 | test02    | 4              | 58            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
        | 18:07:24.456 | test02    | 4              | 57            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
        | 20:03:03.415 | test02    | 2              | 72            | 100           | 4170            | 23                  | 0               | 2              | Li-ion             |
        | 20:05:12.434 | test02    | 2              | 73            | 100           | 4094            | 23                  | 0               | 2              | Li-ion             |
        | 20:07:24.678 | test02    | 2              | 74            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
        | 20:10:34.875 | test02    | 2              | 75            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
        | 21:30:04.415 | test02    | 4              | 74            | 100           | 4170            | 23                  | 0               | 2              | Li-ion             |
        | 21:32:14.434 | test02    | 4              | 73            | 100           | 4094            | 23                  | 0               | 2              | Li-ion             |
        | 21:35:23.678 | test02    | 4              | 72            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
        | 21:37:47.875 | test02    | 4              | 71            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |

        Output results:

        | local_segment | local_segment_label | local_segment_start_datetime | local_segment_end_datetime | phone_battery_rapids_countdischarge | phone_battery_rapids_sumdurationdischarge | phone_battery_rapids_avgconsumptionrate | phone_battery_rapids_maxconsumptionrate | phone_battery_rapids_countcharge | phone_battery_rapids_sumdurationcharge |
        | --- |---|---|---|---|---|---|---|---|---|
        | 18:00:00,23:59:59  | evening | 2020-07-01 18:00:00 | 2020-07-01 23:59:59 | 2 | 75.1306166666666 | 0.0664958369201784 | 0.079525673538274 | 1 | 37.5236666666667 |


    ??? "Frequency (30 mins)"
        
        Input time segments:

        | timestamp | device_id | battery_status | battery_level | battery_scale | battery_voltage | battery_temperature | battery_adaptor | battery_health | battery_technology |
        |---|---|---|---|---|---|---|---|---|---|
        | 20:10:34.875 | test06 | 2 | 75 | 100 | 4157 | 23 | 0 | 2 | Li-ion |
        | 20:20:17.171 | test06 | 4 | 74 | 100 | 4170 | 23 | 0 | 2 | Li-ion |

        Output results

        | local_segment | local_segment_label | local_segment_start_datetime | local_segment_end_datetime | phone_battery_rapids_countdischarge | phone_battery_rapids_sumdurationdischarge | phone_battery_rapids_avgconsumptionrate | phone_battery_rapids_maxconsumptionrate | phone_battery_rapids_countcharge | phone_battery_rapids_sumdurationcharge |
        | --- |---|---|---|---|---|---|---|---|---|
        | 20:00:00,20:29:59 | thirtyminutes0040 | 2020-07-01 20:00:00 | 2020-07-01 20:29:59 | 1 | 14.6351666666667 | 0.0683285693136395 | 0.0683285693136395 | 1 | 12.3074 |


??? "3. A daily segment instance with a charging episode that spans to the next daily instance"

    ??? "Periodic (daily)"

        Input time segments:

        | timestamp | device_id | battery_status | battery_level | battery_scale | battery_voltage | battery_temperature | battery_adaptor | battery_health | battery_technology |
        |---|---|---|---|---|---|---|---|---|---|
        | 11:59:28.434 | test02 | 2 | 63 | 100 | 4094 | 23 | 0 | 2 | Li-ion |
        | 12:04:37.678 | test02 | 2 | 64 | 100 | 4157 | 23 | 0 | 2 | Li-ion |

    ??? "Frequency (30 mins)"

        Input time segements:

        | timestamp | device_id | battery_status | battery_level | battery_scale | battery_voltage | battery_temperature | battery_adaptor | battery_health | battery_technology |
        |---|---|---|---|---|---|---|---|---|---|
        | 11:59:28.434 | test06 | 2 | 63 | 100 | 4094 | 23 | 0 | 2 | Li-ion |
        | 12:04:37.678 | test06 | 2 | 64 | 100 | 4157 | 23 | 0 | 2 | Li-ion |

??? "4. A daily segment instance with a discharge episode that spans to the next daily instance"

    ??? "Periodic (daily)"

        Input time segements:

        | timestamp | device_id | battery_status | battery_level | battery_scale | battery_voltage | battery_temperature | battery_adaptor | battery_health | battery_technology |
        |---|---|---|---|---|---|---|---|---|---|
        | 05:59:49.434 | test02 | 4 | 79 | 100 | 4094 | 23 | 0 | 2 | Li-ion |
        | 06:02:19.321 | test02 | 4 | 78 | 100 | 4157 | 23 | 0 | 2 | Li-ion |

    ??? "Frequency (30 mins)"

        Input time segements:

        | timestamp | device_id | battery_status | battery_level | battery_scale | battery_voltage | battery_temperature | battery_adaptor | battery_health | battery_technology |
        |---|---|---|---|---|---|---|---|---|---|
        | 17:59:41.434 | test06 | 4 | 59 | 100 | 4094 | 23 | 0 | 2 | Li-ion |
        | 18:04:14.321 | test06 | 4 | 58 | 100 | 4157 | 23 | 0 | 2 | Li-ion |
    
??? "5. Three-day segments that repeat everyday"

    [Time segment tested:](../setup/configuration.md#time-segments)

    | label | start_time | length | repeats_on | repeats_value |
    |---|---|---|---|---|
    | daily | 00:00:00 | 23H 59M 59S | every_day | 0 |


    Data tested:

    We test 14 segments, one at the beginning of the first day, one at the end of the last day

    | timestamp             | device_id | battery_status | battery_level | battery_scale | battery_voltage | battery_temperature | battery_adaptor | battery_health | battery_technology |
    |-----------------------|-----------|----------------|---------------|---------------|-----------------|---------------------|-----------------|----------------|--------------------|
    | 2020-07-02 00:03:47.875 | test01    | 3              | 63            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-07-02 00:05:47.875 | test01    | 3              | 62            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-07-02 23:55:47.875 | test01    | 3              | 55            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-07-02 23:59:47.875 | test01    | 3              | 54            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-07-03 00:06:47.875 | test01    | 3              | 53            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-07-03 00:09:47.875 | test01    | 3              | 52            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-07-03 23:47:05.000 | test01    | 3              | 60            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-07-03 23:55:05.000 | test01    | 3              | 59            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-07-04 00:15:05.000 | test01    | 3              | 58            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-07-04 00:18:05.000 | test01    | 3              | 57            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-07-04 23:51:00.000 | test01    | 3              | 41            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-07-04 23:57:00.000 | test01    | 3              | 40            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-07-05 00:21:00.000 | test01    | 3              | 39            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-07-05 00:23:00.000 | test01    | 3              | 38            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |


    Output results:

    | local_segment                                    | local_segment_label | local_segment_start_datetime | local_segment_end_datetime | phone_battery_rapids_countdischarge | phone_battery_rapids_sumdurationdischarge | phone_battery_rapids_avgconsumptionrate | phone_battery_rapids_maxconsumptionrate | phone_battery_rapids_countcharge | phone_battery_rapids_sumdurationcharge |
    |--------------------------------------------------|---------------------|------------------------------|----------------------------|-------------------------------------|-------------------------------------------|-----------------------------------------|-----------------------------------------|----------------------------------|----------------------------------------|
    | threeday#2020-07-02 00:00:00,2020-07-04 23:59:59 | threeday            | 2020-07-02 00:00:00          | 2020-07-04 23:59:59        | 4                                   | 149.7954                                  | 0.0710868450815781                      | 0.111113168762384                       | 0                                | 0                                      |
    | threeday#2020-07-03 00:00:00,2020-07-05 23:59:59 | threeday            | 2020-07-03 00:00:00          | 2020-07-05 23:59:59        | 3                                   | 162.7952                                  | 0.0492745931499224                      | 0.0502547286558745                      | 0                                | 0                                      |
    | threeday#2020-07-04 00:00:00,2020-07-06 23:59:59 | threeday            | 2020-07-04 00:00:00          | 2020-07-06 23:59:59        | 2                                   | 110.0815                                  | 0.0449915246814979                      | 0.0483879032392475                      | 0                                | 0                                      |
    | threeday#2020-07-05 00:00:00,2020-07-07 23:59:59 | threeday            | 2020-07-05 00:00:00          | 2020-07-07 23:59:59        | 1                                   | 52.9991166666667                          | 0.0377364779979038                      | 0.0377364779979038                      | 0                                | 0                                      |

??? "6. A three-day segment that repeats on a fixed day"
    
    [Time segment tested:](../setup/configuration.md#time-segments)

    | label | start_time | length | repeats_on | repeats_value |
    |---|---|---|---|---|
    | weekends | 00:00:00 | 2D 23H 59M 59S | wday | 5 |

    Data tested:

    We test 10 segments, one at the beginning of the first day, one at the end of the last day

    | timestamp             | device_id | battery_status | battery_level | battery_scale | battery_voltage | battery_temperature | battery_adaptor | battery_health | battery_technology |
    |-----------------------|-----------|----------------|---------------|---------------|-----------------|---------------------|-----------------|----------------|--------------------|
    | 2020-07-03 00:06:47.875 | test01    | 3              | 53            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-07-03 00:09:47.875 | test01    | 3              | 52            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-07-03 23:47:05.000 | test01    | 3              | 60            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-07-03 23:55:05.000 | test01    | 3              | 59            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-07-04 00:15:05.000 | test01    | 3              | 58            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-07-04 00:18:05.000 | test01    | 3              | 57            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-07-04 23:51:00.000 | test01    | 3              | 41            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-07-04 23:57:00.000 | test01    | 3              | 40            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-07-05 00:21:00.000 | test01    | 3              | 39            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-07-05 00:23:00.000 | test01    | 3              | 38            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |

    Output results:

    | local_segment                                    | local_segment_label | local_segment_start_datetime | local_segment_end_datetime | phone_battery_rapids_countdischarge | phone_battery_rapids_sumdurationdischarge | phone_battery_rapids_avgconsumptionrate | phone_battery_rapids_maxconsumptionrate | phone_battery_rapids_countcharge | phone_battery_rapids_sumdurationcharge |
    |--------------------------------------------------|---------------------|------------------------------|----------------------------|-------------------------------------|-------------------------------------------|-----------------------------------------|-----------------------------------------|----------------------------------|----------------------------------------|
    | weekends#2020-07-03 00:00:00,2020-07-05 23:59:59 | weekends            | 2020-07-03 00:00:00          | 2020-07-05 23:59:59        | 3                                   | 162.7952                                  | 0.0492745931499224                      | 0.0502547286558745                      | 0                                | 0                                      |

??? "7. Event segements"

    [Time segments tested:](../setup/configuration.md#time-segments)

    | label | event_timestamp | length | shift | shift_direction | device_id |
    |---|---|---|---|---|---|
    | survey1 | 1587661220000 | 10H | 10H | -1 | a748ee1a-1d0b-4ae9-9074-279a2b6ba524 |
    | survey2 | 1587661220000 | 10H | 5H | -1 | a748ee1a-1d0b-4ae9-9074-279a2b6ba524 |
    | survey3 | 1587661220000 | 10H | 0H | 1 | a748ee1a-1d0b-4ae9-9074-279a2b6ba524 |

    Data tested: 

    We test 7 segments, one at the beginning of the first day, one at the end of the last day

    | timestamp             | device_id                            | battery_status | battery_level | battery_scale | battery_voltage | battery_temperature | battery_adaptor | battery_health | battery_technology |
    |-----------------------|--------------------------------------|----------------|---------------|---------------|-----------------|---------------------|-----------------|----------------|--------------------|
    | 2020-04-23 03:15:00.000 | a748ee1a-1d0b-4ae9-9074-279a2b6ba524 | 3              | 90            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-04-23 03:21:00.000 | a748ee1a-1d0b-4ae9-9074-279a2b6ba524 | 3              | 89            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-04-23 07:50:00.000 | a748ee1a-1d0b-4ae9-9074-279a2b6ba524 | 3              | 80            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-04-23 08:05:00.000 | a748ee1a-1d0b-4ae9-9074-279a2b6ba524 | 3              | 79            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-04-23 08:12:00.000 | a748ee1a-1d0b-4ae9-9074-279a2b6ba524 | 3              | 78            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-04-23 22:50:00.000 | a748ee1a-1d0b-4ae9-9074-279a2b6ba524 | 3              | 50            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |
    | 2020-04-23 22:53:00.000 | a748ee1a-1d0b-4ae9-9074-279a2b6ba524 | 3              | 49            | 100           | 4157            | 23                  | 0               | 2              | Li-ion             |

    Output results:

    | local_segment                                   | local_segment_label | local_segment_start_datetime | local_segment_end_datetime | phone_battery_rapids_sumdurationcharge | phone_battery_rapids_countdischarge | phone_battery_rapids_sumdurationdischarge | phone_battery_rapids_maxconsumptionrate | phone_battery_rapids_avgconsumptionrate | phone_battery_rapids_countcharge |
    |-------------------------------------------------|---------------------|------------------------------|----------------------------|----------------------------------------|-------------------------------------|-------------------------------------------|-----------------------------------------|-----------------------------------------|----------------------------------|
    | survey1#2020-04-23 03:00:20,2020-04-23 13:00:20 | survey1             | 2020-04-23 03:00:20          | 2020-04-23 13:00:20        | 0                                      | 2                                   | 87.9985333333333                          | 0.0384621794978634                      | 0.0331202101231602                      | 0                                |
    | survey2#2020-04-23 08:00:20,2020-04-23 18:00:20 | survey2             | 2020-04-23 08:00:20          | 2020-04-23 18:00:20        | 0                                      | 1                                   | 41.6659833333333                          | 0.0480007872129103                      | 0.0480007872129103                      | 0                                |
    | survey3#2020-04-23 13:00:20,2020-04-23 23:00:20 | survey3             | 2020-04-23 13:00:20          | 2020-04-23 23:00:20        | 0                                      | 1                                   | 10.3498                                   | 0.0966202245454018                      | 0.0966202245454018                      | 0                                |
