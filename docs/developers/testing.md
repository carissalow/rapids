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