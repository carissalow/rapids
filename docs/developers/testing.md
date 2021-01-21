# Testing

The following is a simple guide to testing RAPIDS. All files necessary for testing are stored in the `/tests` directory

## Steps for Testing

1.  To begin testing RAPIDS place the fake raw input data `csv` files of each fake participant in
    `tests/data/raw/`. The fake participant files should be placed in
    `tests/data/external/participant_files`. The expected output files of RAPIDS after
    processing the input data should be placed in `tests/data/processesd/frequency` and `tests/data/processesd/periodic` for frequency and periodic respectively.
2.  Edit `tests/settings/frequency/config.yaml` and `tests/settings/periodic/config.yaml` to add and/or remove the rules
    to be run for testing from the `forcerun` list.
3.  Edit `tests/settings/frequency/testing_config.yaml` and `tests/settings/frequency/testing_config.yaml` to configure the        settings and enable/disable sensors to be tested.
4.  Add any additional testscripts in `tests/scripts`.
5.  Run the testing shell script with
    ```bash
    tests/scripts/run_tests.sh
    run_test.sh [-l] [all | periodic | frequency] [test]
    ```
    `[-l]` will delete all the existing files in `/data` before running tests.
    `[all | periodic | frequency]` will generate feature data for all or specific type of features and save in `data/processed`.
    `[test]` will compare the features generated with the precomputed and verified features in `/tests/data/processed`.

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

The results above show that the for periodic both `test_sensors_files_exist` and `test_sensors_features_calculations` passed while for frequency first test `test_sensors_files_exist` passed while `test_sensors_features_calculations` failed. In addition you should get the traceback of the failure (not shown here). For more information on how to implement test scripts and use unittest please see [Unittest Documentation](https://docs.python.org/3.7/library/unittest.html#command-line-interface)

Testing of the RAPIDS sensors and features is a work-in-progress. Please see `test-cases`{.interpreted-text role="ref"} for a list of sensors and features that have testing currently available.

Currently the repository is set up to test a number of sensors out of the box by simply running the `tests/scripts/run_tests.sh` command once the RAPIDS python environment is active.
