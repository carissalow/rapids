# Testing

The following is a simple guide to testing RAPIDS. All files necessary for testing are stored in the `/tests` directory

## Steps for Testing

1.  To begin testing RAPIDS place the fake raw input data `csv` files in
    `tests/data/raw/`. The fake participant files should be placed in
    `tests/data/external/`. The expected output files of RAPIDS after
    processing the input data should be placed in
    `tests/data/processesd/`.
2.  The Snakemake rule(s) that are to be tested must be placed in the
    `tests/Snakemake` file. The current `tests/Snakemake` is a good
    example of how to define them. (At the time of writing this
    documentation the snakefile contains rules messages (SMS), calls and
    screen)
3.  Edit the `tests/settings/config.yaml`. Add and/or remove the rules
    to be run for testing from the `forcerun` list.
4.  Edit the `tests/settings/testing_config.yaml` with the necessary
    configuration settings for running the rules to be tested.
5.  Add any additional testscripts in `tests/scripts`.
6.  Uncomment or comment off lines in the testing shell script
    `tests/scripts/run_tests.sh`.
7.  Run the testing shell script.

    ```bash
    tests/scripts/run_tests.sh
    ```

The following is a snippet of the output you should see after running your test.

```bash
test_sensors_files_exist (test_sensor_features.TestSensorFeatures) ... ok
test_sensors_features_calculations (test_sensor_features.TestSensorFeatures) ... FAIL

======================================================================
FAIL: test_sensors_features_calculations (test_sensor_features.TestSensorFeatures)
----------------------------------------------------------------------
```

The results above show that the first test `test_sensors_files_exist` passed while `test_sensors_features_calculations` failed. In addition you should get the traceback of the failure (not shown here). For more information on how to implement test scripts and use unittest please see [Unittest Documentation](https://docs.python.org/3.7/library/unittest.html#command-line-interface)

Testing of the RAPIDS sensors and features is a work-in-progress. Please see `test-cases`{.interpreted-text role="ref"} for a list of sensors and features that have testing currently available.

Currently the repository is set up to test a number of sensors out of the box by simply running the `tests/scripts/run_tests.sh` command once the RAPIDS python environment is active.
