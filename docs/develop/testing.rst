Testing 
==========

The following is a simple guide to testing RAPIDS. All of documents that are used for testing is stored in the ``tests`` directory. The following is the structure of the ``tests``  directory

::

    ├── tests
    │   ├── data                        <- Replication of the project root data directory for testing.
    │   │   ├── external                <- Contains Data from third party sources used for testing.
    │   │   ├── interim                 <- The expected intermediate data that has been transformed.
    │   │   ├── processed               <- The expected final, canonical data sets for modeling.
    │   │   └── raw                     <- The specially created raw input datasets that will be used for testing.
    │   │   
    │   ├── scripts                     <- Scripts for testing. Add test scripts in this directory.
    │   │   └── utils.py                <- Contains any helper functions and methods.
    │   │
    │   ├── settings                    <- The directory contains the config and settings files for testing snakemake.
    │   │   ├── config.yaml             <- Defines the testing profile configurations for running snakemake.
    │   │   └── testing_config.yaml     <- Contains the actual snakemake configuration settings for testing.
    │   │
    │   └── Snakefile                   <- The Snakefile for testing only. It contains the rules that you would be testing.
    │

To begin testing  RAPIDS place the raw data ``csv`` files in the ``tests/data/raw`` directory. The expected output of applying the rules being tested on the input data should be placed in the ``tests/data/processesd``. In order to test a rule, copy the necessary input data files from the ``tests/data/raw`` directory into the ``data/raw`` directory. The rule(s) that are to be tested must be placed in the ``tests/Snakemake`` file. The current ``tests/Snakemake`` is a good example of how to define the rules that are to be tested. 
Your test scripts are to be placed in the ``tests/scripts`` directory.  The next step is to run the rule that is to be tested. You can run all rules in the ``tests/Snakemake`` file by using the following command 

::

    snakemake --profile tests/settings

Or run a single rule by using a command similar to the following example

:: 

    snakemake --profile tests/settings -R sms_features

The above example runs the ``sms_features`` rule that is defined in the ``tests/Snakemake`` file. This can be changed to the name of the rule that you desire to test. It should be noted that the ``--profile`` option is used in order to specify which configuration setting are to be used in the running of Snakemake. The path to the configurations is passed as an argument for this option. 

Once the rule(s) has been run on the sample data the next step is to test the output. The testing is implemented using Python's Unittest. To run the tests scripts that are stored in the ``tests/scripts`` directory use the following command

::

    python -m unittest discover tests/scripts/ -v

The above command would run all of the test scripts that are in the ``tests/scripts`` directory. The ``discover`` option in the above command finds all of the test scripts within the ``tests/scripts`` directory by matching all of the scripts with that begin with ``test_`` and run them. It should be noted that the test methods in the test scripts also begin with ``test_``. This is how unittest finds the tests to be run. 

The following is an example snippet of the output that you should see after running your test. 

::

    test_sensors_files_exist (test_sensor_features.TestSensorFeatures) ... ok
    test_sensors_features_calculations (test_sensor_features.TestSensorFeatures) ... FAIL

    ======================================================================
    FAIL: test_sensors_features_calculations (test_sensor_features.TestSensorFeatures)
    ----------------------------------------------------------------------

The results above show that the first test ``test_sensors_files_exist`` passed while ``test_sensors_features_calculations`` failed. Note that this is just a snippet of the expected results and the rest of the results (not shown) is the traceback and an explanation of for the failure. For more information on how to implement test scripts and use unittest please see `Unittest Documentation`_

.. _`Unittest Documentation`: https://docs.python.org/3.7/library/unittest.html#command-line-interface
