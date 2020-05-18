Testing 
==========

The following is a simple guide to testing RAPIDS. All files necessary for testing are stored in the ``tests`` directory:

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

To begin testing  RAPIDS place the input data ``csv`` files in the ``tests/data/raw`` directory. The expected output of RAPIDS with the raw input data should be placed in the ``tests/data/processesd``. 

Copy all files from ``tests/data/raw`` directory into the ``data/raw`` directory. The rule(s) that are to be tested must be placed in the ``tests/Snakemake`` file. The current ``tests/Snakemake`` is a good example of how to define the rules that are to be tested. 

Store your test scripts in the ``tests/scripts`` directory. Next, you can run all rules in the ``tests/Snakemake`` with:

::

    snakemake --profile tests/settings

Or run a single rule with

:: 

    snakemake --profile tests/settings -R sms_features

The above example runs the ``sms_features`` rule that is defined in the ``tests/Snakemake`` file. Replace this with the name of the rule you want to test. The ``--profile`` flag is used to run ``Snakemake`` with the ``Snakfile`` and ``confi.yaml`` file stored in ``tests/settings``. 

Once RAPIDS has processed the sample data, the next step is to test the output. Testing is implemented using Python's Unittest. To run all the tests scripts stored in the ``tests/scripts`` directory use the following command:

::

    python -m unittest discover tests/scripts/ -v

The ``discover`` flag finds and runs all of the test scripts within the ``tests/scripts`` directory that start with ``test_``. The name of all test methods in the these scripts should also start with ``test_``.

The following is a snippet of the output you should see after running your test. 

::

    test_sensors_files_exist (test_sensor_features.TestSensorFeatures) ... ok
    test_sensors_features_calculations (test_sensor_features.TestSensorFeatures) ... FAIL

    ======================================================================
    FAIL: test_sensors_features_calculations (test_sensor_features.TestSensorFeatures)
    ----------------------------------------------------------------------

The results above show that the first test ``test_sensors_files_exist`` passed while ``test_sensors_features_calculations`` failed. In addition you should get the traceback of the failure (not shown here). For more information on how to implement test scripts and use unittest please see `Unittest Documentation`_

.. _`Unittest Documentation`: https://docs.python.org/3.7/library/unittest.html#command-line-interface
