Testing 
==========

The following is a simple guide to testing RAPIDS. All files necessary for testing are stored in the ``tests`` directory:

::

    ├── tests
    │   ├── data                        <- Replica of the project root data directory for testing.
    │   │   ├── external                <- Contains the fake testing participant files. 
    │   │   ├── interim                 <- The expected intermediate data that has been transformed.
    │   │   ├── processed               <- The expected final data, canonical data sets for modeling used to test/validate feature calculations.
    │   │   └── raw                     <- The specially created raw input datasets (fake data) that will be used for testing.
    │   │   
    │   ├── scripts                     <- Scripts for testing. Add test scripts in this directory.
    │   │   ├── run_tests.sh            <- The shell script to runs RAPIDS pipeline test data and test the results
    │   │   ├── test_sensor_features.py <- The default test script for testing RAPIDS builting sensor features. 
    │   │   └── utils.py                <- Contains any helper functions and methods.
    │   │
    │   ├── settings                    <- The directory contains the config and settings files for testing snakemake.
    │   │   ├── config.yaml             <- Defines the testing profile configurations for running snakemake.
    │   │   └── testing_config.yaml     <- Contains the actual snakemake configuration settings for testing.
    │   │
    │   └── Snakefile                   <- The Snakefile for testing only. It contains the rules that you would be testing.
    │


Steps for Testing
""""""""""""""""""

#. To begin testing  RAPIDS place the fake raw input data ``csv`` files in ``tests/data/raw/``. The fake participant files should be placed in ``tests/data/external/``. The expected output files of RAPIDS after processing the input data should be placed in ``tests/data/processesd/``. 

#. The Snakemake rule(s) that are to be tested must be placed in the ``tests/Snakemake`` file. The current ``tests/Snakemake`` is a good example of how to define them. (At the time of writing this documentation the snakefile contains rules messages (SMS), calls and screen)

#. Edit the ``tests/settings/config.yaml``. Add and/or remove the rules to be run for testing from the ``forcerun`` list.

#. Edit the ``tests/settings/testing_config.yaml`` with the necessary configuration settings for running the rules to be tested. 

#. Add any additional testscripts in ``tests/scripts``.

#. Uncomment or comment off lines in the testing shell script ``tests/scripts/run_tests.sh``.

#. Run the testing shell script.

::

    $ tests/scripts/run_tests.sh


The following is a snippet of the output you should see after running your test. 

::

    test_sensors_files_exist (test_sensor_features.TestSensorFeatures) ... ok
    test_sensors_features_calculations (test_sensor_features.TestSensorFeatures) ... FAIL

    ======================================================================
    FAIL: test_sensors_features_calculations (test_sensor_features.TestSensorFeatures)
    ----------------------------------------------------------------------

The results above show that the first test ``test_sensors_files_exist`` passed while ``test_sensors_features_calculations`` failed. In addition you should get the traceback of the failure (not shown here). For more information on how to implement test scripts and use unittest please see `Unittest Documentation`_

Testing of the RAPIDS sensors and features is a work-in-progess. Please see :ref:`test-cases` for a list of sensors and features that have testing currently available. 

Currently the repository is set up to test a number of senssors out of the box by simply running the ``tests/scripts/run_tests.sh`` command once the RAPIDS python environment is active. 

.. _`Unittest Documentation`: https://docs.python.org/3.7/library/unittest.html#command-line-interface
