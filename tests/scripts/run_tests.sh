#!/bin/bash
# Commmands necessary to setup and run the tests for RAPIDS

echo Setting up for testing...

# Uncomment the section below if neccessary to remove old files when testing locally
# echo deleting old data...
# rm -rf data/raw/*
# rm -rf data/processed/*
# rm -rf data/interim/*
# rm -rf data/external/test*

echo Copying files... 
cp -r tests/data/raw/* data/raw
cp tests/data/external/* data/external

# Uncomment the section below to backup snakemake file when testing locally
# echo Backing up preprocessing...
# cp rules/preprocessing.smk bak

echo Disabling downloading of dataset...
sed -e '27,39 s/^/#/' -e  's/rules.download_dataset.output/"data\/raw\/\{pid\}\/\{sensor\}_raw\.csv"/' rules/preprocessing.smk > tmp
cp tmp rules/preprocessing.smk

echo Running RAPIDS Pipeline periodic segment on testdata...
snakemake --profile tests/settings/periodic/ 

echo Moving produced data from previous pipeline run ...
# rm -rf data/raw/*
mkdir data/processed/features/periodic
mv data/processed/features/test* data/processed/features/periodic/
rm -rf data/interim/*
# rm -rf data/external/test*

echo Running RAPIDS Pipeline frequnecy segment on testdata...
snakemake --profile tests/settings/frequency/ 

echo Moving produced data from previous pipeline run...
mkdir data/processed/features/frequency
mv data/processed/features/test* data/processed/features/frequency/

echo Running tests on periodic data produced...
python -m unittest discover tests/scripts/ -v 

echo Backing up Testing script...
cp tests/scripts/test_sensor_features.py test_bak

echo Re-writing the config file being loaded for testing
sed -e  's/tests\/settings\/periodic\/testing_config\.yaml/tests\/settings\/frequency\/testing_config\.yaml/' tests/scripts/test_sensor_features.py > test_tmp
cp test_tmp tests/scripts/test_sensor_features.py

echo Running tests on frequency data produced...
python -m unittest discover tests/scripts/ -v

# Uncomment to return snakemake back to the original version when testing locally
# echo Cleaning up...
# mv bak rules/preprocessing.smk
# mv test_bak tests/scripts/test_sensor_features.py
# rm test_bak 
# rm test_tmp
# rm bak
# rm tmp
