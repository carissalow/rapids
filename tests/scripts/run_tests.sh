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
# cp rules/preprocessing.snakefile bak

echo Disabling downloading of dataset...
sed -e '46,58 s/^/#/' -e  's/rules.download_dataset.output/"data\/raw\/\{pid\}\/\{sensor\}_raw\.csv"/' rules/preprocessing.snakefile > tmp
cp tmp rules/preprocessing.snakefile

echo Running RAPIDS Pipeline on testdata...
snakemake --profile tests/settings 

echo Running tests on data produced...
python -m unittest discover tests/scripts/ -v 

# Uncomment to return snakemake back to the original version when testing locally
# echo Cleaning up...
# mv bak rules/preprocessing.snakefile
# rm tmp