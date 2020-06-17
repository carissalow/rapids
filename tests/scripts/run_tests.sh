#!/bin/bash
# Commmands necessary to setup and run the tests for RAPIDS

echo Setting up for testing...

echo Copying files...
cp -r tests/data/raw/* data/raw
cp tests/data/external/* data/external

echo Disabling downloading of dataset...
sed -e '10,20 s/^/#/' -e  's/rules.download_dataset.output/"data\/raw\/\{pid\}\/\{sensor\}_raw\.csv"/' rules/preprocessing.snakefile > tmp
cp tmp rules/preprocessing.snakefile

echo Disabling downloading of dataset...
snakemake --profile tests/settings

echo Running tests on data produced...
python -m unittest discover tests/scripts/ -v 