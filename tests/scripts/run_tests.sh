#!/bin/bash
# Commmands necessary to setup and run the tests for RAPIDS

cp -r tests/data/raw/* data/raw
cp tests/data/external/* data/external 
sed -e '10,20 s/^/#/' -e  's/rules.download_dataset.output/"data\/raw\/\{pid\}\/\{sensor\}_raw\.csv"/' rules/preprocessing.snakefile > tmp
cp tmp rules/preprocessing.snakefile
snakemake --profile tests/settings
python -m unittest discover tests/scripts/ -v 