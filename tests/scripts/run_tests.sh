#!/bin/bash
# Commmands necessary to setup and run the tests for RAPIDS

echo Setting up for testing...

clean_old_data() { 

    echo deleting old data...
    rm -rf data/processed/*
    rm -rf data/interim/*

    echo Backing up preprocessing...
    cp rules/preprocessing.smk bak
} 

run_periodic_pipeline() {
    
    echo Running RAPIDS Pipeline periodic segment on testdata...
    snakemake --profile tests/settings/periodic/ 

    echo Moving produced data from previous pipeline run ...
    mkdir data/processed/features/periodic
    mv data/processed/features/test* data/processed/features/periodic/
    rm -rf data/interim/*
}

run_frequency_pipeline() {

    echo Running RAPIDS Pipeline frequency segment on testdata...
    snakemake --profile tests/settings/frequency/ 

    echo Moving produced data from previous pipeline run...
    mkdir data/processed/features/frequency
    mv data/processed/features/test* data/processed/features/frequency/
}

run_periodic_test() {

    echo Running tests on periodic data produced...
    python -m unittest discover tests/scripts/ -v 

    echo Re-writing the config file being loaded for testing
    sed -e  's/tests\/settings\/[a-z]*\/testing_config\.yaml/tests\/settings\/periodic\/testing_config\.yaml/' tests/scripts/test_sensor_features.py > test_tmp
    mv test_tmp tests/scripts/test_sensor_features.py
}

run_frequency_test() {

    # echo Backing up Testing script...
    # cp tests/scripts/test_sensor_features.py test_bak

    echo Re-writing the config file being loaded for testing
    sed -e  's/tests\/settings\/[a-z]*\/testing_config\.yaml/tests\/settings\/frequency\/testing_config\.yaml/' tests/scripts/test_sensor_features.py > test_tmp
    mv test_tmp tests/scripts/test_sensor_features.py

    echo Running tests on frequency data produced...
    python -m unittest discover tests/scripts/ -v
}

display_usage() {

    echo "Usage: run_test.sh [-l] all | periodic | frequency [test]"
}

echo Copying files... 
cp -r tests/data/raw/* data/raw
cp tests/data/external/* data/external

echo Disabling downloading of dataset...
sed -e '27,39 s/^/#/' -e  's/rules.download_dataset.output/"data\/raw\/\{pid\}\/\{sensor\}_raw\.csv"/' rules/preprocessing.smk > tmp
mv tmp rules/preprocessing.smk

echo $1 
echo $2

if [ $# -eq 1 ]
then
    if [ $1 == '-l' ] || [ $1 == 'test' ]
    then
        display_usage
    elif [ $1 == 'all' ]
    then
        run_periodic_pipeline
        run_frequency_pipeline
    elif [ $1 == 'periodic' ]
    then
        run_periodic_pipeline
    elif [ $1 == 'frequency' ]
    then
        run_frequency_pipeline
    else
        display_usage
    fi
elif [ $# -gt 1 ]
then
    if [ $1 == '-l' ]
    then
        clean_old_data
        if [ $2 == 'all' ]
        then
            run_periodic_pipeline
            run_frequency_pipeline
            if [ $# -gt 2 ] && [ $3 == 'test' ]
            then
                run_periodic_test
                run_frequency_test
            else
                display_usage
            fi
        elif [ $2 == 'periodic' ]
        then
            run_periodic_pipeline
            if [ $# -gt 2 ] && [ $3 == 'test' ]
            then
                run_periodic_test
            else
                display_usage
            fi
        elif [ $2 == 'frequency' ]
        then
            run_frequency_pipeline
            if [ $# -gt 2 ] && [ $3 == 'test' ]
            then
                run_frequency_test
            else
                display_usage
            fi
        else
            display_usage
        fi
        mv bak rules/preprocessing.smk
    elif [ $1 == 'all' ]
    then
        run_periodic_pipeline
        run_frequency_pipeline
        if [ $2 == 'test' ]
        then
            run_periodic_test
            run_frequency_test
        else
            display_usage
        fi
    elif [ $1 == 'periodic' ]
    then
        run_periodic_pipeline
        if [ $2 == 'test' ]
        then
            run_periodic_test
        else
            display_usage
        fi
    elif [ $1 == 'frequency' ]
    then
        run_frequency_pipeline
        if [ $2 == 'test' ]
        then
            run_frequency_test
        else
            display_usage
        fi
    else
        display_usage
    fi
else
    display_usage
fi
