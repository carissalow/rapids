#!/bin/bash

run_pipeline() {
    if [ "$TYPE" == 'frequency' ]
    then
        CONFIG_FILE="./tests/settings/frequency_config.yaml"
    elif [ "$TYPE" == 'event' ]
    then 
        CONFIG_FILE="./tests/settings/event_config.yaml"
    else
        CONFIG_FILE="./tests/settings/periodic_config.yaml"
    fi

    echo "Copying participant files"
    mkdir -p data/external/participant_files/
    cp -r tests/data/external/participant_files/* data/external/participant_files/

    echo $TYPE
    echo "Deleting old outputs"
    snakemake --configfile=$(echo $CONFIG_FILE) --delete-all-output -j1 

    echo "Running RAPIDS"
    snakemake --configfile=$(echo $CONFIG_FILE) -R pull_phone_data -j1
}

display_usage() {
    echo "Usage: run_test.sh [-t|--type] [periodic | frequency | event] [-a|--action] [ all | run | both]"
    exit 1
}

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -t|--type)
    TYPE="$2"
    shift # past argument
    shift # past value
    ;;
    -a|--action)
    ACTION="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done

if { [ "$TYPE" == 'all' ]; }
then
    TYPE="frequency"
    run_pipeline
    python tests/scripts/run_tests.py frequency
    TYPE="periodic"
    run_pipeline
    python tests/scripts/run_tests.py periodic
    TYPE="event"
    run_pipeline
    python tests/scripts/run_tests.py event
else
    if { [ "$ACTION" != 'test' ] && [ "$ACTION" != 'run' ] && [ "$ACTION" != 'both' ]; }
    then
        display_usage
    fi

    if { [ "$TYPE" != 'frequency' ] && [ "$TYPE" != 'periodic' ] && [ "$TYPE" != 'event' ]; }
    then
        display_usage
    fi

    if { [ "$ACTION" == 'run' ] || [ "$ACTION" == 'both' ]; }
    then
        run_pipeline
    fi
    if { [ "$ACTION" == 'test' ] || [ "$ACTION" == 'both' ]; }
    then
        python tests/scripts/run_tests.py $(echo $TYPE)
    fi
fi
