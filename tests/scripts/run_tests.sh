#!/bin/bash

run_pipeline() {
    if [ "$TYPE" == 'stz_frequency' ]
    then
        CONFIG_FILE="./tests/settings/stz_frequency_config.yaml"
    elif [ "$TYPE" == 'mtz_frequency' ]
    then
        CONFIG_FILE="./tests/settings/mtz_frequency_config.yaml"
    elif [ "$TYPE" == 'stz_event' ]
    then 
        CONFIG_FILE="./tests/settings/stz_event_config.yaml"
    elif [ "$TYPE" == 'mtz_event' ]
    then
        CONFIG_FILE="./tests/settings/mtz_event_config.yaml"
    elif [ "$TYPE" == 'stz_periodic' ]
    then
        CONFIG_FILE="./tests/settings/stz_periodic_config.yaml"
    else
        CONFIG_FILE="./tests/settings/mtz_periodic_config.yaml"
    fi

    echo "Copying participant files"
    mkdir -p data/external/participant_files/
    cp -r tests/data/external/participant_files/* data/external/participant_files/

    echo $TYPE
    echo "Deleting old outputs"
    snakemake --configfile=$(echo $CONFIG_FILE) --delete-all-output -j1 || exit

    echo "Running RAPIDS"
    snakemake --configfile=$(echo $CONFIG_FILE) -R pull_phone_data -j2 || exit
}

display_usage() {
    echo "Usage: run_test.sh [-t|--type] [all | stz_periodic | mtz_periodic | stz_frequency | mtz_frequency | stz_event | mtz_event] [-a|--action] [test | run | both]"
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
    TYPE="stz_frequency"
    run_pipeline
    python tests/scripts/run_tests.py stz_frequency || exit
    TYPE="mtz_frequency"
    run_pipeline
    python tests/scripts/run_tests.py mtz_frequency || exit
    TYPE="stz_periodic"
    run_pipeline
    python tests/scripts/run_tests.py stz_periodic || exit
    TYPE="mtz_periodic"
    run_pipeline
    python tests/scripts/run_tests.py mtz_periodic || exit
    TYPE="stz_event"
    run_pipeline
    python tests/scripts/run_tests.py stz_event || exit
    TYPE="mtz_event"
    run_pipeline
    python tests/scripts/run_tests.py mtz_event || exit
else
    if { [ "$ACTION" != 'test' ] && [ "$ACTION" != 'run' ] && [ "$ACTION" != 'both' ]; }
    then
        display_usage
    fi

    if { [ "$TYPE" != 'stz_frequency' ] && [ "$TYPE" != 'mtz_frequency' ] && [ "$TYPE" != 'stz_periodic' ] && [ "$TYPE" != 'mtz_periodic' ] && [ "$TYPE" != 'stz_event' ] && [ "$TYPE" != 'mtz_event' ]; }
    then
        display_usage
    fi

    if { [ "$ACTION" == 'run' ] || [ "$ACTION" == 'both' ]; }
    then
        run_pipeline
    fi
    if { [ "$ACTION" == 'test' ] || [ "$ACTION" == 'both' ]; }
    then
        python tests/scripts/run_tests.py $(echo $TYPE) || exit
    fi
fi
