# Behavioral Features Introduction

Every device sensor has a corresponding config section in `config.yaml`, these sections follow a similar structure and we'll use `PHONE_ACCELEROMETER` as an example to explain this structure.

!!! hint
    - We recommend reading this page if you are using RAPIDS for the first time
    - All computed sensor features are stored under `/data/processed/features` on files per sensor, per participant and per study (all participants).
    - Every time you change any sensor parameters, provider parameters or provider features, all the necessary files will be updated as soon as you execute RAPIDS.


!!! example "Config section example for `PHONE_ACCELEROMETER`"

    ```yaml
    # 1) Config section
    PHONE_ACCELEROMETER:
        # 2) Parameters for PHONE_ACCELEROMETER
        TABLE: accelerometer

        # 3) Providers for PHONE_ACCELEROMETER
        PROVIDERS:
            # 4) RAPIDS provider
            RAPIDS:
                # 4.1) Parameters of RAPIDS provider of PHONE_ACCELEROMETER
                COMPUTE: False
                # 4.2) Features of RAPIDS provider of PHONE_ACCELEROMETER
                FEATURES: ["maxmagnitude", "minmagnitude", "avgmagnitude", "medianmagnitude", "stdmagnitude"]
                SRC_FOLDER: "rapids" # inside src/features/phone_accelerometer
                SRC_LANGUAGE: "python"
            
            # 5) PANDA provider
            PANDA:
                # 5.1) Parameters of PANDA provider of PHONE_ACCELEROMETER
                COMPUTE: False
                VALID_SENSED_MINUTES: False
                # 5.2) Features of PANDA provider of PHONE_ACCELEROMETER
                FEATURES:
                    exertional_activity_episode: ["sumduration", "maxduration", "minduration", "avgduration", "medianduration", "stdduration"]
                    nonexertional_activity_episode: ["sumduration", "maxduration", "minduration", "avgduration", "medianduration", "stdduration"]
                SRC_FOLDER: "panda" # inside src/features/phone_accelerometer
                SRC_LANGUAGE: "python"
    ```

## Sensor Parameters
Each sensor configuration section has a "parameters" subsection (see `#2` in the example). These are parameters that affect different aspects of how the raw data is downloaded, and processed. The `TABLE` parameter exists for every sensor, but some sensors will have extra parameters like [`[PHONE_LOCATIONS]`](../phone-locations/). We explain these parameters in a table at the top of each sensor documentation page.

## Sensor Providers
Each sensor configuration section can have zero, one or more behavioral feature **providers** (see `#3` in the example). A provider is a script created by the core RAPIDS team or other researchers that extracts behavioral features for that sensor. In this example, accelerometer has two providers: RAPIDS (see `#4`) and PANDA (see `#5`).

### Provider Parameters
Each provider has parameters that affect the computation of the behavioral features it offers (see `#4.1` or `#5.1` in the example). These parameters will include at least a `[COMPUTE]` flag that you switch to `True` to extract a provider's behavioral features. 

We explain every provider's parameter in a table under the `Parameters description` heading on each provider documentation page.

### Provider Features
Each provider offers a set of behavioral features (see `#4.2` or `#5.2` in the example). For some providers these features are grouped in an array (like those for `RAPIDS` provider in `#4.2`) but for others they are grouped in a collection of arrays depending on the meaning and purpose of those features (like those for `PANDAS` provider in `#5.2`). In either case, you can delete the features you are not interested in and they will not be included in the sensor's output feature file. 

We explain each behavioral feature in a table under the `Features description` heading on each provider documentation page.
