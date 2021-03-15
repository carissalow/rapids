# Behavioral Features Introduction

Every device sensor has a corresponding config section in `config.yaml`, these sections follow a similar structure and we'll use `PHONE_ACCELEROMETER` as an example to explain this structure.

!!! hint
    - We recommend reading this page if you are using RAPIDS for the first time
    - All computed sensor features are stored under `/data/processed/features` on files per sensor, per participant and per study (all participants).
    - Every time you change any sensor parameters, provider parameters or provider features, all the necessary files will be updated as soon as you execute RAPIDS.
    - In short, to extract features offered by a provider, you need to set its `[COMPUTE]` flag to `TRUE`, configure any of its parameters, and [execute](../../setup/execution) RAPIDS.


### Explaining the config.yaml sensor sections with an example

Each sensor section follows the same structure. Click on the numbered markers to know more.

``` { .yaml .annotate }
PHONE_ACCELEROMETER: # (1)

    CONTAINER: accelerometer # (2)

    PROVIDERS: # (3)
        RAPIDS:
            COMPUTE: False # (4)
            FEATURES: ["maxmagnitude", "minmagnitude", "avgmagnitude", "medianmagnitude", "stdmagnitude"]

            SRC_SCRIPT: src/features/phone_accelerometer/rapids/main.py
        
        PANDA:
            COMPUTE: False
            VALID_SENSED_MINUTES: False
            FEATURES: # (5)
                exertional_activity_episode: ["sumduration", "maxduration", "minduration", "avgduration", "medianduration", "stdduration"]
                nonexertional_activity_episode: ["sumduration", "maxduration", "minduration", "avgduration", "medianduration", "stdduration"]

                        # (6)
            SRC_SCRIPT: src/features/phone_accelerometer/panda/main.py
```

--8<--- "docs/snippets/feature_introduction_example.md"

These are descriptions of each marker for accessibility:

--8<--- "docs/snippets/feature_introduction_example.md"