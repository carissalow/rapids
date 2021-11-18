Data Cleaning
=============

The goal of this module is to perform basic clean tasks on the behavioral features that RAPIDS computes. You might need to do further processing depending on your analysis objectives. This module can clean features at the individual level and at the study level. If you are interested in creating individual models (using each participant's features independently of the others) use [`ALL_CLEANING_INDIVIDUAL`]. If you are interested in creating population models (using everyone's data in the same model) use [`ALL_CLEANING_OVERALL`]
    
## Clean sensor features for individual participants

!!! info "File Sequence"
    ```bash
    - data/processed/features/{pid}/all_sensor_features.csv
    - data/processed/features/{pid}/all_sensor_features_cleaned_{provider_key}.csv
    ```

### RAPIDS provider

Parameters description for `[ALL_CLEANING_INDIVIDUAL][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]` | Set to `True` to execute the cleaning tasks described below. You can use the parameters of each task to tweak them or deactivate them|
|`[IMPUTE_SELECTED_EVENT_FEATURES]`     | Fill NAs with 0 only for event-based features, see table below
|`[COLS_NAN_THRESHOLD]`                 | Discard columns with missing value ratios higher than `[COLS_NAN_THRESHOLD]`. Set to 1 to disable
|`[COLS_VAR_THRESHOLD]`                 | Set to `True` to discard columns with zero variance
|`[ROWS_NAN_THRESHOLD]`                 | Discard rows with missing value ratios higher than `[ROWS_NAN_THRESHOLD]`. Set to 1 to disable
|`[DATA_YIELD_FEATURE]`                 | `RATIO_VALID_YIELDED_HOURS` or `RATIO_VALID_YIELDED_MINUTES`
|`[DATA_YIELD_RATIO_THRESHOLD]`         | Discard rows with `ratiovalidyieldedhours` or `ratiovalidyieldedminutes` feature less than `[DATA_YIELD_RATIO_THRESHOLD]`. The feature name is determined by `[DATA_YIELD_FEATURE]` parameter. Set to 0 to disable
|`DROP_HIGHLY_CORRELATED_FEATURES`      | Discard highly correlated features, see table below

Parameters description for `[ALL_CLEANING_INDIVIDUAL][PROVIDERS][RAPIDS][IMPUTE_SELECTED_EVENT_FEATURES]`:

|Parameters                             | Description                                                    |
|-------------------------------------- |----------------------------------------------------------------|
|`[COMPUTE]`                            | Set to `True` to fill NAs with 0 for phone event-based features
|`[MIN_DATA_YIELDED_MINUTES_TO_IMPUTE]` | Any feature value in a time segment instance with phone data yield > `[MIN_DATA_YIELDED_MINUTES_TO_IMPUTE]` will be replaced with a zero. See below for an explanation. |

Parameters description for `[ALL_CLEANING_INDIVIDUAL][PROVIDERS][RAPIDS][DROP_HIGHLY_CORRELATED_FEATURES]`:

|Parameters                             | Description                                                    |
|-------------------------------------- |----------------------------------------------------------------|
|`[COMPUTE]`                            | Set to `True` to drop highly correlated features
|`[MIN_OVERLAP_FOR_CORR_THRESHOLD]`     | Minimum ratio of observations required per pair of columns (features) to be considered as a valid correlation. 
|`[CORR_THRESHOLD]` | The absolute values of pair-wise correlations are calculated. If two variables have a valid correlation higher than `[CORR_THRESHOLD]`, we looks at the mean absolute correlation of each variable and removes the variable with the largest mean absolute correlation.

Steps to clean sensor features for individual participants. It only considers the **phone sensors** currently.

??? info "1. Fill NA with 0 for the selected event features."
    Some event features should be zero instead of NA. In this step, we fill those missing features with 0 when the `phone_data_yield_rapids_ratiovalidyieldedminutes` column is higher than the `[IMPUTE_SELECTED_EVENT_FEATURES][MIN_DATA_YIELDED_MINUTES_TO_IMPUTE]` parameter. Plugins such as Activity Recognition sensor are not considered. You can skip this step by setting `[IMPUTE_SELECTED_EVENT_FEATURES][COMPUTE]` to `False`.
    
    Take phone calls sensor as an example. If there are no calls records during a time segment for a participant, then (1) the calls sensor was not working during that time segment; or (2) the calls sensor was working and the participant did not have any calls during that time segment. To differentiate these two situations, we assume the selected sensors are working when `phone_data_yield_rapids_ratiovalidyieldedminutes > [MIN_DATA_YIELDED_MINUTES_TO_IMPUTE]`.

    The following phone event-based features are considered currently:

      - Application foreground: countevent, countepisode, minduration, maxduration, meanduration, sumduration.
      - Battery: all features.
      - Calls: count, distinctcontacts, sumduration, minduration, maxduration, meanduration, modeduration.
      - Keyboard: sessioncount, averagesessionlength, changeintextlengthlessthanminusone, changeintextlengthequaltominusone, changeintextlengthequaltoone, changeintextlengthmorethanone, maxtextlength, totalkeyboardtouches.
      - Messages: count, distinctcontacts.
      - Screen: sumduration, maxduration, minduration, avgduration, countepisode.
      - WiFi: all connected and visible features.

??? info "2. Discard unreliable rows."
    Extracted features might be not reliable if the sensor only works for a short period during a time segment. In this step, we discard rows when the `phone_data_yield_rapids_ratiovalidyieldedminutes` column or the `phone_data_yield_rapids_ratiovalidyieldedhours` column is less than the `[DATA_YIELD_RATIO_THRESHOLD]` parameter. We recommend using `phone_data_yield_rapids_ratiovalidyieldedminutes` column (set `[DATA_YIELD_FEATURE]` to `RATIO_VALID_YIELDED_MINUTES`) on time segments that are shorter than two or three hours and `phone_data_yield_rapids_ratiovalidyieldedhours` (set `[DATA_YIELD_FEATURE]` to `RATIO_VALID_YIELDED_HOURS`) for longer segments. We do not recommend you to skip this step, but you can do it by setting `[DATA_YIELD_RATIO_THRESHOLD]` to 0.

??? info "3. Discard columns (features) with too many missing values."
    In this step, we discard columns with missing value ratios higher than `[COLS_NAN_THRESHOLD]`. We do not recommend you to skip this step, but you can do it by setting `[COLS_NAN_THRESHOLD]` to 1.

??? info "4. Discard columns (features) with zero variance."
    In this step, we discard columns with zero variance. We do not recommend you to skip this step, but you can do it by setting `[COLS_VAR_THRESHOLD]` to `False`.

??? info "5. Drop highly correlated features."
    As highly correlated features might not bring additional information and will increase the complexity of a model, we drop them in this step. The absolute values of pair-wise correlations are calculated. Each correlation vector between two variables is regarded as valid only if the ratio of valid value pairs (i.e. non NA pairs) is greater than or equal to `[DROP_HIGHLY_CORRELATED_FEATURES][MIN_OVERLAP_FOR_CORR_THRESHOLD]`. If two variables have a correlation coefficient higher than `[DROP_HIGHLY_CORRELATED_FEATURES][CORR_THRESHOLD]`, we look at the mean absolute correlation of each variable and remove the variable with the largest mean absolute correlation. This step can be skipped by setting `[DROP_HIGHLY_CORRELATED_FEATURES][COMPUTE]` to False.

??? info "6. Discard rows with too many missing values."
    In this step, we discard rows with missing value ratios higher than `[ROWS_NAN_THRESHOLD]`. We do not recommend you to skip this step, but you can do it by setting `[ROWS_NAN_THRESHOLD]` to 1. In other words, we are discarding time segments (e.g. days) that did not have enough data to be considered reliable. This step is similar to step 2 except the ratio is computed based on NA values instead of a phone data yield threshold.




## Clean sensor features for all participants

!!! info "File Sequence"
    ```bash
    - data/processed/features/all_participants/all_sensor_features.csv
    - data/processed/features/all_participants/all_sensor_features_cleaned_{provider_key}.csv
    ```


### RAPIDS provider

Parameters description and the steps are the same as the above [RAPIDS provider](#rapids-provider) section for individual participants.


