Data Cleaning
=============

This module is to clean the extracted sensor features before merging it with the target labels.
    
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
|`[COMPUTE]` | Set to `True` to clean sensor features for individual participants from the `RAPIDS` provider|
|`[IMPUTE_SELECTED_EVENT_FEATURES]`     | Fill NA with 0 for the selected event features, see table below
|`[COLS_NAN_THRESHOLD]`                 | Discard columns with missing value ratios higher than `[COLS_NAN_THRESHOLD]`. Set to 1 to disable
|`[COLS_VAR_THRESHOLD]`                 | Set to `True` to discard columns with zero variance
|`[ROWS_NAN_THRESHOLD]`                 | Discard rows with missing value ratios higher than `[ROWS_NAN_THRESHOLD]`. Set to 1 to disable
|`[DATA_YIELDED_HOURS_RATIO_THRESHOLD]` | Discard rows with `phone_data_yield_rapids_ratiovalidyieldedhours` feature less than `[DATA_YIELDED_HOURS_RATIO_THRESHOLD]`. Set to 0 to disable
|`DROP_HIGHLY_CORRELATED_FEATURES`      | Discard highly correlated features, see table below

Parameters description for `[ALL_CLEANING_INDIVIDUAL][PROVIDERS][RAPIDS][IMPUTE_SELECTED_EVENT_FEATURES]`:

|Parameters                             | Description                                                    |
|-------------------------------------- |----------------------------------------------------------------|
|`[COMPUTE]`                            | Set to `True` to fill NA with 0 for the selected event features
|`[MIN_DATA_YIELDED_MINUTES_TO_IMPUTE]` | Assume the selected event sensor is working when phone_data_yield_rapids_ratiovalidyieldedminutes > `[MIN_DATA_YIELDED_MINUTES_TO_IMPUTE]`. |

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

??? info "2. Discard unreliable rows."
    Extracted features might be not reliable if the sensor only works for a short period during a time segment. In this step, we discard rows when the `phone_data_yield_rapids_ratiovalidyieldedhours` column is less than the `[DATA_YIELDED_HOURS_RATIO_THRESHOLD]` parameter. We do not recommend you to skip this step, but you can do it by setting `[DATA_YIELDED_HOURS_RATIO_THRESHOLD]` to 0.

??? info "3. Discard columns (features) with too many missing values."
    In this step, we discard columns with missing value ratios higher than `[COLS_NAN_THRESHOLD]`. We do not recommend you to skip this step, but you can do it by setting `[COLS_NAN_THRESHOLD]` to 1.

??? info "4. Discard columns (features) with zero variance."
    In this step, we discard columns with zero variance. We do not recommend you to skip this step, but you can do it by setting `[COLS_VAR_THRESHOLD]` to `False`.

??? info "5. Drop highly correlated features."
    As highly correlated features might not bring additional information and will increase the complexity of our model, we drop them in this step. The absolute values of pair-wise correlations are calculated. It is regarded as valid only if the ratio of this pair of columns (features) are less than `[DROP_HIGHLY_CORRELATED_FEATURES][MIN_OVERLAP_FOR_CORR_THRESHOLD]`. If two variables have a valid correlation higher than `[DROP_HIGHLY_CORRELATED_FEATURES][CORR_THRESHOLD]`, we looks at the mean absolute correlation of each variable and removes the variable with the largest mean absolute correlation. This step can be skip by setting `[DROP_HIGHLY_CORRELATED_FEATURES][COMPUTE]` to `False`.

??? info "6. Discard rows with too many missing values."
    In this step, we discard rows with missing value ratios higher than `[ROWS_NAN_THRESHOLD]`. We do not recommend you to skip this step, but you can do it by setting `[ROWS_NAN_THRESHOLD]` to 1.




## Clean sensor features for all participants.

!!! info "File Sequence"
    ```bash
    - data/processed/features/all_participants/all_sensor_features.csv
    - data/processed/features/all_participants/all_sensor_features_cleaned_{provider_key}.csv
    ```


### RAPIDS provider

Parameters description and the steps are the same as the above [RAPIDS provider](#rapids-provider) section for individual participants.


