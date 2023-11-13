# Fitbit Steps Intraday

Sensor parameters description for `[FITBIT_STEPS_INTRADAY]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Container where your steps intraday data is stored, depending on the data stream you are using this can be a database table, a CSV file, etc. |
|`[EXCLUDE_SLEEP]` | Step data will be excluded if it was logged during sleep periods when at least one `[EXCLUDE]` flag is set to `True`. Sleep can be delimited by (1) a fixed period that repeats on every day if `[TIME_BASED][EXCLUDE]` is True or (2) by Fitbit summary sleep episodes if `[FITBIT_BASED][EXCLUDE]` is True. If both are True (3), we use all Fitbit sleep episodes as well as the time-based episodes that do not overlap with any Fitbit episodes. If `[TIME_BASED][EXCLUDE]` is True, make sure Fitbit sleep summary container points to a valid table or file.

## RAPIDS provider

!!! info "Available time segments"
    - Available for all time segments

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/fitbit_steps_intraday_raw.csv
    - data/raw/{pid}/fitbit_steps_intraday_with_datetime.csv
    - data/raw/{pid}/fitbit_sleep_summary_raw.csv (Only when [EXCLUDE_SLEEP][EXCLUDE]=True and [EXCLUDE_SLEEP][TYPE]=FITBIT_BASED)
    - data/interim/{pid}/fitbit_steps_intraday_with_datetime_exclude_sleep.csv (Only when [EXCLUDE_SLEEP][EXCLUDE]=True)
    - data/interim/{pid}/fitbit_steps_intraday_features/fitbit_steps_intraday_{language}_{provider_key}.csv
    - data/processed/features/{pid}/fitbit_steps_intraday.csv
    ```


Parameters description for `[FITBIT_STEPS_INTRADAY][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`                            | Set to `True` to extract `FITBIT_STEPS_INTRADAY` features from the `RAPIDS` provider|
|`[FEATURES]`                           |         Features to be computed from steps intraday data, see table below           |
|`[REFERENCE_HOUR]`                     | The reference point from which `firststeptime` or `laststeptime` is to be computed, default is midnight |
|`[THRESHOLD_ACTIVE_BOUT]`              | Every minute with Fitbit steps data wil be labelled as `sedentary` if its step count is below this threshold, otherwise, `active`.    |  
|`[THRESHOLD_DEVICE_NONWEAR_TIME]`      | Number of consecutive minutes with zero step counts at or above which to label as a period of device non-wear time. Used only to calculate `uncensoredmeancadence` feature. Default is 1, which corresponds to a measure of mean non-zero walking cadence. Set to 60 for consistency with Tudor-Locke et al.'s definition. Set to 0 to disable.  |   
|`[THRESHOLD_MINUTE_LEVEL_STEP_COUNT]`  | Rows with minute-level step counts above this threshold will be excluded from feature computation. Default is 1000; set to 0 to disable. |  
|`[INCLUDE_ZERO_STEP_ROWS]`             | Whether or not to include time segments with a 0 step count during the whole day.                          |


Features description for `[FITBIT_STEPS_INTRADAY][PROVIDERS][RAPIDS]`:

|Feature                                    |Units          |Description                                                  |
|------------------------------------------ |-------------- |-------------------------------------------------------------|
|sumsteps                                   |steps          |The total step count during a time segment.
|maxsteps                                   |steps          |The maximum step count during a time segment.
|minsteps                                   |steps          |The minimum step count during a time segment.
|avgsteps                                   |steps          |The average step count during a time segment.
|stdsteps                                   |steps          |The standard deviation of step count during a time segment.
|firststeptime                              |minutes        |Minutes until the first non-zero step count.
|laststeptime                               |minutes        |Minutes until the last non-zero step count.
|countepisodesedentarybout                  |bouts          |Number of sedentary bouts during a time segment.
|sumdurationsedentarybout                   |minutes        |Total duration of all sedentary bouts during a time segment.
|maxdurationsedentarybout                   |minutes        |The maximum duration of any sedentary bout during a time segment.
|mindurationsedentarybout                   |minutes        |The minimum duration of any sedentary bout during a time segment.
|avgdurationsedentarybout                   |minutes        |The average duration of sedentary bouts during a time segment.
|stddurationsedentarybout                   |minutes        |The standard deviation of the duration of sedentary bouts during a time segment.
|countepisodeactivebout                     |bouts          |Number of active bouts during a time segment.
|sumdurationactivebout                      |minutes        |Total duration of all active bouts during a time segment.
|maxdurationactivebout                      |minutes        |The maximum duration of any active bout during a time segment.
|mindurationactivebout                      |minutes        |The minimum duration of any active bout during a time segment.
|avgdurationactivebout                      |minutes        |The average duration of active bouts during a time segment.
|stddurationactivebout                      |minutes        |The standard deviation of the duration of active bouts during a time segment.  
|activetosedentarytransitionprobability     |unitless       |Active-to-sedentary transition probability (ASTP). The reciprocal of the mean active bout length during a time segment. Bounded by 0 and 1.    
|sumdurationactivitylessthan5minutes        |minutes        |Sum of duration of active bouts less than 5 minutes in length during a time segment.   
|sumdurationactivity5to105minutes           |minutes        |Sum of duration of active bouts between 5 and 10 minutes in length, inclusive, during a time segment.  
|sumdurationactivitygreaterthan10minutes    |minutes        |Sum of duration of active bouts more than 10 minutes in length during a time segment.  
|ginicoefficient                            |unitless       |Measure of (absolute, not squared) variability of active bout durations normalized by the average active bout duration during a time segment. When the Gini coefficient is close to 1, it indicates that total time is accumulated via a small number of longer bouts. Conversely, when the Gini coefficient is close to 0, it indicates that all bouts contribute equally to total time.  
|meancadence                                |steps/minute   |The average steps/minute during a time segment.  
|uncensoredmeancadence                      |steps/minute   |The total raw steps accumulated divided by device wear time during a time segment.  
|peak1minutecadence                         |steps/minute   |The steps/minute recorded for the highest single minute during a time segment; this is also equivalent to the maximum 1 minute cadence.  
|peak30minutecadence                        |steps/minute   |The average steps/minute recorded for the 30 highest, but not necessarily consecutive, minutes during a time segment.  
|peak60minutecadence                        |steps/minute   |The average steps/minute recorded for the 60 highest, but not necessarily consecutive, minutes during a time segment.  
|max5minutecadence                          |steps/minute   |The average steps/minute of the maximum number of steps obtained over 5 continuous minutes during a time segment.  
|max20minutecadence                         |steps/minute   |The average steps/minute of the maximum number of steps obtained over 20 continuous minutes during a time segment.  
|max30minutecadence                         |steps/minute   |The average steps/minute of the maximum number of steps obtained over 30 continuous minutes during a time segment.  
|max60minutecadence                         |steps/minute   |The average steps/minute of the maximum number of steps obtained over 60 continuous minutes during a time segment.  
|totalminutes0cadence                       |minutes        |Total minutes with no movement (0 steps/minute) during a time segment; an indicator of non-ambulatory or sedentary time.  
|totalminutes1to19cadence                   |minutes        |Total minutes with incidental movement ([1, 20) steps/minute) during a time segment.    
|totalminutes20to39cadence                  |minutes        |Total minutes with sporadic movement ([20, 40) steps/minute) during a time segment.  
|totalminutes40to59cadence                  |minutes        |Total minutes with purposeful steps ([40, 60) steps/minute) during a time segment.  
|totalminutes60to79cadence                  |minutes        |Total minutes with slow walking ([60, 80) steps/minute) during a time segment.  
|totalminutes80to99cadence                  |minutes        |Total minutes with medium walking ([80, 100) steps/minute) during a time segment.  
|totalminutes100to119cadence                |minutes        |Total minutes with brisk walking ([100, 120) steps/minute) during a time segment.  
|totalminutes120pluscadence                 |minutes        |Total minutes with all faster ambulation ([120, Inf) steps/minute) during a time segment.  
|totalminutesabove0cadence                  |minutes        |Total minutes with any movement (>0 steps/minute) during a time segment.     
|totalminutesabove19cadence                 |minutes        |Total minutes with non-incidental movement (>19 steps/minute) during a time segment.  
|totalminutesabove100cadence                |minutes        |Total minutes with movement of at least moderate intensity (>100 steps/minute) during a time segment; this threshold was based on findings from ostensibly health individuals without movement impairments.    

!!! note "Assumptions/Observations"
    
    1. _Active and sedentary bouts_. If the step count per minute is smaller than `THRESHOLD_ACTIVE_BOUT` (default value is 10), that minute is labelled as sedentary, otherwise, is labelled as active. Active and sedentary bouts are periods of consecutive minutes labelled as `active` or `sedentary`.  
    2. Activity fragmentation features are based on the following papers: [Di et al., 2017](https://doi.org/10.1101/182337), [Schrack et al., 2019](https://doi.org/10.1093/gerona/gly243), and [Wanigatunga et al., 2019](https://doi.org/10.1001/jamanetworkopen.2019.12352).   
    3. Walking cadence features are based on the following review paper: [Tudor-Locke et al., 2018](https://doi.org/10.1136/bjsports-2017-097628).    

