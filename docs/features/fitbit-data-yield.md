# Fitbit Data Yield

We use Fitbit **heart rate intraday** data to extract data yield features. Fitbit data yield features can be used to remove rows ([time segments](../../setup/configuration/#time-segments)) that do not contain enough Fitbit data. You should decide what is your "enough" threshold depending on the time a participant was supposed to be wearing their Fitbit, the length of your study, and the rates of missing data that your analysis could handle.

!!! hint "Why is Fitbit data yield important?"
    Imagine that you want to extract `FITBIT_STEPS_SUMMARY` features on daily segments (`00:00` to `23:59`). Let's say that on day 1 the Fitbit logged 6k as the total step count and the heart rate sensor logged 24 hours of data and on day 2 the Fitbit logged 101 as the total step count and the heart rate sensor logged 2 hours of data. It’s very likely that on day 2 you walked during the other 22 hours so including this day in your analysis could bias your results.
Sensor parameters description for `[FITBIT_DATA_YIELD]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;          | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[SENSORS]`| The Fitbit sensor we considered for calculating the Fitbit data yield features. We only support `FITBIT_HEARTRATE_INTRADAY` since sleep data is commonly collected only overnight, and step counts are 0 even when not wearing the Fitbit device.

## RAPIDS provider

Before explaining the data yield features, let's define the following relevant concepts:

- A valid minute is any 60 second window when Fitbit heart rate intraday sensor logged at least 1 row of data
- A valid hour is any 60 minute window with at least X valid minutes. The X or threshold is given by `[MINUTE_RATIO_THRESHOLD_FOR_VALID_YIELDED_HOURS]`

!!! info "Available time segments and platforms"
    - Available for all time segments

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/fitbit_heartrate_intraday_raw.csv
    - data/raw/{pid}/fitbit_heartrate_intraday_with_datetime.csv
    - data/interim/{pid}/fitbit_data_yield_features/fitbit_data_yield_{language}_{provider_key}.csv
    - data/processed/features/{pid}/fitbit_data_yield.csv
    ```


Parameters description for `[FITBIT_DATA_YIELD][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`| Set to `True` to extract `FITBIT_DATA_YIELD` features from the `RAPIDS` provider|
|`[FEATURES]` |  Features to be computed, see table below
|`[MINUTE_RATIO_THRESHOLD_FOR_VALID_YIELDED_HOURS]` | The proportion `[0.0 ,1.0]` of valid minutes in a 60-minute window necessary to flag that window as valid.


Features description for `[FITBIT_DATA_YIELD][PROVIDERS][RAPIDS]`:

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
|ratiovalidyieldedminutes   |-          | The ratio between the number of valid minutes and the duration in minutes of a time segment.
|ratiovalidyieldedhours     |-          | The ratio between the number of valid hours and the duration in hours of a time segment. If the time segment is shorter than 1 hour this feature will always be 1.


!!! note "Assumptions/Observations"
    
    1. We recommend using `ratiovalidyieldedminutes` on time segments that are shorter than two or three hours and `ratiovalidyieldedhours` for longer segments. This is because relying on yielded minutes only can be misleading when a big chunk of those missing minutes are clustered together. 
    
        For example, let's assume we are working with a 24-hour time segment that is missing 12 hours of data. Two extreme cases can occur: 

        <ol type="A">
        <li>the 12 missing hours are from the beginning of the segment or </li>
        <li>30 minutes could be missing from every hour (24 * 30 minutes = 12 hours).</li>
        </ol>
        
        `ratiovalidyieldedminutes` would be 0.5 for both `a` and `b` (hinting the missing circumstances are similar). However, `ratiovalidyieldedhours` would be 0.5 for `a` and 1.0 for `b` if `[MINUTE_RATIO_THRESHOLD_FOR_VALID_YIELDED_HOURS]` is between [0.0 and 0.49] (hinting that the missing circumstances might be more favorable for `b`. In other words, sensed data for `b` is more evenly spread compared to `a`.
    
    2. We assume your Fitbit intraday data was sampled (requested form the Fitbit API) at 1 minute intervals, if the interval is longer, for example 15 minutes, you need to take into account that valid minutes and valid hours ratios are going to be small (for example you would have at most 4 “minutes” of data per hour because you would have four 15-minute windows) and so you should adjust your thresholds to include and exclude rows accordingly. If you are in this situation, get in touch with us, we could implement this use case but we are not sure there is enough demand for it at the moment since you can control the sampling rate of the data you request from Fitbit API.
