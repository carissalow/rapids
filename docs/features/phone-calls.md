# Phone Calls

Sensor parameters description for `[PHONE_CALLS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Data stream [container](../../datastreams/data-streams-introduction/) (database table, CSV file, etc.) where the calls data is stored

## RAPIDS Provider

!!! info "Available time segments and platforms"
    - Available for all time segments
    - Available for Android and iOS

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/phone_calls_raw.csv
    - data/raw/{pid}/phone_calls_with_datetime.csv
    - data/interim/{pid}/phone_calls_features/phone_calls_{language}_{provider_key}.csv
    - data/processed/features/{pid}/phone_calls.csv
    ```


Parameters description for `[PHONE_CALLS][PROVIDERS][RAPIDS]`:

| Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;        | Description |
|-------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|`[COMPUTE]`| Set to `True` to extract `PHONE_CALLS` features from the `RAPIDS` provider|
|`[FEATURES_TYPE]`| Set to `EPISODES` to extract features based on call episodes or `EVENTS` to extract features based on events.|
| `[CALL_TYPES]`   | The particular call_type that will be analyzed. The options for this parameter are incoming, outgoing or missed.                                                                                                                                                 |
| `[FEATURES]`    | Features to be computed for `outgoing`, `incoming`, and `missed` calls. Note that the same features are available for both incoming and outgoing calls, while missed calls has its own set of features. See the tables below. |


Features description for `[PHONE_CALLS][PROVIDERS][RAPIDS]` incoming and outgoing calls:

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
|count                    |calls      |Number of calls of a particular `call_type` occurred during a particular `time_segment`.
|distinctcontacts         |contacts   |Number of distinct contacts that are associated with a particular `call_type` for a particular `time_segment`
|meanduration             |seconds    |The mean duration of all calls of a particular `call_type` during a particular `time_segment`.
|sumduration              |seconds    |The sum of the duration of all calls of a particular `call_type` during a particular `time_segment`.
|minduration              |seconds    |The duration of the shortest call of a particular `call_type` during a particular `time_segment`.
|maxduration              |seconds    |The duration of the longest call of a particular `call_type` during a particular `time_segment`.
|stdduration              |seconds    |The standard deviation of the duration of all the calls of a particular `call_type` during a particular `time_segment`.
|modeduration             |seconds    |The mode of the duration of all the calls of a particular `call_type` during a particular `time_segment`.
|entropyduration          |nats       |The estimate of the Shannon entropy for the the duration of all the calls of a particular `call_type` during a particular `time_segment`.
|timefirstcall            |minutes    |The time in minutes between 12:00am (midnight) and the first call of `call_type`.
|timelastcall             |minutes    |The time in minutes between 12:00am (midnight) and the last call of `call_type`.
|countmostfrequentcontact |calls      |The number of calls of a particular `call_type` during a particular `time_segment` of the most frequent contact throughout the monitored period.

Features description for `[PHONE_CALLS][PROVIDERS][RAPIDS]` missed calls:

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
|count                      |calls      |Number of `missed` calls that occurred during a particular `time_segment`.
|distinctcontacts           |contacts   |Number of distinct contacts that are associated with `missed` calls for a particular `time_segment`
|timefirstcall              |minutes    |The time in hours from 12:00am (Midnight) that the first `missed` call occurred.
|timelastcall               |minutes    |The time in hours from 12:00am (Midnight) that the last `missed` call occurred.
|countmostfrequentcontact   |calls      |The number of `missed` calls during a particular `time_segment` of the most frequent contact throughout the monitored period.

!!! note "Assumptions/Observations"
    1. Traces for iOS calls are unique even for the same contact calling a participant more than once which renders `countmostfrequentcontact` meaningless and `distinctcontacts` equal to the total number of traces. 
    2. `[CALL_TYPES]` and `[FEATURES]` keys in `config.yaml` need to match. For example, `[CALL_TYPES]` `outgoing` matches the `[FEATURES]` key `outgoing`
    3. iOS calls data is transformed to match Android calls data format.
