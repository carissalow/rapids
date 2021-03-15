# Phone Messages

Sensor parameters description for `[PHONE_MESSAGES]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Data stream [container](../../datastreams/data-streams-introduction/) (database table, CSV file, etc.) where the messages data is stored

## RAPIDS provider

!!! info "Available time segments and platforms"
    - Available for all time segments
    - Available for Android only

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/phone_messages_raw.csv
    - data/raw/{pid}/phone_messages_with_datetime.csv
    - data/interim/{pid}/phone_messages_features/phone_messages_{language}_{provider_key}.csv
    - data/processed/features/{pid}/phone_messages.csv
    ```


Parameters description for `[PHONE_MESSAGES][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`| Set to `True` to extract `PHONE_MESSAGES` features from the `RAPIDS` provider|
|`[MESSAGES_TYPES]` |  The `messages_type` that will be analyzed. The options for this parameter are `received` or `sent`.
|`[FEATURES]` |         Features to be computed, see table below for `[MESSAGES_TYPES]` `received` and `sent`


Features description for `[PHONE_MESSAGES][PROVIDERS][RAPIDS]`:

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
|count                      |messages   |Number of messages of type `messages_type` that occurred during a particular `time_segment`.
|distinctcontacts           |contacts   |Number of distinct contacts that are associated with a particular `messages_type` during a particular `time_segment`.
|timefirstmessages          |minutes    |Number of minutes between 12:00am (midnight) and the first `message` of a particular `messages_type` during a particular `time_segment`.
|timelastmessages           |minutes    |Number of minutes between 12:00am (midnight) and the last `message` of a particular `messages_type` during a particular `time_segment`.
|countmostfrequentcontact   |messages   |Number of messages from the contact with the most messages of `messages_type` during a `time_segment` throughout the whole dataset of each participant.

!!! note "Assumptions/Observations"
    1. `[MESSAGES_TYPES]` and `[FEATURES]` keys in `config.yaml` need to match. For example, `[MESSAGES_TYPES]` `sent` matches the `[FEATURES]` key `sent`


