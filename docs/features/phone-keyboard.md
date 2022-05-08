# Phone Keyboard

Sensor parameters description for `[PHONE_KEYBOARD]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Data stream [container](../../datastreams/data-streams-introduction/) (database table, CSV file, etc.) where the keyboard data is stored

## RAPIDS provider

!!! info "Available time segments and platforms"
    - Available for all time segments
    - Available for Android only

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/phone_keyboard_raw.csv
    - data/raw/{pid}/phone_keyboard_with_datetime.csv
    - data/interim/{pid}/phone_keyboard_features/phone_keyboard_{language}_{provider_key}.csv
    - data/processed/features/{pid}/phone_keyboard.csv
    ```

Parameters description for `[PHONE_KEYBOARD][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`     | Set to `True` to extract `PHONE_KEYBOARD` features from the `RAPIDS` provider|
|`[FEATURES]`    | Features to be computed, see table below
|`[TYPING_SESSION_DURATION]` | Minimum seconds to detect the end of a typing session. A session begins with any keypress and finishes until `TYPING_SESSION_DURATION` seconds (by default, 5 seconds) have elapsed since the last key was pressed or the application that the user was typing on changes.

Features description for `[PHONE_KEYBOARD][PROVIDERS][RAPIDS]`:

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
|sessioncount                                            | -    |Number of typing sessions in a time segment. Type sesssions are detected based on `TYPING_SESSION_DURATION` parameter.
|averagesessionlength                                           | milliseconds          | Average length of all sessions in a time segment instance
|averageinterkeydelay                                                |milliseconds        |The average time between keystrokes measured in milliseconds.
|changeintextlengthlessthanminusone                                                 |         | Number of times a keyboard typing or swiping event changed the length of the current text to less than one fewer character.
|changeintextlengthequaltominusone                                                 |         | Number of times a keyboard typing or swiping event changed the length of the current text in exactly one fewer character.
|changeintextlengthequaltoone                                                 |         | Number of times a keyboard typing or swiping event changed the length of the current text in exactly one more character.
|changeintextlengthmorethanone                                                 |         | Number of times a keyboard typing or swiping event changed the length of the current text to more than one character.
|maxtextlength                                                      |        | Length in characters of the longest sentence(s) contained in the typing text box of any app during the time segment.
|lastmessagelength                                                  |       | Length of the last text in characters of the sentence(s) contained in the typing text box of any app during the time segment.
|totalkeyboardtouches                                               |       | Average number of typing events across all sessions in a time segment instance.

!!! note
    1. We did not find a reliable way to distinguish between AutoCorrect or AutoComplete changes, since both can be applied with a single touch or swipe event and can decrease or increase the length of the text by an arbitrary number of characters.
    2. The default value (5 seconds) of `TYPING_SESSION_DURATION` parameter is based on empirical tests with our datasets. It could be updated as needed. Vesel et al. (check this [paper](https://academic.oup.com/jamia/article/27/7/1007/5848291)) used 8 seconds instead.
    