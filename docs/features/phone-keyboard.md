# Phone Keyboard

Sensor parameters description for `[PHONE_KEYBOARD]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Data stream [container](../../datastreams/data-streams-introduction/) (database table, CSV file, etc.) where the keyboard data is stored

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/phone_keyboard_raw.csv
    - data/raw/{pid}/phone_keyboard_with_datetime.csv
    - data/interim/{pid}/phone_keyboard_features/phone_keyboard_{language}_{provider_key}.csv
    - data/processed/features/{pid}/phone_keyboard.csv
    ```

Features description for `[PHONE_KEYBOARD]`:

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
|sessioncount                                            | -    |Number of sessions: A session begins when a keypress is initiated and ≥5 s has elapsed since the last key was pressed. A session ends when ≥5 s has elapsed since the last key was pressed.
|averagesessionlength                                           | milliseconds          | Length of sessions in milliseconds averaged over the segment.
|averageinterkeydelay                                                |milliseconds        |The average time between keystrokes measured in milliseconds.
|changeintextlengthlessthanminusone                                                 |         | Number of times the keyboard touch changed the length of the text to less than -1.
|changeintextlengthequaltominusone                                                 |         | Number of times the keyboard touch changed the length of the text equal to -1.
|changeintextlengthequaltoone                                                 |         | Number of times the keyboard touch changed the length of the text equal to 1.
|changeintextlengthmorethanone                                                 |         | Number of times the keyboard touch changed the length of the text to more than 1.
|maxtextlength                                                      |        | Length of the biggest text in a session averaged over the sessions.
|lastmessagelength                                                  |       | Length of the last text in a session averaged over the sessions.
|uniqueapplications                                                 |       | Number of distinct applications in a session averaged over the sessions.
|totalkeyboardtouches                                               |       | Number of times keyboard was touched in a session averaged over the sessions.

!!! note
    Given the set of features, it was difficult to distinguish between AutoCorrect or AutoComplete event. Since both can be applied using a single touch and can decrease or increase the length of the text, there was no visible pattern in the raw data.
    