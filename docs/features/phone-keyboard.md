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
|sessioncount                                            | -    |Number of typing sessions in a time segment. A session begins with any keypress and finishes until 5 seconds have elapsed since the last key was pressed or the application that the user was typing on changes.
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
    We did not find a reliable way to distinguish between AutoCorrect or AutoComplete changes, since both can be applied with a single touch or swipe event and can decrease or increase the length of the text by an arbitrary number of characters.
    