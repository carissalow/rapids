# Phone Keyboard

Sensor parameters description for `[PHONE_KEYBOARD]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Data stream [container](../../datastreams/data-streams-introduction/) (database table, CSV file, etc.) where the keyboard data is stored
|`[NEW_SESSION_KEYSTROKE_MIN_PAUSE_LENGTH]`| Minimum length of time between keystrokes (in ms) to indicate that a new session should be created.  Converted to a parameter after literature became split on how long this duration should be. <sup>(1)</sup> Original Rapids hardcode value was 5000.
|`[MIN_KEYSTROKE_PER_SESSION]`| Minimum number of keystrokes to be considered a session.  Converted to a parameter after literature split on the minimum number of required keystrokes. <sup>(2)</sup>  Rapids originally did not have a minimum (now represented by a 0 value)

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
    

<sup>(1)</sup> Vesel, C., Rashidisabet, H., Zulueta, J., Stange, J. P., Duffecy, J., Hussain, F., Piscitello, A., Bark, J., Langenecker, S. A., Young, S., Mounts, E., Omberg, L., Nelson, P. C., Moore, R. C., Koziol, D., Bourne, K., Bennett, C. C., Ajilore, O., Demos, A. P., &amp; Leow, A. (2020). Effects of mood and aging on keystroke dynamics metadata and their diurnal patterns in a large open-science sample: A BIAFFECT IOS study. Journal of the American Medical Informatics Association, 27(7), 1007–1018. https://doi.org/10.1093/jamia/ocaa057 <br>
 vs <br>
Zulueta1, J., Piscitello1, A., Rasic1, M., Easter1, R., Babu2, P., Langenecker1, S. A., McInnis2, M., Ajilore1, O., Nelson1, P. C., Ryan2, K., &amp; Leow1, A. (n.d.). Predicting mood disturbance severity with mobile phone keystroke metadata: A BIAFFECT Digital Phenotyping Study. Journal of Medical Internet Research. Retrieved June 22, 2022, from https://www.jmir.org/2018/7/e241/ 

<sup>(2)</sup> Mastoras, R.-E., Iakovakis, D., Hadjidimitriou, S., Charisis, V., Kassie, S., Alsaadi, T., Khandoker, A., &amp; Hadjileontiadis, L. J. (2019, September 16). Touchscreen typing pattern analysis for remote detection of the depressive tendency. Scientific reports. Retrieved June 22, 2022, from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6746713/  <br>
vs <br>
Cao, B., Zheng, L., Zhang, C., Yu, P. S., Piscitello, A., Zulueta, J., Ajilore, O., Ryan, K., &amp; Leow, A. D. (2017, August 1). Deepmood: Proceedings of the 23rd ACM SIGKDD International Conference on Knowledge Discovery and data mining. ACM Conferences. Retrieved June 22, 2022, from https://dl.acm.org/doi/10.1145/3097983.3098086  