# Phone Conversation

Sensor parameters description for `[PHONE_CONVERSATION]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER][ANDROID]`| Data stream [container](../../datastreams/data-streams-introduction/) (database table, CSV file, etc.) where the conversation data from Android devices is stored (the AWARE client saves this data on different tables for Android and iOS)
|`[CONTAINER][IOS]`| Data stream [container](../../datastreams/data-streams-introduction/) (database table, CSV file, etc.) where the conversation data from iOS devices is stored (the AWARE client saves this data on different tables for Android and iOS)

## RAPIDS provider

!!! info "Available time segments and platforms"
    - Available for all time segments
    - Available for Android only

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/phone_conversation_raw.csv
    - data/raw/{pid}/phone_conversation_with_datetime.csv
    - data/interim/{pid}/phone_conversation_features/phone_conversation_{language}_{provider_key}.csv
    - data/processed/features/{pid}/phone_conversation.csv
    ```


Parameters description for `[PHONE_CONVERSATION][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`| Set to `True` to extract `PHONE_CONVERSATION` features from the `RAPIDS` provider|
|`[FEATURES]` |         Features to be computed, see table below
|`[RECORDING_MINUTES]` | Minutes the plugin was recording audio (default 1 min)
|`[PAUSED_MINUTES]` |  Minutes the plugin was NOT recording audio (default 3 min)


Features description for `[PHONE_CONVERSATION][PROVIDERS][RAPIDS]`:

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
| minutessilence          | minutes | Minutes labeled as silence                                                                                                                                                                 |
| minutesnoise            | minutes | Minutes labeled as noise                                                                                                                                                                   |
| minutesvoice            | minutes | Minutes labeled as voice                                                                                                                                                                   |
| minutesunknown          | minutes | Minutes labeled as unknown                                                                                                                                                                 |
| sumconversationduration | minutes | Total duration of all conversations                                                                                                                                                        |
| maxconversationduration | minutes | Longest duration of all conversations                                                                                                                                                      |
| minconversationduration | minutes | Shortest duration of all conversations                                                                                                                                                     |
| avgconversationduration | minutes | Average duration of all conversations                                                                                                                                                      |
| sdconversationduration  | minutes | Standard Deviation of the duration of all conversations                                                                                                                                    |
| timefirstconversation   | minutes | Minutes since midnight when the first conversation for a time segment was detected                                                                                                          |
| timelastconversation    | minutes | Minutes since midnight when the last conversation for a time segment was detected                                                                                                           |
| noisesumenergy          | L2-norm | Sum of all energy values when inference is noise                                                                                                                                           |
| noiseavgenergy          | L2-norm | Average of all energy values when inference is noise                                                                                                                                       |
| noisesdenergy           | L2-norm | Standard Deviation of all energy values when inference is noise                                                                                                                            |
| noiseminenergy          | L2-norm | Minimum of all energy values when inference is noise                                                                                                                                       |
| noisemaxenergy          | L2-norm | Maximum of all energy values when inference is noise                                                                                                                                       |
| voicesumenergy          | L2-norm | Sum of all energy values when inference is voice                                                                                                                                           |
| voiceavgenergy          | L2-norm | Average of all energy values when inference is voice                                                                                                                                       |
| voicesdenergy           | L2-norm | Standard Deviation of all energy values when inference is voice                                                                                                                            |
| voiceminenergy          | L2-norm | Minimum of all energy values when inference is voice                                                                                                                                       |
| voicemaxenergy          | L2-norm | Maximum of all energy values when inference is voice                                                                                                                                       |
| silencesensedfraction   |   -      | Ratio between minutessilence and the sum of (minutessilence, minutesnoise, minutesvoice, minutesunknown)                                                                                   |
| noisesensedfraction     |   -      | Ratio between minutesnoise and the sum of (minutessilence, minutesnoise, minutesvoice, minutesunknown)                                                                                     |
| voicesensedfraction     |   -      | Ratio between minutesvoice and the sum of (minutessilence, minutesnoise, minutesvoice, minutesunknown)                                                                                     |
| unknownsensedfraction   |   -      | Ratio between minutesunknown and the sum of (minutessilence, minutesnoise, minutesvoice, minutesunknown)                                                                                   |
| silenceexpectedfraction |   -      | Ration between minutessilence and the number of minutes that in  theory should have been sensed based on the record and pause cycle of  the plugin (1440 / recordingMinutes+pausedMinutes) |
| noiseexpectedfraction   |   -      | Ration between minutesnoise and the number of minutes that in theory  should have been sensed based on the record and pause cycle of the  plugin (1440 / recordingMinutes+pausedMinutes)   |
| voiceexpectedfraction   |   -      | Ration between minutesvoice and the number of minutes that in theory  should have been sensed based on the record and pause cycle of the  plugin (1440 / recordingMinutes+pausedMinutes)   |
| unknownexpectedfraction |   -      | Ration between minutesunknown and the number of minutes that in  theory should have been sensed based on the record and pause cycle of  the plugin (1440 / recordingMinutes+pausedMinutes) |

!!! note "Assumptions/Observations"
    1. The timestamp of conversation rows in iOS is in seconds so we convert it to milliseconds to match Android's format
