# Sentiment Data Analysis

Sensor parameters description for `[PHONE_SENTIMENT]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Data stream [container](../../datastreams/data-streams-introduction/) (database table, CSV file, etc.) where the messages data is stored

## RAPIDS provider

!!! info "Available time segments and platforms"
    - Available for all time segments
    - Available for Android only

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/phone_sentiment_raw.csv
    - data/raw/{pid}/phone_sentiment_with_datetime.csv
    - data/interim/{pid}/phone_sentiment_features/phone_sentiment_{language}_{provider_key}.csv
    - data/processed/features/{pid}/phone_sentiment.csv
    - ?? data/processed/features/{pid}/all_sensor_features.csv
    ```


Parameters description for `[PHONE_SENTIMENT][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`| Set to `True` to extract `PHONE_SENTIMENT` features from the `RAPIDS` provider|
|`[SCOPE]` |  List of app designations which should be included in the analysis. (Ex - `com.mozilla.firefox` or `com.android.systemui`.) Special designation `per_app` indicates all apps in the data should be done individually (so listing all desired apps is not necessary.) Special designation `all_apps` indcates that along with any app specific processing indicated, a run of the data not divided up by app should also be done.
|`[LEXICA_INCLUDED]` | List of lexicons that are used in the file.  Current logic will normalize vs. the total words for the entry, however, further work will change AWARE to collect lexica specific data to normalize against individual lexicon words.
|`[FEATURES]` | `score` is the normalized value for each lexical item.


Features description for `[PHONE_SENTIMENT][PROVIDERS][RAPIDS]`:

|Feature  |Units       |Description|
|---------|-------------|-----------|
|score    | real number | normalized value for that particular item across a particular `time_segment`.

Lexicon Citations:

|Lexicon            |Citation|
|-------------------|-------------------------------------------------------------------------------|
|life satisfaction|Jaidka, K., Giorgi, S., Schwartz, H. A., Kern, M. L., Ungar, L. H., &amp; Eichstaedt, J. C. (2020). Estimating geographic subjective well-being from Twitter: A comparison of dictionary and data-driven language methods. Proceedings of the National Academy of Sciences, 117(19), 10165–10171. https://doi.org/10.1073/pnas.1906364117 
|happiness|Giorgi, S., Guntuku, S. C., Eichstaedt, J. C., Pajot, C., Schwartz, H. A., &amp; Ungar, L. H. (n.d.). Well-being depends on social comparison: Hierarchical models of Twitter language suggest that richer neighbors make you less happy. Proceedings of the International AAAI Conference on Web and Social Media. Retrieved June 16, 2022, from https://ojs.aaai.org/index.php/ICWSM/article/view/18132 
|stress|Guntuku, S. C., Buffone, A., Jaidka, K., Eichstaedt, J. C., &amp; Ungar, L. H. (n.d.). Understanding and measuring psychological stress using social media. Proceedings of the International AAAI Conference on Web and Social Media. Retrieved June 16, 2022, from https://ojs.aaai.org/index.php/ICWSM/article/view/3223 
|loneliness|Guntuku, S. C., Schneider, R., Pelullo, A., Young, J., Wong, V., Ungar, L., Polsky, D., Volpp, K. G., &amp; Merchant, R. (2019, November 1). Studying expressions of loneliness in individuals using Twitter: An observational study. BMJ Open. Retrieved June 16, 2022, from https://bmjopen.bmj.com/content/9/11/e030355 (http://dx.doi.org/10.1136/bmjopen-2019-030355)
|affect|Preoţiuc-Pietro, D., Schwartz, H. A., Park, G., Eichstaedt, J., Kern, M., Ungar, L., &amp; Shulman, E. (n.d.). Modelling valence and arousal in Facebook posts. ACL Anthology. Retrieved June 16, 2022, from https://aclanthology.org/W16-0404/ 
|politeness|Li, M., Hickman, L., Tay, L., Ungar, L., &amp; Guntuku, S. C. (2020, October 1). Studying politeness across cultures using English Twitter and Mandarin weibo. Proceedings of the ACM on Human-Computer Interaction. Retrieved June 16, 2022, from https://dl.acm.org/doi/abs/10.1145/3415190  