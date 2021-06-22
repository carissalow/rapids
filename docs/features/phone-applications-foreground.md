# Phone Applications Foreground

Sensor parameters description for `[PHONE_APPLICATIONS_FOREGROUND]` (these parameters are used by the only provider available at the moment, RAPIDS):

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Data stream [container](../../datastreams/data-streams-introduction/) (database table, CSV file, etc.) where the applications foreground data is stored
|`[APPLICATION_CATEGORIES][CATALOGUE_SOURCE]` | `FILE` or `GOOGLE`. If `FILE`, app categories (genres) are read from `[CATALOGUE_FILE]`. If `[GOOGLE]`, app categories (genres) are scrapped from the Play Store
|`[APPLICATION_CATEGORIES][CATALOGUE_FILE]` | CSV file with a `package_name` and `genre` column. By default we provide the catalogue created by [Stachl et al](../../citation#stachl-applications-foreground) in `data/external/stachl_application_genre_catalogue.csv`
|`[APPLICATION_CATEGORIES][UPDATE_CATALOGUE_FILE]` | if `[CATALOGUE_SOURCE]` is equal to `FILE`, this flag signals whether or not to update `[CATALOGUE_FILE]`, if `[CATALOGUE_SOURCE]` is equal to `GOOGLE` all scraped genres will be saved to `[CATALOGUE_FILE]`
|`[APPLICATION_CATEGORIES][SCRAPE_MISSING_CATEGORIES]` | This flag signals whether or not to scrape categories (genres) missing from the `[CATALOGUE_FILE]`. If `[CATALOGUE_SOURCE]` is equal to `GOOGLE`, all genres are scraped anyway (this flag is ignored)

## RAPIDS provider

The app category (genre) catalogue used in these features was originally created by [Stachl et al](../../citation#stachl-applications-foreground).

!!! info "Available time segments and platforms"
    - Available for all time segments
    - Available for Android only

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/phone_applications_foreground_raw.csv
    - data/raw/{pid}/phone_applications_foreground_with_datetime.csv
    - data/raw/{pid}/phone_applications_foreground_with_datetime_with_categories.csv
    - data/interim/{pid}/phone_applications_foreground_features/phone_applications_foreground_{language}_{provider_key}.csv
    - data/processed/features/{pid}/phone_applications_foreground.csv
    ```


Parameters description for `[PHONE_APPLICATIONS_FOREGROUND][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`| Set to `True` to extract `PHONE_APPLICATIONS_FOREGROUND` features from the `RAPIDS` provider|
|`[INCLUDE_EPISODE_FEATURES]`| Set to `True` to extract application episode features as well from the `RAPIDS` provider|
|`[FEATURES]` |         Features to be computed, see table below
|`[SINGLE_CATEGORIES]`     | An array of app categories to be *included* in the feature extraction computation. The special keyword `all` represents a category with all the apps from each participant. By default we use the category catalogue pointed by `[APPLICATION_CATEGORIES][CATALOGUE_FILE]` (see the Sensor parameters description table above)
|`[OWN_CATEGORIES]`   | An array of collections representing your own app categories. They key of each element is the name of the in-house category and the value is an array of the package names (apps) included in that category.
|`[MULTIPLE_CATEGORIES]`   | An array of collections representing meta-categories (a group of categories). They key of each element is the name of the `meta-category` and the value is an array of member app categories. By default we use the category catalogue pointed by `[APPLICATION_CATEGORIES][CATALOGUE_FILE]` (see the Sensor parameters description table above)
|`[SINGLE_APPS]`           | An array of apps to be *included* in the feature extraction computation. Use their package name (e.g. `com.google.android.youtube`) or the reserved keyword `top1global` (the most used app by a participant over the whole monitoring study)
|`[EXCLUDED_CATEGORIES]`   | An array of app categories to be *excluded* from the feature extraction computation. By default we use the category catalogue pointed by `[APPLICATION_CATEGORIES][CATALOGUE_FILE]` (see the Sensor parameters description table above)
|`[EXCLUDED_APPS]`         | An array of apps to be excluded from the feature extraction computation. Use their package name, for example: `com.google.android.youtube`

Features description for `[PHONE_APPLICATIONS_FOREGROUND][PROVIDERS][RAPIDS]`:

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
|count              |apps      | Number of times a single app or apps within a category were used (i.e. they were brought to the foreground either by tapping their icon or switching to it from another app)
|timeoffirstuse     |minutes   | The time in minutes between 12:00am (midnight) and the first use of a single app or apps within a category during a `time_segment`
|timeoflastuse      |minutes   | The time in minutes between 12:00am (midnight) and the last use of a single app or apps within a category during a `time_segment`
|frequencyentropy   |nats      | The entropy of the used apps within a category during a `time_segment` (each app is seen as a unique event, the more apps were used, the higher the entropy). This is especially relevant when computed over all apps. Entropy cannot be obtained for a single app
|minduration        |minutes   | For a `time_segment`, the minimum duration an application was used in minutes
|maxduration        |minutes   | For a `time_segment`, the maximum duration an application was used in minutes
|meanduration       |minutes   | For a `time_segment`, the mean duration of all the applications used in minutes
|sumduration        |minutes   | For a `time_segment`, the sum duration of all the applications used in minutes

!!! note "Assumptions/Observations"
    Features can be computed by app, by apps grouped under a single category (genre) and by multiple categories grouped together (meta-categories). For example, we can get features for `Facebook` (single app), for `Social Network` apps (a category including Facebook and other social media apps) or for `Social` (a meta-category formed by `Social Network` and `Social Media Tools` categories).

    Apps installed by default like YouTube are considered systems apps on some phones. We do an exact match to exclude apps where "genre" == `EXCLUDED_CATEGORIES` or "package_name" == `EXCLUDED_APPS`.

    We provide three ways of classifying and app within a category (genre): a) by automatically scraping its official category from the Google Play Store, b) by using the catalogue created by Stachl et al. which we provide in RAPIDS (`data/external/stachl_application_genre_catalogue.csv`), or c) by manually creating a personalized catalogue. You can choose a, b or c by modifying `[APPLICATION_GENRES]` keys and values (see the Sensor parameters description table above).
    
    The application episodes are calculated using the application foreground and screen unlock episode data. An application episode starts when the application is launched and ends when either, a new application is launched or the screen is locked.
