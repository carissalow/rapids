# Phone Applications Crashes

Sensor parameters description for `[PHONE_APPLICATIONS_CRASHES]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Data stream [container](../../datastreams/data-streams-introduction/) (database table, CSV file, etc.) where the applications crashes data is stored
|`[APPLICATION_CATEGORIES][CATALOGUE_SOURCE]` | `FILE` or `GOOGLE`. If `FILE`, app categories (genres) are read from `[CATALOGUE_FILE]`. If `[GOOGLE]`, app categories (genres) are scrapped from the Play Store
|`[APPLICATION_CATEGORIES][CATALOGUE_FILE]` | CSV file with a `package_name` and `genre` column. By default we provide the catalogue created by [Stachl et al](../../citation#stachl-applications-crashes) in `data/external/stachl_application_genre_catalogue.csv`
|`[APPLICATION_CATEGORIES][UPDATE_CATALOGUE_FILE]` | if `[CATALOGUE_SOURCE]` is equal to `FILE`, this flag signals whether or not to update `[CATALOGUE_FILE]`, if `[CATALOGUE_SOURCE]` is equal to `GOOGLE` all scraped genres will be saved to `[CATALOGUE_FILE]`
|`[APPLICATION_CATEGORIES][SCRAPE_MISSING_CATEGORIES]` | This flag signals whether or not to scrape categories (genres) missing from the `[CATALOGUE_FILE]`. If `[CATALOGUE_SOURCE]` is equal to `GOOGLE`, all genres are scraped anyway (this flag is ignored)

!!! note
    No feature providers have been implemented for this sensor yet, however you can use its key (`PHONE_APPLICATIONS_CRASHES`) to improve [`PHONE_DATA_YIELD`](../phone-data-yield) or you can [implement your own features](../add-new-features).