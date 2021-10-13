# Phone Locations

Sensor parameters description for `[PHONE_LOCATIONS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Data stream [container](../../datastreams/data-streams-introduction/) (database table, CSV file, etc.) where the location data is stored
|`[LOCATIONS_TO_USE]`| Type of location data to use, one of `ALL`, `GPS`, `ALL_RESAMPLED` or `FUSED_RESAMPLED`. This filter is based on the `provider` column of the locations table, `ALL` includes every row, `GPS` only includes rows where the provider is gps, `ALL_RESAMPLED` includes all rows after being resampled, and `FUSED_RESAMPLED` only includes rows where the provider is fused after being resampled.
|`[FUSED_RESAMPLED_CONSECUTIVE_THRESHOLD]`| If `ALL_RESAMPLED` or `FUSED_RESAMPLED` is used, the original fused data has to be resampled. A location row is resampled to the next valid timestamp (see the Assumptions/Observations below) only if the time difference between them is less or equal than this threshold (in minutes).
|`[FUSED_RESAMPLED_TIME_SINCE_VALID_LOCATION]`| If `ALL_RESAMPLED` or `FUSED_RESAMPLED` is used, the original fused data has to be resampled. A location row is resampled at most for this long (in minutes).
|`[ACCURACY_LIMIT]` | An integer in meters, any location rows with an accuracy higher or equal than this is dropped. This number means there's a 68% probability the actual location is within this radius.

!!! note "Assumptions/Observations"
    **Types of location data to use**
    Android and iOS clients can collect location coordinates through the phone's GPS, the network cellular towers around the phone, or Google's fused location API. 
    
    - If you want to use only the GPS provider, set `[LOCATIONS_TO_USE]` to `GPS`
    - If you want to use all providers, set `[LOCATIONS_TO_USE]` to `ALL`
    - If you collected location data from different providers, including the fused API, use `ALL_RESAMPLED`
    - If your mobile client was configured to use fused location only or want to focus only on this provider, set `[LOCATIONS_TO_USE]` to `FUSED_RESAMPLED`.
    
    `ALL_RESAMPLED` and `FUSED_RESAMPLED` take the original location coordinates and replicate each pair forward in time as long as the phone was sensing data as indicated by the joined timestamps of [`[PHONE_DATA_YIELD][SENSORS]`](../phone-data-yield/). This is done because Google's API only logs a new location coordinate pair when it is sufficiently different in time or space from the previous one and because GPS and network providers can log data at variable rates.

    There are two parameters associated with resampling fused location.
    
    1. `FUSED_RESAMPLED_CONSECUTIVE_THRESHOLD` (in minutes, default 30) controls the maximum gap between any two coordinate pairs to replicate the last known pair. For example, participant A's phone did not collect data between 10.30 am and 10:50 am and between 11:05am and 11:40am, the last known coordinate pair is replicated during the first period but not the second. In other words, we assume that we cannot longer guarantee the participant stayed at the last known location if the phone did not sense data for more than 30 minutes. 
    2. `FUSED_RESAMPLED_TIME_SINCE_VALID_LOCATION` (in minutes, default 720 or 12 hours) stops the last known fused location from being replicated longer than this threshold even if the phone was sensing data continuously. For example, participant A went home at 9 pm, and their phone was sensing data without gaps until 11 am the next morning, the last known location is replicated until 9 am. 
    
    If you have suggestions to modify or improve this resampling, let us know.

## BARNETT provider

These features are based on the original open-source implementation by [Barnett et al](../../citation#barnett-locations) and some features created by [Canzian et al](../../citation#barnett-locations).


!!! info "Available time segments and platforms"
    - Available only for segments that start at 00:00:00 and end at 23:59:59 of the same or a different day (daily, weekly, weekend, etc.)
    - Available for Android and iOS

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/phone_locations_raw.csv
    - data/interim/{pid}/phone_locations_processed.csv
    - data/interim/{pid}/phone_locations_processed_with_datetime.csv
    - data/interim/{pid}/phone_locations_barnett_daily.csv
    - data/interim/{pid}/phone_locations_features/phone_locations_{language}_{provider_key}.csv
    - data/processed/features/{pid}/phone_locations.csv
    ```


Parameters description for `[PHONE_LOCATIONS][PROVIDERS][BARNETT]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`| Set to `True` to extract `PHONE_LOCATIONS` features from the `BARNETT` provider|
|`[FEATURES]` |         Features to be computed, see table below
|`[IF_MULTIPLE_TIMEZONES]` |    Currently, `USE_MOST_COMMON` is the only value supported. If the location data for a participant belongs to multiple time zones, we select the most common because Barnett's algorithm can only handle one time zone 
|`[MINUTES_DATA_USED]` |    Set to `True` to include an extra column in the final location feature file containing the number of minutes used to compute the features on each time segment. Use this for quality control purposes; the more data minutes exist for a period, the more reliable its features should be. For fused location, a single minute can contain more than one coordinate pair if the participant is moving fast enough.



Features description for `[PHONE_LOCATIONS][PROVIDERS][BARNETT]` adapted from [Beiwe Summary Statistics](http://wiki.beiwe.org/wiki/Summary_Statistics):

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
|hometime                              |minutes     | Time at home. Time spent at home in minutes. Home is the most visited significant location between 8 pm and 8 am, including any pauses within a 200-meter radius.
|disttravelled                         |meters      | Total distance traveled over a day (flights).
|rog                                   |meters      | The Radius of Gyration (rog) is a measure in meters of the area covered by a person over a day. A centroid is calculated for all the places (pauses) visited during a day, and a weighted distance between all the places and that centroid is computed. The weights are proportional to the time spent in each place.
|maxdiam                               |meters      | The maximum diameter is the largest distance between any two pauses.
|maxhomedist                           |meters      | The maximum distance from home in meters.
|siglocsvisited                        |locations   | The number of significant locations visited during the day. Significant locations are computed using k-means clustering over pauses found in the whole monitoring period. The number of clusters is found iterating k from 1 to 200 stopping until the centroids of two significant locations are within 400 meters of one another.
|avgflightlen                          |meters      | Mean length of all flights.
|stdflightlen                          |meters      | Standard deviation of the length of all flights.
|avgflightdur                          |seconds     | Mean duration of all flights.
|stdflightdur                           |seconds     | The standard deviation of the duration of all flights. 
|probpause                              |     -      | The fraction of a day spent in a pause (as opposed to a flight)
|siglocentropy                          |nats        | Shannon's entropy measurement is based on the proportion of time spent at each significant location visited during a day.
|circdnrtn                              |      -     |   A continuous metric quantifying a person's circadian routine that can take any value between 0 and 1, where 0 represents a daily routine completely different from any other sensed days and 1 a routine the same as every other sensed day.
|wkenddayrtn                            |       -    | Same as circdnrtn but computed separately for weekends and weekdays.

!!! note "Assumptions/Observations"
    **Multi day segment features**
    Barnett's features are only available on time segments that span entire days (00:00:00 to 23:59:59). Such segments can be one-day long (daily) or multi-day (weekly, for example). Multi-day segment features are computed based on daily features summarized the following way:

    - sum for `hometime`, `disttravelled`, `siglocsvisited`, and `minutes_data_used`
    - max for `maxdiam`, and `maxhomedist`
    - mean for `rog`, `avgflightlen`, `stdflightlen`, `avgflightdur`, `stdflightdur`, `probpause`, `siglocentropy`, `circdnrtn`, `wkenddayrtn`, and `minsmissing`

    **Computation speed**
    The process to extract these features can be slow compared to other sensors and providers due to the required simulation.

    **How are these features computed?**
    These features are based on a Pause-Flight model. A pause is defined as a mobility trace (location pings) within a certain duration and distance (by default, 300 seconds and 60 meters). A flight is any mobility trace between two pauses. Data is resampled and imputed before the features are computed. See [Barnett et al](../../citation#barnett-locations) for more information. In RAPIDS, we only expose one parameter for these features (accuracy limit). You can change other parameters in `src/features/phone_locations/barnett/library/MobilityFeatures.R`.

    **Significant Locations**
    Significant locations are determined using K-means clustering on pauses longer than 10 minutes. The number of clusters (K) is increased until no two clusters are within 400 meters from each other. After this, pauses within a certain range of a cluster (200 meters by default) count as a visit to that significant location. This description was adapted from the Supplementary Materials of [Barnett et al](../../citation#barnett-locations).

    **The Circadian Calculation**
    For a detailed description of how this is calculated, see [Canzian et al](../../citation#barnett-locations).

## DORYAB provider

These features are based on the original implementation by [Doryab et al.](../../citation#doryab-locations).


!!! info "Available time segments and platforms"
    - Available for all time segments
    - Available for Android and iOS

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/phone_locations_raw.csv
    - data/interim/{pid}/phone_locations_processed.csv
    - data/interim/{pid}/phone_locations_processed_with_datetime.csv
    - data/interim/{pid}/phone_locations_processed_with_datetime_with_doryab_columns_episodes.csv
    - data/interim/{pid}/phone_locations_processed_with_datetime_with_doryab_columns_episodes_resampled.csv
    - data/interim/{pid}/phone_locations_processed_with_datetime_with_doryab_columns_episodes_resampled_with_datetime.csv
    - data/interim/{pid}/phone_locations_features/phone_locations_{language}_{provider_key}.csv
    - data/processed/features/{pid}/phone_locations.csv
    ```


Parameters description for `[PHONE_LOCATIONS][PROVIDERS][DORYAB]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`| Set to `True` to extract `PHONE_LOCATIONS` features from the `DORYAB` provider|
|`[FEATURES]` |         Features to be computed, see table below
| `[DBSCAN_EPS]`             | The maximum distance in meters between two samples for one to be considered as in the neighborhood of the other. This is not a maximum bound on the distances of points within a cluster. This is the most important DBSCAN parameter to choose appropriately for your data set and distance function.
| `[DBSCAN_MINSAMPLES]`      | The number of samples (or total weight) in a neighborhood for a point to be considered as a core point of a cluster. This includes the point itself.
| `[THRESHOLD_STATIC]`       | It is the threshold value in km/hr which labels a row as Static or Moving.
| `[MAXIMUM_ROW_GAP]`   | The maximum gap (in seconds) allowed between any two consecutive rows for them to be considered part of the same displacement. If this threshold is too high, it can throw speed and distance calculations off for periods when the phone was not sensing. This value must be larger than your GPS sampling interval when `[LOCATIONS_TO_USE]` is `ALL` or `GPS`, otherwise all the stationary-related features will be NA. If `[LOCATIONS_TO_USE]` is `ALL_RESAMPLED` or `FUSED_RESAMPLED`, you can use the default value as every row will be resampled at 1-minute intervals.
| `[MINUTES_DATA_USED]`     | Set to `True` to include an extra column in the final location feature file containing the number of minutes used to compute the features on each time segment. Use this for quality control purposes; the more data minutes exist for a period, the more reliable its features should be. For fused location, a single minute can contain more than one coordinate pair if the participant is moving fast enough.
| `[CLUSTER_ON]`             | Set this flag to `PARTICIPANT_DATASET` to create clusters based on the entire participant's dataset or to `TIME_SEGMENT` to create clusters based on all the instances of the corresponding time segment (e.g. all mornings) or to `TIME_SEGMENT_INSTANCE` to create clusters based on a single instance (e.g. 2020-05-20's morning).
|`[INFER_HOME_LOCATION_STRATEGY]`          | The strategy applied to infer home locations. Set to `DORYAB_STRATEGY` to infer one home location for the entire dataset of each participant or to `SUN_LI_VEGA_STRATEGY` to infer one home location per day per participant. See Observations below to know more.
|`[MINIMUM_DAYS_TO_DETECT_HOME_CHANGES]`   | The minimum number of consecutive days a new home location candidate has to repeat before it is considered the participant's new home. This parameter will be used only when `[INFER_HOME_LOCATION_STRATEGY]` is set to `SUN_LI_VEGA_STRATEGY`.
| `[CLUSTERING_ALGORITHM]`   | The original Doryab et al. implementation uses `DBSCAN`, `OPTICS` is also available with similar (but not identical) clustering results and lower memory consumption.
| `[RADIUS_FOR_HOME]`        | All location coordinates within this distance (meters) from the home location coordinates are considered a homestay (see `timeathome` feature).


Features description for `[PHONE_LOCATIONS][PROVIDERS][DORYAB]`:

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
|locationvariance                                            |$meters^2$    |The sum of the variances of the latitude and longitude columns. 
|loglocationvariance                                           | -          | Log of the sum of the variances of the latitude and longitude columns.
|totaldistance                                                |meters        |Total distance traveled in a time segment using the haversine formula.
|avgspeed                                                 |km/hr         |Average speed in a time segment considering only the instances labeled as Moving. This feature is 0 when the participant is stationary during a time segment.
|varspeed                                                      |km/hr         |Speed variance in a time segment considering only the instances labeled as Moving. This feature is 0 when the participant is stationary during a time segment.
|{--circadianmovement--}                                      |-             | Deprecated, see Observations below. \ "It encodes the extent to which a person's location patterns follow a 24-hour circadian cycle.\" [Doryab et al.](../../citation#doryab-locations).
|numberofsignificantplaces                                    |places        |Number of significant locations visited. It is calculated using the DBSCAN/OPTICS clustering algorithm which takes in EPS and MIN_SAMPLES as parameters to identify clusters. Each cluster is a significant place.
|numberlocationtransitions                                    |transitions   |Number of movements between any two clusters in a time segment.
|radiusgyration                                               |meters        |Quantifies the area covered by a participant
|timeattop1location                                           |minutes       |Time spent at the most significant location.
|timeattop2location                                           |minutes       |Time spent at the 2nd most significant location.
|timeattop3location                                           |minutes       |Time spent at the 3rd most significant location. 
|movingtostaticratio                                          | -   |  Ratio between stationary time and total location sensed time. A lat/long coordinate pair is labeled as stationary if its speed (distance/time) to the next coordinate pair is less than 1km/hr. A higher value represents a more stationary routine.
|outlierstimepercent                                          | -   | Ratio between the time spent in non-significant clusters divided by the time spent in all clusters (stationary time. Only stationary samples are clustered). A higher value represents more time spent in non-significant clusters.
|maxlengthstayatclusters                                      |minutes       |Maximum time spent in a cluster (significant location).
|minlengthstayatclusters                                      |minutes       |Minimum time spent in a cluster (significant location).
|avglengthstayatclusters                                      |minutes       |Average time spent in a cluster (significant location).
|stdlengthstayatclusters                                      |minutes       |Standard deviation of time spent in a cluster (significant location).
|locationentropy                                              |nats          |Shannon Entropy computed over the row count of each cluster (significant location), it is higher the more rows belong to a cluster (i.e., the more time a participant spent at a significant location).
|normalizedlocationentropy                                    |nats          |Shannon Entropy computed over the row count of each cluster (significant location) divided by the number of clusters; it is higher the more rows belong to a cluster (i.e., the more time a participant spent at a significant location).
|timeathome                                                   |minutes       | Time spent at home (see Observations below for a description on how we compute home).
|homelabel                                                    |-             | An integer that represents a different home location. It will be a constant number (1) for all participants when `[INFER_HOME_LOCATION_STRATEGY]` is set to `DORYAB_STRATEGY` or an incremental index if the strategy is set to `SUN_LI_VEGA_STRATEGY`.

!!! note "Assumptions/Observations"
    **Significant Locations Identified**
    Significant locations are determined using `DBSCAN` or `OPTICS` clustering on locations that a participant visited over the course of the period of data collection. The most significant location is the place where the participant stayed for the longest time.

    **Circadian Movement Calculation**
    Note Feb 3 2021. It seems the implementation of this feature is not correct; we suggest not to use this feature until a fix is in place. For a detailed description of how this should be calculated, see [Saeb et al](https://pubmed.ncbi.nlm.nih.gov/28344895/).

    **Fine-Tuning Clustering Parameters**
    Based on an experiment where we collected fused location data for 7 days with a mean accuracy of 86 & SD of 350.874635, we determined that `EPS/MAX_EPS`=100 produced closer clustering results to reality. Higher values (>100) missed out on some significant places, like a short grocery visit, while lower values (<100) picked up traffic lights and stop signs while driving as significant locations. We recommend you set `EPS` based on your location data's accuracy (the more accurate your data is, the lower you should be able to set EPS).

    **Duration Calculation**
    To calculate the time duration component for our features, we compute the difference between consecutive rows' timestamps to take into account sampling rate variability. If this time difference is larger than a threshold (300 seconds by default), we replace it with NA and label that row as Moving.

    **Home location**

    - `DORYAB_STRATEGY`: home is calculated using all location data of a participant between 12 am and 6 am, then applying a clustering algorithm (`DBSCAN` or `OPTICS`) and considering the center of the biggest cluster home for that participant.
    
    - `SUN_LI_VEGA_STRATEGY`: home is calculated using all location data of a participant between 12 am and 6 am, then applying a clustering algorithm (`DBSCAN` or `OPTICS`). The following steps are used to infer the home location per day for that participant:
        
        1.  if there are records within [03:30:00, 04:30:00] for that night:<br>
                &nbsp;&nbsp;&nbsp;&nbsp;we choose the most common cluster during that period as a home candidate for that day.<br>
            elif there are records within [midnight, 03:30:00) for that night:<br>
                &nbsp;&nbsp;&nbsp;&nbsp;we choose the last valid cluster during that period as a home candidate for that day.<br>
            elif there are records within (04:30:00, 06:00:00] for that night:<br>
                &nbsp;&nbsp;&nbsp;&nbsp;we choose the first valid cluster during that period as a home candidate for that day.<br>
            else:<br>
                &nbsp;&nbsp;&nbsp;&nbsp;the home location is NA (missing) for that day.

        2. If the count of consecutive days with the same candidate home location cluster label is larger or equal to `[MINIMUM_DAYS_TO_DETECT_HOME_CHANGES]`,
            the candidate will be regarded as the home cluster; otherwise, the home cluster will be the last valid day's cluster.
            If there are no valid clusters before that day, the first home location in the days after is used.

    **Clustering algorithms**
    [`DBSCAN`](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html) and [`OPTICS`](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.OPTICS.html#r2c55e37003fe-1) algorithms are available currently. Duplicated locations are discarded while clustering. The `DBSCAN` algorithm takes the time spent at each location into consideration. However, the `OPTICS` algorithm ignores it as it is not supported in the current [scikit-learn](https://github.com/scikit-learn/scikit-learn/issues/12394) implementation.
