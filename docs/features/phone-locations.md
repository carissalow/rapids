# Phone Locations

Sensor parameters description for `[PHONE_LOCATIONS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[TABLE]`| Database table where the location data is stored
|`[LOCATIONS_TO_USE]`| Type of location data to use, one of `ALL`, `GPS` or `FUSED_RESAMPLED`. This filter is based on the `provider` column of the AWARE locations table, `ALL` includes every row, `GPS` only includes rows where provider is gps, and `FUSED_RESAMPLED` only includes rows where provider is fused after being resampled.
|`[FUSED_RESAMPLED_CONSECUTIVE_THRESHOLD]`| if `FUSED_RESAMPLED` is used, the original fused data has to be resampled, a location row will be resampled to the next valid timestamp (see the Assumptions/Observations below) only if the time difference between them is less or equal than this threshold (in minutes).
|`[FUSED_RESAMPLED_TIME_SINCE_VALID_LOCATION]`| if `FUSED_RESAMPLED` is used, the original fused data has to be resampled, a location row will be resampled at most for this long (in minutes)

!!! note "Assumptions/Observations"
    **Types of location data to use**
    AWARE Android and iOS clients can collect location coordinates through the phone\'s GPS, the network cellular towers around the phone or Google\'s fused location API. If you want to use only the GPS provider set `[LOCATIONS_TO_USE]` to `GPS`, if you want to use all providers (not recommended due to the difference in accuracy) set `[LOCATIONS_TO_USE]` to `ALL`, if your AWARE client was configured to use fused location only or want to focus only on this provider, set `[LOCATIONS_TO_USE]` to `RESAMPLE_FUSED`. `RESAMPLE_FUSED` takes the original fused location coordinates and replicates each pair forward in time as long as the phone was sensing data as indicated by the joined timestamps of [`[PHONE_DATA_YIELD][SENSORS]`](../phone-data-yield/), this is done because Google\'s API only logs a new location coordinate pair when it is sufficiently different in time or space from the previous one.

    There are two parameters associated with resampling fused location. `FUSED_RESAMPLED_CONSECUTIVE_THRESHOLD` (in minutes, default 30) controls the maximum gap between any two coordinate pairs to replicate the last known pair (for example, participant A\'s phone did not collect data between 10.30am and 10:50am and between 11:05am and 11:40am, the last known coordinate pair will be replicated during the first period but not the second, in other words, we assume that we cannot longer guarantee the participant stayed at the last known location if the phone did not sense data for more than 30 minutes). `FUSED_RESAMPLED_TIME_SINCE_VALID_LOCATION` (in minutes, default 720 or 12 hours) stops the last known fused location from being replicated longer that this threshold even if the phone was sensing data continuously (for example, participant A went home at 9pm and their phone was sensing data without gaps until 11am the next morning, the last known location will only be replicated until 9am). If you have suggestions to modify or improve this resampling, let us know.

## BARNETT provider

These features are based on the original open-source implementation by [Barnett et al](../../citation#barnett-locations) and some features created by [Canzian et al](../../citation#barnett-locations).


!!! info "Available time segments and platforms"
    - Available only for segments that start at 00:00:00 and end at 23:59:59 of the same day (daily segments)
    - Available for Android and iOS

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/phone_locations_raw.csv
    - data/interim/{pid}/phone_locations_processed.csv
    - data/interim/{pid}/phone_locations_processed_with_datetime.csv
    - data/interim/{pid}/phone_locations_features/phone_locations_{language}_{provider_key}.csv
    - data/processed/features/{pid}/phone_locations.csv
    ```


Parameters description for `[PHONE_LOCATIONS][PROVIDERS][BARNETT]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`| Set to `True` to extract `PHONE_LOCATIONS` features from the `BARNETT` provider|
|`[FEATURES]` |         Features to be computed, see table below
|`[ACCURACY_LIMIT]` |   An integer in meters, any location rows with an accuracy higher than this will be dropped. This number means there's a 68% probability the true location is within this radius
|`[TIMEZONE]` |    Timezone where the location data was collected. By default points to the one defined in the [Configuration](../../setup/configuration#timezone-of-your-study)
|`[MINUTES_DATA_USED]` |    Set to `True` to include an extra column in the final location feature file containing the number of minutes used to compute the features on each time segment. Use this for quality control purposes, the more data minutes exist for a period, the more reliable its features should be. For fused location, a single minute can contain more than one coordinate pair if the participant is moving fast enough.



Features description for `[PHONE_LOCATIONS][PROVIDERS][BARNETT]` adapted from [Beiwe Summary Statistics](http://wiki.beiwe.org/wiki/Summary_Statistics):

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
|hometime                              |minutes     | Time at home. Time spent at home in minutes. Home is the most visited significant location between 8 pm and 8 am including any pauses within a 200-meter radius.
|disttravelled                         |meters      | Total distance travelled over a day (flights).
|rog                                   |meters      | The Radius of Gyration (rog) is a measure in meters of the area covered by a person over a day. A centroid is calculated for all the places (pauses) visited during a day and a weighted distance between all the places and that centroid is computed. The weights are proportional to the time spent in each place.
|maxdiam                               |meters      | The maximum diameter is the largest distance between any two pauses.
|maxhomedist                           |meters      | The maximum distance from home in meters.
|siglocsvisited                        |locations   | The number of significant locations visited during the day. Significant locations are computed using k-means clustering over pauses found in the whole monitoring period. The number of clusters is found iterating k from 1 to 200 stopping until the centroids of two significant locations are within 400 meters of one another.
|avgflightlen                          |meters      | Mean length of all flights.
|stdflightlen                          |meters      | Standard deviation of the length of all flights.
|avgflightdur                          |seconds     | Mean duration of all flights.
|stdflightdur                           |seconds     | The standard deviation of the duration of all flights. 
|probpause                              |     -      | The fraction of a day spent in a pause (as opposed to a flight)
|siglocentropy                          |nats        | Shannon’s entropy measurement based on the proportion of time spent at each significant location visited during a day.
|circdnrtn                              |      -     |   A continuous metric quantifying a person’s circadian routine that can take any value between 0 and 1, where 0 represents a daily routine completely different from any other sensed days and 1 a routine the same as every other sensed day.
|wkenddayrtn                            |       -    | Same as circdnrtn but computed separately for weekends and weekdays.

!!! note "Assumptions/Observations"
    **Barnett\'s et al features**
    These features are based on a Pause-Flight model. A pause is defined as a mobiity trace (location pings) within a certain duration and distance (by default 300 seconds and 60 meters). A flight is any mobility trace between two pauses. Data is resampled and imputed before the features are computed. See [Barnett et al](../../citation#barnett-locations) for more information. In RAPIDS we only expose two parameters for these features (timezone and accuracy limit). You can change other parameters in `src/features/phone_locations/barnett/library/MobilityFeatures.R`.

    **Significant Locations**
    Significant locations are determined using K-means clustering on pauses longer than 10 minutes. The number of clusters (K) is increased until no two clusters are within 400 meters from each other. After this, pauses within a certain range of a cluster (200 meters by default) will count as a visit to that significant location. This description was adapted from the Supplementary Materials of [Barnett et al](../../citation#barnett-locations).

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
    - data/interim/{pid}/phone_locations_features/phone_locations_{language}_{provider_key}.csv
    - data/processed/features/{pid}/phone_locations.csv
    ```


Parameters description for `[PHONE_LOCATIONS][PROVIDERS][DORYAB]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`| Set to `True` to extract `PHONE_LOCATIONS` features from the `BARNETT` provider|
|`[FEATURES]` |         Features to be computed, see table below
| `[DBSCAN_EPS]`             | The maximum distance in meters between two samples for one to be considered as in the neighborhood of the other. This is not a maximum bound on the distances of points within a cluster. This is the most important DBSCAN parameter to choose appropriately for your data set and distance function.
| `[DBSCAN_MINSAMPLES]`      | The number of samples (or total weight) in a neighborhood for a point to be considered as a core point of a cluster. This includes the point itself.
| `[THRESHOLD_STATIC]`       | It is the threshold value in km/hr which labels a row as Static or Moving.
| `[MAXIMUM_GAP_ALLOWED]`   | The maximum gap (in seconds) allowed between any two consecutive rows for them to be considered part of the same displacement. If this threshold is too high, it can throw speed and distance calculations off for periods when the the phone was not sensing.
| `[MINUTES_DATA_USED]`     | Set to `True` to include an extra column in the final location feature file containing the number of minutes used to compute the features on each time segment. Use this for quality control purposes, the more data minutes exist for a period, the more reliable its features should be. For fused location, a single minute can contain more than one coordinate pair if the participant is moving fast enough.
| `[SAMPLING_FREQUENCY]`     | Expected time difference between any two location rows in minutes. If set to `0`, the sampling frequency will be inferred automatically as the median of all the differences between any two consecutive row timestamps (recommended if you are using `FUSED_RESAMPLED` data). This parameter impacts all the time calculations.
| `[CLUSTER_ON]`             | Set this flag to `PARTICIPANT_DATASET` to create clusters based on the entire participant's dataset or to `TIME_SEGMENT` to create clusters based on all the instances of the corresponding time segment (e.g. all mornings).


Features description for `[PHONE_LOCATIONS][PROVIDERS][DORYAB]`:

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
|locationvariance                                            |$meters^2$    |The sum of the variances of the latitude and longitude columns. 
|loglocationvariance                                           | -          | Log of the sum of the variances of the latitude and longitude columns.
|totaldistance                                                |meters        |Total distance travelled in a time segment using the haversine formula.
|averagespeed                                                 |km/hr         |Average speed in a time segment considering only the instances labeled as Moving.
|varspeed                                                      |km/hr         |Speed variance in a time segment considering only the instances labeled as Moving. 
|circadianmovement                                              |-             | \"It encodes the extent to which a person's location patterns follow a 24-hour circadian cycle.\" [Doryab et al.](../../citation#doryab-locations).
|numberofsignificantplaces                                    |places        |Number of significant locations visited. It is calculated using the DBSCAN clustering algorithm which takes in EPS and MIN_SAMPLES as parameters to identify clusters. Each cluster is a significant place.
|numberlocationtransitions                                    |transitions   |Number of movements between any two clusters in a time segment.
|radiusgyration                                               |meters        |Quantifies the area covered by a participant
|timeattop1location                                           |minutes       |Time spent at the most significant location.
|timeattop2location                                           |minutes       |Time spent at the 2nd most significant location.
|timeattop3location                                           |minutes       |Time spent at the 3rd most significant location. 
|movingtostaticratio                                          | -   |  Ratio between the number of rows labeled Moving versus Static
|outlierstimepercent                                          | -   | Ratio between the number of rows that belong to non-significant clusters divided by the total number of rows in a time segment.
|maxlengthstayatclusters                                      |minutes       |Maximum time spent in a cluster (significant location).
|minlengthstayatclusters                                      |minutes       |Minimum time spent in a cluster (significant location).
|meanlengthstayatclusters                                     |minutes       |Average time spent in a cluster (significant location).
|stdlengthstayatclusters                                      |minutes       |Standard deviation of time spent in a cluster (significant location).
|locationentropy                                              |nats          |Shannon Entropy computed over the row count of each cluster (significant location), it will be higher the more rows belong to a cluster (i.e. the more time a participant spent at a significant location).
|normalizedlocationentropy                                    |nats          |Shannon Entropy computed over the row count of each cluster (significant location) divided by the number of clusters, it will be higher the more rows belong to a cluster (i.e. the more time a participant spent at a significant location).


!!! note "Assumptions/Observations"
    **Significant Locations Identified**
    Significant locations are determined using DBSCAN clustering on locations that a patient visit over the course of the period of data collection.

    **The Circadian Calculation**
    For a detailed description of how this is calculated, see [Canzian et al](../../citation#doryab-locations).