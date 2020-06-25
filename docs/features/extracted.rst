.. _rapids_features:

RAPIDS Features
===============

*How do I compute any of these features?* In your ``config.yaml``, go to the sensor section you are interested in and set the corresponding ``COMPUTE`` option to ``TRUE`` as well as ``DB_TABLE`` to the senor's table name in your database (the default table name is the one assigned by Aware), for example
::

    MESSAGES:
        COMPUTE: True
        DB_TABLE: messages
        ...
 
 
 If you want to extract phone_valid_sensed_days.csv, screen features or locaton features based on fused location data don't forget to configure ``TABLES_FOR_SENSED_BINS`` (see below).

.. _global-sensor-doc:

Global Parameters
"""""""""""""""""

.. _sensor-list:

- ``TABLES_FOR_SENSED_BINS`` - Add as many sensor tables as you have in your database. All sensors included are used to compute ``phone_sensed_bins.csv`` (bins of time when the smartphone was sensing data). In turn, these bins are used to compute ``PHONE_VALID_SENSED_DAYS`` (see below), ``episodepersensedminutes`` feature of :ref:`Screen<screen-sensor-doc>` and to resample fused location data if you configure Barnett's location features to use ``RESAMPLE_FUSED``. See TABLES_FOR_SENSED_BINS_ variable in ``config`` file (therefore, when you are extracting screen or Barnett's location features, screen and locations tables are mandatory).  

.. _pid: 

- ``PID`` - The list of participant ids to be included in the analysis. These should match the names of the files created in the ``data/external`` directory  (:ref:`see more details<db-configuration>`).

.. _day-segments: 

- ``DAY_SEGMENTS`` - The list of day epochs that features can be segmented into: ``daily``, ``morning`` (6am-12pm), ``afternnon`` (12pm-6pm), ``evening`` (6pm-12am) and ``night`` (12am-6am). This list can be modified globally or on a per sensor basis. See DAY_SEGMENTS_ in ``config`` file.

.. _timezone:

- ``TIMEZONE`` - The time zone where data was collected. Use the timezone names from this `List of Timezones`_. Double check your chosen name is correct, for example US Eastern Time is called New America/New_York, not EST.

.. _database_group:

- ``DATABASE_GROUP`` - The name of your database credentials group, it should match the one in ``.env`` (:ref:`see the datbase configuration<db-configuration>`). 

.. _download-dataset:

- ``DOWNLOAD_DATASET``

    - ``GROUP``. Credentials group to connect to the database containing ``SENSORS``. By default it points to ``DATABASE_GROUP``.

.. _readable-datetime:

- ``READABLE_DATETIME`` - Configuration to convert UNIX timestamps into readbale date time strings.

    - ``FIXED_TIMEZONE``. See ``TIMEZONE`` above. This assumes that all data of all participants was collected within one time zone.
    - Support for multiple time zones for each participant coming soon based on the ``timezone`` table collected by Aware.

.. _phone-valid-sensed-days:

- ``PHONE_VALID_SENSED_DAYS``.
    
    Contains three attributes: ``BIN_SIZE``, ``MIN_VALID_HOURS``, ``MIN_BINS_PER_HOUR``. 

    On any given day, Aware could have sensed data only for a few minutes or for 24 hours. Daily estimates of features should be considered more reliable the more hours Aware was running and logging data (for example, 10 calls logged on a day when only one hour of data was recorded is a less reliable feature compared to 10 calls on a day when 23 hours of data were recorded. 

    Therefore, we define a valid hour as those that contain a minimum number of valid bins. In turn, a valid bin are those that contain at least one row of data from any sensor logged within that period. We divide an hour into N bins of size ``BIN_SIZE`` (in minutes) and we mark an hour as valid if contains at least ``MIN_BINS_PER_HOUR`` (out of the total possible number of bins that can be captured in an hour based on their length i.e. 60min/``BIN_SIZE`` bins). Days with valid sensed hours less than ``MIN_VALID_HOURS`` will be excluded form the output of this file. See PHONE_VALID_SENSED_DAYS_ in ``config.yaml``.

    Note that RAPIDS *DOES NOT* filter your feature files automatically, you need to do this manually based on ``"data/interim/{pid}/phone_valid_sensed_days.csv"``.

    You can get access to every phone's sensed bins matrix (days x bins) in ``data/interim/{pid}/phone_sensed_bins.csv``. As mentioned above, RAPIDS uses this file to compute ``phone_valid_sensed_days.csv``, ``episodepersensedminutes`` feature of :ref:`Screen<screen-sensor-doc>` and to resample fused location data if you configure Barnett's location features to use ``RESAMPLE_FUSED``.

.. _individual-sensor-settings:


.. _messages-sensor-doc:

Messages (SMS)
"""""""""""""""

See `Messages Config Code`_

**Available Epochs (day_segment) :** daily, morning, afternoon, evening, night

**Available Platforms:** Android

**Rule Chain:**

- Rule ``rules/preprocessing.snakefile/download_dataset``
- Rule ``rules/preprocessing.snakefile/readable_datetime``
- Rule ``rules/features.snakefile/messages_features``

.. _messages-parameters:

**Messages Rule Parameters (messages_features):**

==============    ===================
Name	          Description
==============    ===================
messages_type     The particular ``messages_type`` that will be analyzed. The options for this parameter are ``received`` or ``sent``.
day_segment       The particular ``day_segment`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
features          Features to be computed, see table below
==============    ===================

.. _messages-available-features:

**Available Message Features**

=========================   =========     =============
Name                        Units         Description
=========================   =========     =============
count                       messages      Number of messages of type ``messages_type`` that occurred during a particular ``day_segment``.
distinctcontacts            contacts      Number of distinct contacts that are associated with a particular ``messages_type`` during a particular ``day_segment``.
timefirstsms                minutes       Number of minutes between 12:00am (midnight) and the first ``message`` of a particular ``messages_type``.
timelastsms                 minutes       Number of minutes between 12:00am (midnight) and the last ``message`` of a particular ``messages_type``.
countmostfrequentcontact    messages      Number of messages from the contact with the most messages of ``messages_type`` during a ``day_segment`` throughout the whole dataset of each participant.
=========================   =========     =============

**Assumptions/Observations:** 

``TYPES`` and ``FEATURES`` keys in ``config.yaml`` need to match. For example, below the ``TYPE`` ``sent`` matches the ``FEATURES`` key ``sent``::

        MESSAGES:
            ...
            TYPES: [sent]
            FEATURES: 
                sent: [count, distinctcontacts, timefirstsms, timelastsms, countmostfrequentcontact]


.. _call-sensor-doc:

Calls
""""""

See `Call Config Code`_

**Available Epochs (day_segment) :** daily, morning, afternoon, evening, night

**Available Platforms:** Android and iOS

**Rule Chain:**

- Rule ``rules/preprocessing.snakefile/download_dataset``
- Rule ``rules/preprocessing.snakefile/readable_datetime``
- Rule ``rules/features.snakefile/call_features``
    
.. _calls-parameters:

**Call Rule Parameters (call_features):**

============    ===================
Name	        Description
============    ===================
call_type       The particular ``call_type`` that will be analyzed. The options for this parameter are ``incoming``, ``outgoing`` or ``missed``.
day_segment     The particular ``day_segment`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
features        Features to be computed. Note that the same features are available for both ``incoming`` and ``outgoing`` calls, while ``missed`` calls has its own set of features. See :ref:`Available Incoming and Outgoing Call Features <available-in-and-out-call-features>` Table and :ref:`Available Missed Call Features <available-missed-call-features>` Table below.
============    ===================

.. _available-in-and-out-call-features:

**Available Incoming and Outgoing Call Features**

=========================   =========     =============
Name                        Units         Description
=========================   =========     =============
count                       calls         Number of calls of a particular ``call_type`` occurred during a particular ``day_segment``.
distinctcontacts            contacts      Number of distinct contacts that are associated with a particular ``call_type`` for a particular ``day_segment``
meanduration                seconds       The mean duration of all calls of a particular ``call_type`` during a particular ``day_segment``.
sumduration                 seconds       The sum of the duration of all calls of a particular ``call_type`` during a particular ``day_segment``.
minduration                 seconds       The duration of the shortest call of a particular ``call_type`` during a particular ``day_segment``.
maxduration                 seconds       The duration of the longest call of a particular ``call_type`` during a particular ``day_segment``.
stdduration                 seconds       The standard deviation of the duration of all the calls of a particular ``call_type`` during a particular ``day_segment``.
modeduration                seconds       The mode of the duration of all the calls of a particular ``call_type`` during a particular ``day_segment``.
entropyduration             nats          The estimate of the Shannon entropy for the the duration of all the calls of a particular ``call_type`` during a particular ``day_segment``.
timefirstcall               minutes       The time in minutes between 12:00am (midnight) and the first call of ``call_type``.
timelastcall                minutes       The time in minutes between 12:00am (midnight) and the last call of ``call_type``.
countmostfrequentcontact    calls         The number of calls of a particular ``call_type`` during a particular ``day_segment`` of the most frequent contact throughout the monitored period.
=========================   =========     =============

.. _available-missed-call-features:

**Available Missed Call Features**

=========================   =========     =============
Name                        Units         Description
=========================   =========     =============
count                       calls         Number of ``missed`` calls that occurred during a particular ``day_segment``.
distinctcontacts            contacts      Number of distinct contacts that are associated with ``missed`` calls for a particular ``day_segment``
timefirstcall               minutes       The time in hours from 12:00am (Midnight) that the first ``missed`` call occurred.
timelastcall                minutes       The time in hours from 12:00am (Midnight) that the last ``missed`` call occurred.
countmostfrequentcontact    calls         The number of ``missed`` calls during a particular ``day_segment`` of the most frequent contact throughout the monitored period.
=========================   =========     =============

**Assumptions/Observations:** 

``TYPES`` and ``FEATURES`` keys in ``config.yaml`` need to match. For example, below the ``TYPE`` ``missed`` matches the ``FEATURES`` key ``missed``::

    CALLS:
        ...
        TYPES: [missed]
        FEATURES: 
            missed: [count, distinctcontacts, timefirstcall, timelastcall, countmostfrequentcontact]

Aware Android client stores call types 1=incoming, 2=outgoing, 3=missed while Aware iOS client stores call status 1=incoming, 2=connected, 3=dialing, 4=disconnected. We extract iOS call types based on call status sequences: (1,2,4)=incoming=1, (3,2,4)=outgoing=2, (1,4) or (3,4)=missed=3. Sometimes (due to a possible bug in Aware) sequences get logged on the exact same timestamp, thus 3-item sequences can be 2,3,4 or 3,2,4. Although iOS stores the duration of ringing/dialing stages for missed calls, we set it to 0 to match Android.


.. _bluetooth-sensor-doc:

Bluetooth
""""""""""

See `Bluetooth Config Code`_

**Available Epochs (day_segment) :** daily, morning, afternoon, evening, night

**Available Platforms:** Android and iOS

**Snakemake rule chain:**

- Rule ``rules/preprocessing.snakefile/download_dataset``
- Rule ``rules/preprocessing.snakefile/readable_datetime``
- Rule ``rules/features.snakefile/bluetooth_features``
    
.. _bluetooth-parameters:

**Bluetooth Rule Parameters (bluetooth_features):**

============    ===================
Name	        Description
============    ===================
day_segment     The particular ``day_segment`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
features        Features to be computed, see table below
============    ===================

.. _bluetooth-available-features:

**Available Bluetooth Features**

===========================   =========     =============
Name                          Units         Description
===========================   =========     =============
countscans                    devices       Number of scanned devices during a ``day_segment``, a device can be detected multiple times over time and these appearances are counted separately
uniquedevices                 devices       Number of unique devices during a ``day_segment`` as identified by their hardware address
countscansmostuniquedevice    scans         Number of scans of the most scanned device during a ``day_segment`` across the whole monitoring period
===========================   =========     =============

**Assumptions/Observations:** N/A 


.. _wifi-sensor-doc:

WiFi
""""""""""

See `WiFi Config Code`_

**Available Epochs (day_segment) :** daily, morning, afternoon, evening, night

**Available Platforms:** Android and iOS

**Snakemake rule chain:**

- Rule ``rules/preprocessing.snakefile/download_dataset``
- Rule ``rules/preprocessing.snakefile/readable_datetime``
- Rule ``rules/features.snakefile/wifi_features``
    
.. _wifi-parameters:

**WiFi Rule Parameters (wifi_features):**

============    ===================
Name	        Description
============    ===================
day_segment     The particular ``day_segment`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
features        Features to be computed, see table below
============    ===================

.. _wifi-available-features:

**Available WiFi Features**

===========================   =========     =============
Name                          Units         Description
===========================   =========     =============
countscans                    devices       Number of scanned WiFi access points during a ``day_segment``, an access point can be detected multiple times over time and these appearances are counted separately
uniquedevices                 devices       Number of unique access point during a ``day_segment`` as identified by their hardware address
countscansmostuniquedevice    scans         Number of scans of the most scanned access point during a ``day_segment`` across the whole monitoring period
===========================   =========     =============

**Assumptions/Observations:** N/A 


.. _accelerometer-sensor-doc:

Accelerometer
""""""""""""""

See `Accelerometer Config Code`_

**Available Epochs (day_segment) :** daily, morning, afternoon, evening, night

**Available Platforms:** Android and iOS

**Rule chain:**

- Rule ``rules/preprocessing.snakefile/download_dataset``
- Rule ``rules/preprocessing.snakefile/readable_datetime``
- Rule ``rules/features.snakefile/accelerometer_features``
    
.. _Accelerometer-parameters:

**Accelerometer Rule Parameters (accelerometer_features):**

============    ===================
Name	        Description
============    ===================
day_segment     The particular ``day_segment`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
features        Features to be computed, see table below
============    ===================

.. _accelerometer-available-features:

**Available Accelerometer Features**

======================    ==============    =============
Name                      Units             Description
======================    ==============    =============
maxmagnitude              m/s\ :sup:`2`     The maximum magnitude of acceleration (:math:`\|acceleration\| = \sqrt{x^2 + y^2 + z^2}`).
minmagnitude              m/s\ :sup:`2`     The minimum magnitude of acceleration.
avgmagnitude              m/s\ :sup:`2`     The average magnitude of acceleration.
medianmagnitude           m/s\ :sup:`2`     The median magnitude of acceleration.
stdmagnitude              m/s\ :sup:`2`     The standard deviation of acceleration.
sumduration               minutes           Total duration of all exertional or non-exertional activity episodes.
maxduration               minutes           Longest duration of any exertional or non-exertional activity episode.
minduration               minutes           Shortest duration of any exertional or non-exertional activity episode.
avgduration               minutes           Average duration of any exertional or non-exertional activity episode.
medianduration            minutes           Median duration of any exertional or non-exertional activity episode.
stdduration               minutes           Standard deviation of the duration of all exertional or non-exertional activity episodes.
======================    ==============    =============

**Assumptions/Observations:**

Exertional activity episodes are based on this paper: Panda N, Solsky I, Huang EJ, et al. Using Smartphones to Capture Novel Recovery Metrics After Cancer Surgery. JAMA Surg. 2020;155(2):123–129. doi:10.1001/jamasurg.2019.4702


.. _applications-foreground-sensor-doc:

Applications Foreground
""""""""""""""""""""""""

See `Applications Foreground Config Code`_

**Available Epochs (day_segment) :** daily, morning, afternoon, evening, night

**Available Platforms:** Android

**Snakemake rule chain:**

- Rule ``rules/preprocessing.snakefile/download_dataset`` 
- Rule ``rules/preprocessing.snakefile/readable_datetime`` 
- Rule ``rules/preprocessing.snakefile/application_genres``
- Rule ``rules/features.snakefile/applications_foreground_features`` 
   
.. _applications-foreground-parameters:

**Applications Foreground Rule Parameters (applications_foreground_features):**

====================    ===================
Name	                Description
====================    ===================
day_segment             The particular ``day_segment`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
single_categories       App categories to be included in the feature extraction computation. See ``APPLICATION_GENRES`` in this file to add new categories or use the catalogue we provide and read :ref:`Assumtions and Observations <applications-foreground-observations>` for more information.
multiple_categories     You can group multiple categories into meta categories, for example ``social: ["socialnetworks", "socialmediatools"]``.
single_apps             Apps to be included in the feature extraction computation. Use their package name, for example, ``com.google.android.youtube`` or the reserved word ``top1global`` (the most used app by a participant over the whole monitoring study).
excluded_categories     App categories to be excluded in the feature extraction computation. See ``APPLICATION_GENRES`` in this file to add new categories or use the catalogue we provide and read :ref:`Assumtions and Observations <applications-foreground-observations>` for more information.
excluded_apps           Apps to be excluded in the feature extraction computation. Use their package name, for example: ``com.google.android.youtube``
features                Features to be computed, see table below
====================    ===================

.. _applications-foreground-available-features:

**Available Applications Foreground Features**

==================   =========   =============
Name                 Units       Description
==================   =========   =============
count                apps        Number of times a single app or apps within a category were used (i.e. they were brought to the foreground either by tapping their icon or switching to it from another app).
timeoffirstuse       minutes     The time in minutes between 12:00am (midnight) and the first use of a single app or apps within a category during a ``day_segment``.
timeoflastuse        minutes     The time in minutes between 12:00am (midnight) and the last use of a single app or apps within a category during a ``day_segment``.
frequencyentropy     nats        The entropy of the used apps within a category during a ``day_segment`` (each app is seen as a unique event, the more apps were used, the higher the entropy). This is especially relevant when computed over all apps. Entropy cannot be obtained for a single app.
==================   =========   =============

.. _applications-foreground-observations:

**Assumptions/Observations:** 

Features can be computed by app, by apps grouped under a single category (genre) and by multiple categories grouped together (meta categories). For example, we can get features for Facebook, for Social Network Apps (including Facebook and others) or for a meta category called Social formed by Social Network and Social Media Tools categories. 

Apps installed by default like YouTube are considered systems apps on some phones. We do an exact match to exclude apps where "genre" == ``EXCLUDED_CATEGORIES`` or "package_name" == ``EXCLUDED_APPS``.

We provide three ways of classifying and app within a category (genre): a) by automatically scraping its official category from the Google Play Store, b) by using the catalogue created by Stachl et al. which we provide in RAPIDS (``data/external/``), or c) by manually creating a personalized catalogue.

The way you choose strategy a, b or c is by modifying ``APPLICATION_GENRES`` keys and values. Set ``CATALOGUE_SOURCE`` to ``FILE`` if you want to use a CSV file as catalogue (strategy b and c) or to ``GOOGLE`` if you want to scrape the genres from the Play Store (strategy a). By default ``CATALOGUE_FILE`` points to the catalogue created by  Stachl et al. (strategy b) and you can change this path to your own catalogue that follows the same format (strategy c). In addition, set ``SCRAPE_MISSING_GENRES`` to true if you are using a FILE catalogue and you want to scrape from the Play Store any missing genres and ``UPDATE_CATALOGUE_FILE`` to true if you want to save those scrapped genres back into the FILE.

The genre catalogue we provide was shared as part of the Supplemental Materials of Stachl, C., Au, Q., Schoedel, R., Buschek, D., Völkel, S., Schuwerk, T., … Bühner, M. (2019, June 12). Behavioral Patterns in Smartphone Usage Predict Big Five Personality Traits. https://doi.org/10.31234/osf.io/ks4vd 

.. _battery-sensor-doc:

Battery
"""""""""

See `Battery Config Code`_

**Available Epochs (day_segment) :** daily, morning, afternoon, evening, night

**Available Platforms:** Android and iOS

**Snakemake rule chain:**

- Rule ``rules/preprocessing.snakefile/download_dataset`` 
- Rule ``rules/preprocessing.snakefile/readable_datetime`` 
- Rule ``rules/features.snakefile/battery_deltas`` 
- Rule ``rules/features.snakefile/battery_features``
    
.. _battery-parameters:

**Battery Rule Parameters (battery_features):**

============    ===================
Name	        Description
============    ===================
day_segment     The particular ``day_segment`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
features        Features to be computed, see table below
============    ===================

.. _battery-available-features:

**Available Battery Features**

=====================   =================   =============
Name                    Units               Description
=====================   =================   =============
countdischarge          episodes            Number of discharging episodes.
sumdurationdischarge    minutes             The total duration of all discharging episodes.
countcharge             episodes            Number of battery charging episodes.
sumdurationcharge       minutes             The total duration of all charging episodes.
avgconsumptionrate      episodes/minutes    The average of all episodes’ consumption rates. An episode’s consumption rate is defined as the ratio between its battery delta and duration
maxconsumptionrate      episodes/minutes    The highest of all episodes’ consumption rates. An episode’s consumption rate is defined as the ratio between its battery delta and duration
=====================   =================   =============

**Assumptions/Observations:** 

For Aware iOS client V1 we swap battery status 3 to 5 and 1 to 3, client V2 does not have this problem.

.. _activity-recognition-sensor-doc:


Activity Recognition
""""""""""""""""""""""""""""

See `Activity Recognition Config Code`_

**Available Epochs:** daily, morning, afternoon, evening, night

**Available Platforms:** Android and iOS

**Snakemake rule chain:**

- Rule ``rules/preprocessing.snakefile/download_dataset`` 
- Rule ``rules/preprocessing.snakefile/readable_datetime`` 
- Rule ``rules/preprocessing.snakefile/unify_ios_android`` 
- Rule ``rules/features.snakefile/google_activity_recognition_deltas``
- Rule ``rules/features.snakefile/ios_activity_recognition_deltas``
- Rule ``rules/features.snakefile/activity_features``
    
.. _activity-recognition-parameters:

**Rule Parameters (activity_features):**

============    ===================
Name	        Description
============    ===================
day_segment     The particular ``day_segment`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
features        Features to be computed, see table below
============    ===================

.. _activity-recognition-available-features:

**Available Activity Recognition Features**

======================   ============    =============
Name                     Units           Description
======================   ============    =============
count                    rows            Number of detect activity events (rows).
mostcommonactivity       factor          The most common activity.
countuniqueactivities    activities      Number of unique activities.
activitychangecount      transitions     Number of transitions between two different activities; still to running for example.
sumstationary            minutes         The total duration of episodes of still and tilting (phone) activities.
summobile                minutes         The total duration of episodes of on foot, running, and on bicycle activities
sumvehicle               minutes         The total duration of episodes of on vehicle activity
======================   ============    =============

**Assumptions/Observations:**

iOS Activity Recognition data labels are unified with Google Activity Recognition labels: "automotive" to "in_vehicle", "cycling" to "on_bicycle", "walking" and "running" to "on_foot", "stationary" to "still". In addition, iOS activity pairs formed by "stationary" and "automotive" labels (driving but stopped at a traffic light) are transformed to "automotive" only.

In AWARE, Activity Recognition data for Google (Android) and iOS are stored in two different database tables, RAPIDS (via Snakemake) automatically infers what platform each participant belongs to based on their participant file (``data/external/``) which in turn takes this information from the ``aware_device`` table (see ``optional_ar_input`` function in ``rules/features.snakefile``). 

.. _light-doc:

Light
"""""""

See `Light Config Code`_

**Available Epochs (day_segment) :** daily, morning, afternoon, evening, night

**Available Platforms:** Android

**Rule Chain:**

- Rule: ``rules/preprocessing.snakefile/download_dataset``
- Rule: ``rules/preprocessing.snakefile/readable_datetime``
- Rule: ``rules/features.snakefile/light_features``

.. _light-parameters:

**Light Rule Parameters (light_features):**

============    ===================
Name	        Description
============    ===================
day_segment     The particular ``day_segment`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
features        Features to be computed, see table below
============    ===================

.. _light-available-features:

**Available Light Features**

===========   =========     =============
Name          Units         Description
===========   =========     =============
count         rows          Number light sensor rows recorded.
maxlux        lux           The maximum ambient luminance.
minlux        lux           The minimum ambient luminance.
avglux        lux           The average ambient luminance.
medianlux     lux           The median ambient luminance.
stdlux        lux           The standard deviation of ambient luminance.
===========   =========     =============

**Assumptions/Observations:** N/A


.. _location-sensor-doc:

Location (Barnett’s) Features
""""""""""""""""""""""""""""""
Barnett’s location features are based on the concept of flights and pauses. GPS coordinates are converted into a 
sequence of flights (straight line movements) and pauses (time spent stationary). Data is imputed before features 
are computed. See Ian Barnett, Jukka-Pekka Onnela, Inferring mobility measures from GPS traces with missing data, Biostatistics, Volume 21, Issue 2, April 2020, Pages e98–e112, https://doi.org/10.1093/biostatistics/kxy059. The code for these features was made open source by Ian Barnett (https://scholar.harvard.edu/ibarnett/software/gpsmobility).

See `Location (Barnett’s) Config Code`_

**Available Epochs (day_segment) :** daily

**Available Platforms:** Android and iOS

**Snakemake rule chain:**

- Rule ``rules/preprocessing.snakefile/download_dataset``
- Rule ``rules/preprocessing.snakefile/readable_datetime``
- Rule ``rules/preprocessing.snakefile/phone_sensed_bins``
- Rule ``rules/preprocessing.snakefile/resample_fused_location`` (only relevant if setting ``location_to_use`` to ````RESAMPLE_FUSED``.
- Rule ``rules/features.snakefile/location_barnett_features``
    
.. _location-parameters:

**Location Rule Parameters (location_barnett_features):**

=================    ===================
Name	             Description
=================    ===================
location_to_use      *Read the Observations section below*. The specifies what type of location data will be use in the analysis. Possible options are ``ALL``, ``ALL_EXCEPT_FUSED`` OR ``RESAMPLE_FUSED``
accuracy_limit       This is in meters. The sensor drops location coordinates with an accuracy higher than this. This number means there's a 68% probability the true location is within this radius specified.
timezone             The timezone used to calculate location.
minutes_data_used    This is NOT a feature. This is just a quality control check, and if set to TRUE, a new column is added to the output file with the number of minutes containing location data that were used to compute all features. The more data minutes exist for a period, the more reliable its features should be. For fused location, a single minute can contain more than one coordinate pair if the participant is moving fast enough.
features             Features to be computed, see table below
=================    ===================

.. _location-available-features:

**Available Location Features**

Description taken from `Beiwe Summary Statistics`_.

================   =========     =============
Name               Units         Description
================   =========     =============
hometime           minutes       Time at home. Time spent at home in minutes. Home is the most visited significant location between 8 pm and 8 am including any pauses within a 200-meter radius.
disttravelled      meters        Total distance travelled over a day (flights).
rog                meters        The Radius of Gyration (rog) is a measure in meters of the area covered by a person over a day. A centroid is calculated for all the places (pauses) visited during a day and a weighted distance between all the places and that centroid is computed. The weights are proportional to the time spent in each place.
maxdiam            meters        The maximum diameter is the largest distance between any two pauses.
maxhomedist        meters        The maximum distance from home in meters.
siglocsvisited     locations     The number of significant locations visited during the day. Significant locations are computed using k-means clustering over pauses found in the whole monitoring period. The number of clusters is found iterating k from 1 to 200 stopping until the centroids of two significant locations are within 400 meters of one another.
avgflightlen       meters        Mean length of all flights.
stdflightlen       meters        Standard deviation of the length of all flights.
avgflightdur       seconds       Mean duration of all flights.
stdflightdur       seconds       The standard deviation of the duration of all flights.
probpause                        The fraction of a day spent in a pause (as opposed to a flight)
siglocentropy      nats          Shannon’s entropy measurement based on the proportion of time spent at each significant location visited during a day.
circdnrtn           	         A continuous metric quantifying a person’s circadian routine that can take any value between 0 and 1, where 0 represents a daily routine completely different from any other sensed days and 1 a routine the same as every other sensed day.
wkenddayrtn                      Same as circdnrtn but computed separately for weekends and weekdays.
================   =========     =============

**Assumptions/Observations:** 

*Types of location data to use*

Aware Android and iOS clients can collect location coordinates through the phone's GPS or Google's fused location API. If your Aware client was ONLY configured to use GPS set ``location_to_use`` to ``ALL``, if your client was configured to use BOTH GPS and fused location you can use ``ALL`` or set ``location_to_use`` to  ``ALL_EXCEPT_FUSED`` to ignore fused coordinates, if your client was configured to use fused location only,  set ``location_to_use`` to ``RESAMPLE_FUSED``. ``RESAMPLE_FUSED`` takes the original fused location coordinates and replicates each pair forward in time as long as the phone was sensing data as indicated by ``phone_sensed_bins`` (see :ref:`Phone valid sensed days <phone-valid-sensed-days>`), this is done because Google's API only logs a new location coordinate pair when it is sufficiently different from the previous one. 

There are two parameters associated with resampling fused location in the ``RESAMPLE_FUSED_LOCATION`` section of the ``config.yaml`` file. ``CONSECUTIVE_THRESHOLD`` (in minutes, default 30) controls the maximum gap between any two coordinate pairs to replicate the last known pair (for example, participant A's phone did not collect data between 10.30am and 10:50am and between 11:05am and 11:40am, the last known coordinate pair will be replicated during the first period but not the second, in other words, we assume that we cannot longer guarantee the participant stayed at the last known location if the phone did not sense data for more than 30 minutes). ``TIME_SINCE_VALID_LOCATION`` (in minutes, default 720 or 12 hours) the last known fused location won't be carried over longer that this threshold even if the phone was sensing data continuously (for example, participant A went home at 9pm and their phone was sensing data without gaps until 11am the next morning, the last known location will only be replicated until 9am). If you have suggestions to modify or improve this imputation, let us know.

*Significant Locations Identified*

(i.e. The clustering method used)
Significant locations are determined using K-means clustering on locations that a patient visit over the course of the period of data collection. By setting K=K+1 and repeat clustering until two significant locations are within 100 meters of one another, the results from the previous step (K-1) can   be used as the total number of significant locations. Taken from `Beiwe Summary Statistics`_. 

*Definition of Stationarity*

(i.e., The length of time and distance a person has to be around the same place to be labelled as a pause)
This is based on a Pause-Flight model, The parameters used are a minimum pause duration of 300sec and a minimum pause distance of 60m. See the `Pause-Flight Model`_.

*The Circadian Calculation*

For a detailed description of how this is calculated, see Canzian, L., & Musolesi, M. (2015, September). Trajectories of depression: unobtrusive monitoring of depressive states by means of smartphone mobility traces analysis. In Proceedings of the 2015 ACM international joint conference on pervasive and ubiquitous computing (pp. 1293-1304). Their procedure was followed using 30-min increments as a bin size. Taken from `Beiwe Summary Statistics`_.

.. _screen-sensor-doc:

Screen
""""""""

See `Screen Config Code`_

**Available Epochs (day_segment) :** daily, morning, afternoon, evening, night

**Available Platforms:** Android and iOS

**Snakemake rule chain:**

- Rule ``rules/preprocessing.snakefile/download_dataset``
- Rule ``rules/preprocessing.snakefile/readable_datetime``
- Rule ``rules/features.snakefile/screen_deltas``
- Rule ``rules/features.snakefile/screen_features``

.. _screen-parameters:

**Screen Rule Parameters (screen_features):**

=========================    ===================
Name	                     Description
=========================    ===================
day_segment                  The particular ``day_segments`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
reference_hour_first_use     The reference point from which ``firstuseafter`` is to be computed, default is midnight
features_deltas              Features to be computed, see table below
episode_types                Currently we only support unlock episodes (from when the phone is unlocked until the screen is off)
=========================    ===================

.. _screen-episodes-available-features:

**Available Screen Episodes Features**

=========================   =================   =============
Name                        Units               Description
=========================   =================   =============
sumduration                 minutes             Total duration of all unlock episodes.
maxduration                 minutes             Longest duration of any unlock episode.
minduration                 minutes             Shortest duration of any unlock episode.
avgduration                 minutes             Average duration of all unlock episodes.
stdduration                 minutes             Standard deviation duration of all unlock episodes.
countepisode                episodes            Number of all unlock episodes
episodepersensedminutes     episodes/minute     The ratio between the total number of episodes in an epoch divided by the total time (minutes) the phone was sensing data.
firstuseafter               minutes             Seconds until the first unlock episode.
=========================   =================   =============

**Assumptions/Observations:** 

An ``unlock`` episode is considered as the time between an ``unlock`` event and a ``lock`` event. iOS recorded these episodes reliably (albeit some duplicated ``lock`` events within milliseconds from each other). However, in Android there are some events unrelated to the screen state because of multiple consecutive ``unlock``/``lock`` events, so we keep the closest pair. In our experiments these cases are less than 10% of the screen events collected. This happens because ``ACTION_SCREEN_OFF`` and ``ON`` are "sent when the device becomes non-interactive which may have nothing to do with the screen turning off". Additionally, in Android it is possible to measure the time spent on the ``lock`` screen before an ``unlock`` event as well as the total screen time (i.e. ``ON`` to ``OFF``) but we are only keeping ``unlock`` episodes (``unlock`` to ``OFF``) to be consistent with iOS. 

.. _conversation-sensor-doc:

Conversation
""""""""""""""

See `Conversation Config Code`_

**Available Epochs (day_segment) :** daily, morning, afternoon, evening, night

**Available Platforms:** Android and iOS

**Snakemake rule chain:**

- Rule ``rules/preprocessing.snakefile/download_dataset``
- Rule ``rules/preprocessing.snakefile/readable_datetime``
- Rule ``rules/features.snakefile/conversation_features``

.. _conversation-parameters:

**Conversation Rule Parameters (conversation_features):**

=========================    ===================
Name	                     Description
=========================    ===================
day_segment                  The particular ``day_segments`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
recordingMinutes             Minutes the plugin was recording audio (default 1 min)
pausedMinutes                Minutes the plugin was NOT recording audio (default 3 min)
features                     Features to be computed, see table below
=========================    ===================

.. _conversation-available-features:

**Available Conversation Features**

=========================   =================   =============
Name                        Units               Description
=========================   =================   =============
minutessilence              minutes             Minutes labeled as silence
minutesnoise                minutes             Minutes labeled as noise
minutesvoice                minutes             Minutes labeled as voice
minutesunknown              minutes             Minutes labeled as unknown
sumconversationduration     minutes             Total duration of all conversations
maxconversationduration     minutes             Longest duration of all conversations
minconversationduration     minutes             Shortest duration of all conversations
avgconversationduration     minutes             Average duration of all conversations
sdconversationduration      minutes             Standard Deviation of the duration of all conversations
timefirstconversation       minutes             Minutes since midnight when the first conversation for a day segment was detected
timelastconversation        minutes             Minutes since midnight when the last conversation for a day segment was detected
sumenergy                   L2-norm             Sum of all energy values
avgenergy                   L2-norm             Average of all energy values
sdenergy                    L2-norm             Standard Deviation of all energy values
minenergy                   L2-norm             Minimum of all energy values
maxenergy                   L2-norm             Maximum of all energy values
silencesensedfraction                           Ratio between minutessilence and the sum of (minutessilence, minutesnoise, minutesvoice, minutesunknown)
noisesensedfraction                             Ratio between minutesnoise and the sum of (minutessilence, minutesnoise, minutesvoice, minutesunknown)
voicesensedfraction                             Ratio between minutesvoice and the sum of (minutessilence, minutesnoise, minutesvoice, minutesunknown)
unknownsensedfraction                           Ratio between minutesunknown and the sum of (minutessilence, minutesnoise, minutesvoice, minutesunknown)
silenceexpectedfraction                         Ration between minutessilence and the number of minutes that in theory should have been sensed based on the record and pause cycle of the plugin (1440 / recordingMinutes+pausedMinutes)
noiseexpectedfraction                           Ration between minutesnoise and the number of minutes that in theory should have been sensed based on the record and pause cycle of the plugin (1440 / recordingMinutes+pausedMinutes)
voiceexpectedfraction                           Ration between minutesvoice and the number of minutes that in theory should have been sensed based on the record and pause cycle of the plugin (1440 / recordingMinutes+pausedMinutes)
unknownexpectedfraction                         Ration between minutesunknown and the number of minutes that in theory should have been sensed based on the record and pause cycle of the plugin (1440 / recordingMinutes+pausedMinutes)
=========================   =================   =============

**Assumptions/Observations:** 


.. ------------------------------- Begin Fitbit Section ----------------------------------- ..

.. _fitbit-sleep-sensor-doc:

Fitbit: Sleep
"""""""""""""""""""

See `Fitbit: Sleep Config Code`_

**Available Epochs (day_segment) :** daily

**Available Platforms:**: Fitbit
    
**Snakemake rule chain:**

- Rule ``rules/preprocessing.snakefile/download_dataset``
- Rule ``rules/preprocessing.snakefile/fitbit_with_datetime``
- Rule ``rules/features.snakefile/fitbit_sleep_features``
    
.. _fitbit-sleep-parameters:

**Fitbit: Sleep Rule Parameters (fitbit_sleep_features):**

==================================    ===================
Name	                              Description
==================================    ===================
day_segment                           The particular ``day_segment`` that will be analyzed. For this sensor only ``daily`` is used.
sleep_types                           The types of sleep provided by Fitbit: ``main``, ``nap``, ``all``.
daily_features_from_summary_data      The sleep features that can be computed based on Fitbit's summary data. See :ref:`Available Fitbit: Sleep Features <fitbit-sleep-available-features>` Table below
==================================    ===================

.. _fitbit-sleep-available-features:

**Available Fitbit: Sleep Features**

========================   ===========    =============
Name                       Units          Description
========================   ===========    =============
sumdurationtofallasleep    minutes        Time it took the user to fall asleep for ``sleep_type`` during ``day_segment``.
sumdurationawake           minutes        Time the user was awake but still in bed for ``sleep_type`` during ``day_segment``.
sumdurationasleep          minutes        Sleep duration for ``sleep_type`` during ``day_segment``.
sumdurationafterwakeup     minutes        Time the user stayed in bed after waking up for ``sleep_type`` during ``day_segment``.
sumdurationinbed           minutes        Total time the user stayed in bed (sumdurationtofallasleep + sumdurationawake + sumdurationasleep + sumdurationafterwakeup) for ``sleep_type`` during ``day_segment``.
avgefficiency              scores         Sleep efficiency average for ``sleep_type`` during ``day_segment``.
countepisode               episodes       Number of sleep episodes for ``sleep_type`` during ``day_segment``.
========================   ===========    =============

**Assumptions/Observations:** 

Only features from summary data are available at the momement.

The `fitbit_with_datetime` rule will extract Summary data (`fitbit_sleep_summary_with_datetime.csv`) and Intraday data (`fitbit_sleep_intraday_with_datetime.csv`). There are two versions of Fitbit's sleep API (`version 1`_ and `version 1.2`_), and each provides raw sleep data in a different format:
    
    - Sleep level. In ``v1``, sleep level is an integer with three possible values (1, 2, 3) while in ``v1.2`` is a string. We convert integer levels to strings, ``asleep``,``restless`` or ``awake`` respectively.
    - Count summaries. For Summary data, ``v1`` contains ``count_awake``, ``duration_awake``, ``count_awakenings``, ``count_restless``, and ``duration_restless`` fields for every sleep record while ``v1.2`` does not.
    - Types of sleep records. ``v1.2`` has two types of sleep records: ``classic`` and ``stages``. The ``classic`` type contains three sleep levels: ``awake``, ``restless`` and ``asleep``. The ``stages`` type contains four sleep levels: ``wake``, ``deep``, ``light``, and ``rem``. Sleep records from ``v1`` will have the same sleep levels as `v1.2` classic type; therefore we set their type to ``classic``.
    - Unified level of sleep. For intraday data, we unify sleep levels of each sleep record with a column named ``unified_level``. Based on `this Fitbit forum post`_ , we merge levels into two categories:
      - For the ``classic`` type unified_level is one of {0, 1} where 0 means awake and groups ``awake`` + ``restless``, while 1 means asleep and groups ``asleep``.
      - For the ``stages`` type, unified_level is one of {0, 1} where 0 means awake and groups ``wake`` while 1 means asleep and groups ``deep`` + ``light`` + ``rem``.
    - Short Data. In ``v1.2``, records of type ``stages`` contain ``shortData`` in addition to ``data``. We merge both to extract intraday data. 
      - ``data`` contains sleep stages and any wake periods > 3 minutes (180 seconds).
      - ``shortData`` contains short wake periods representing physiological awakenings that are <= 3 minutes (180 seconds).
    - The following columns of Summary data are not computed by RAPIDS but taken directly from columns with a similar name provided by Fitbit's API: ``efficiency``, ``minutes_after_wakeup``, ``minutes_asleep``, ``minutes_awake``, ``minutes_to_fall_asleep``, ``minutes_in_bed``, ``is_main_sleep`` and ``type``
    - The following columns of Intraday data are not computed by RAPIDS but taken directly from columns with a similar name provided by Fitbit's API: ``original_level``, ``is_main_sleep`` and ``type``. We compute ``unified_level`` as explained above.

These are examples of intraday and summary data:

- Intraday data (at 30-second intervals for ``stages`` type or 60-second intervals for ``classic`` type)

=========    ==============    =============    =============    ======    ===================    ==========    ===========    =========    =================    ==========    ==========    ============    =================
device_id    original_level    unified_level    is_main_sleep    type      local_date_time        local_date    local_month    local_day    local_day_of_week    local_time    local_hour    local_minute    local_day_segment
=========    ==============    =============    =============    ======    ===================    ==========    ===========    =========    =================    ==========    ==========    ============    =================
did          wake              0                1                stages    2020-05-20 22:13:30    2020-05-20    5              20           2                    22:13:30      22            13              evening
did          wake              0                1                stages    2020-05-20 22:14:00    2020-05-20    5              20           2                    22:14:00      22            14              evening
did          light             1                1                stages    2020-05-20 22:14:30    2020-05-20    5              20           2                    22:14:30      22            14              evening
did          light             1                1                stages    2020-05-20 22:15:00    2020-05-20    5              20           2                    22:15:00      22            15              evening
did          light             1                1                stages    2020-05-20 22:15:30    2020-05-20    5              20           2                    22:15:30      22            15              evening
=========    ==============    =============    =============    ======    ===================    ==========    ===========    =========    =================    ==========    ==========    ============    =================

- Summary data

=========    ==========    ====================    ==============    =============    ======================    ==============    =============    ======    =====================    ===================    ================    ==============    =======================    =====================
device_id    efficiency    minutes_after_wakeup    minutes_asleep    minutes_awake    minutes_to_fall_asleep    minutes_in_bed    is_main_sleep    type      local_start_date_time    local_end_date_time    local_start_date    local_end_date    local_start_day_segment    local_end_day_segment
=========    ==========    ====================    ==============    =============    ======================    ==============    =============    ======    =====================    ===================    ================    ==============    =======================    =====================
did          90            0                       381               54               0                         435               1                stages    2020-05-20 22:12:00      2020-05-21 05:27:00    2020-05-20          2020-05-21        evening                    night
did          88            0                       498               86               0                         584               1                stages    2020-05-22 22:03:00      2020-05-23 07:47:03    2020-05-22          2020-05-23        evening                    morning
=========    ==========    ====================    ==============    =============    ======================    ==============    =============    ======    =====================    ===================    ================    ==============    =======================    =====================


.. _fitbit-heart-rate-sensor-doc:

Fitbit: Heart Rate
"""""""""""""""""""

See `Fitbit: Heart Rate Config Code`_

**Available Epochs (day_segment) :** daily, morning, afternoon, evening, night

**Available Platforms:**: Fitbit

**Snakemake rule chain:**

- Rule ``rules/preprocessing.snakefile/download_dataset``
- Rule ``rules/preprocessing.snakefile/fitbit_with_datetime``
- Rule ``rules/features.snakefile/fitbit_heartrate_features``
    
.. _fitbit-heart-rate-parameters:

**Fitbit: Heart Rate Rule Parameters (fitbit_heartrate_features):**

============    ===================
Name	        Description
============    ===================
day_segment     The particular ``day_segment`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
features        The heartrate features that can be computed. See :ref:`Available Fitbit: Heart Rate Features <fitbit-heart-rate-available-features>` Table below
============    ===================

.. _fitbit-heart-rate-available-features:

**Available Fitbit: Heart Rate Features**

==================   ===========    =============
Name                 Units          Description
==================   ===========    =============
restingheartrate     beats/mins     The number of times your heart beats per minute when participant is still and well rested for ``daily`` epoch.
calories             cals           Calories burned during ``heartrate_zone`` for ``daily`` epoch.
maxhr                beats/mins     The maximum heart rate during ``day_segment`` epoch.
minhr                beats/mins     The minimum heart rate during ``day_segment`` epoch.
avghr                beats/mins     The average heart rate during ``day_segment`` epoch.
medianhr             beats/mins     The median of heart rate during ``day_segment`` epoch.
modehr               beats/mins     The mode of heart rate during ``day_segment`` epoch.
stdhr                beats/mins     The standard deviation of heart rate during ``day_segment`` epoch.
diffmaxmodehr        beats/mins     The difference between the maximum and mode heart rate during ``day_segment`` epoch.
diffminmodehr        beats/mins     The difference between the mode and minimum heart rate during ``day_segment`` epoch.
entropyhr            nats           Shannon’s entropy measurement based on heart rate during ``day_segment`` epoch.
minutesonZONE        minutes        Number of minutes the user's heartrate fell within each ``heartrate_zone`` during ``day_segment`` epoch.
==================   ===========    =============

**Assumptions/Observations:** 

There are four heart rate zones: ``out_of_range``, ``fat_burn``, ``cardio``, and ``peak``. Please refer to `Fitbit documentation`_ for more information about the way they are computed.

Calories' accuracy depends on the users’ Fitbit profile (weight, height, etc.).


.. _fitbit-steps-sensor-doc:

Fitbit: Steps
"""""""""""""""

See `Fitbit: Steps Config Code`_

**Available Epochs (day_segment) :** daily, morning, afternoon, evening, night

**Available Platforms:**: Fitbit

**Snakemake rule chain:**

- Rule ``rules/preprocessing.snakefile/download_dataset``
- Rule ``rules/preprocessing.snakefile/fitbit_with_datetime``
- Rule ``rules/features.snakefile/fitbit_step_features``
    
.. _fitbit-steps-parameters:

**Fitbit: Steps Rule Parameters (fitbit_step_features):**

=======================    ===================
Name	                   Description
=======================    ===================
day_segment                The particular ``day_segment`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
features                   The features that can be computed. See :ref:`Available Fitbit: Steps Features <fitbit-steps-available-features>` Table below
threshold_active_bout      Every minute with Fitbit step data wil be labelled as ``sedentary`` if its step count is below this threshold, otherwise, ``active``. 
include_zero_step_rows     Whether or not to include day segments with a 0 step count
exclude_sleep              Whether or not to exclude step rows that happen during sleep
exclude_sleep_type         If ``exclude_sleep`` is True, then you can choose between ``FIXED`` or ``FITBIT_BASED``. ``FIXED`` will exclude all step rows that happen between a start and end time (see below). ``FITBIT_BASED`` will exclude step rows that happen during main sleep segments as measured by the Fitbit device (``config[SLEEP][DB_TABLE]`` should be a valid table in your database, it usually is the same table that contains your STEP data)
exclude_sleep_fixed_start  Start time of the fixed sleep period to exclude. Only relevant if ``exclude_sleep`` is True and ``exclude_sleep_type`` is ``FIXED``
exclude_sleep_fixed_end    Start time of the fixed sleep period to exclude. Only relevant if ``exclude_sleep`` is True and ``exclude_sleep_type`` is ``FIXED``
=======================    ===================

.. _fitbit-steps-available-features:

**Available Fitbit: Steps Features**

=========================   =========     =============
Name                        Units         Description
=========================   =========     =============
sumallsteps                 steps         The total step count during ``day_segment`` epoch.
maxallsteps                 steps         The maximum step count during ``day_segment`` epoch.
minallsteps                 steps         The minimum step count during ``day_segment`` epoch.
avgallsteps                 steps         The average step count during ``day_segment`` epoch.
stdallsteps                 steps         The standard deviation of step count during ``day_segment`` epoch.
countsedentarybout          bouts         Number of sedentary bouts during ``day_segment`` epoch.
maxdurationsedentarybout    minutes       The maximum duration of any sedentary bout during ``day_segment`` epoch.
mindurationsedentarybout    minutes       The minimum duration of any sedentary bout during ``day_segment`` epoch.
avgdurationsedentarybout    minutes       The average duration of sedentary bouts during ``day_segment`` epoch.
stddurationsedentarybout    minutes       The standard deviation of the duration of sedentary bouts during ``day_segment`` epoch.
countactivebout             bouts         Number of active bouts during ``day_segment`` epoch.
maxdurationactivebout       minutes       The maximum duration of any active bout during ``day_segment`` epoch.
mindurationactivebout       minutes       The minimum duration of any active bout during ``day_segment`` epoch.
avgdurationactivebout       minutes       The average duration of active bouts during ``day_segment`` epoch.
stddurationactivebout       minutes       The standard deviation of the duration of active bouts during ``day_segment`` epoch.
=========================   =========     =============

**Assumptions/Observations:** 

Active and sedentary bouts. If the step count per minute is smaller than ``THRESHOLD_ACTIVE_BOUT`` (default value is 10), that minute is labelled as sedentary, otherwise, is labelled as active. Active and sedentary bouts are periods of consecutive minutes labelled as ``active`` or ``sedentary``.
	

.. -------------------------Links ------------------------------------ ..

.. _TABLES_FOR_SENSED_BINS: https://github.com/carissalow/rapids/blob/0c53fd275e628819cf79cf5b87006ce1ad9e597c/config.yaml#L3
.. _`Messages Config Code`: https://github.com/carissalow/rapids/blob/0c53fd275e628819cf79cf5b87006ce1ad9e597c/config.yaml#L35
.. _AWARE: https://awareframework.com/what-is-aware/
.. _`List of Timezones`: https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
.. _DAY_SEGMENTS: https://github.com/carissalow/rapids/blob/0c53fd275e628819cf79cf5b87006ce1ad9e597c/config.yaml#L10
.. _PHONE_VALID_SENSED_DAYS: https://github.com/carissalow/rapids/blob/0c53fd275e628819cf79cf5b87006ce1ad9e597c/config.yaml#L61
.. _`Call Config Code`: https://github.com/carissalow/rapids/blob/0c53fd275e628819cf79cf5b87006ce1ad9e597c/config.yaml#L45
.. _`WiFi Config Code`: https://github.com/carissalow/rapids/blob/0c53fd275e628819cf79cf5b87006ce1ad9e597c/config.yaml#L169
.. _`Bluetooth Config Code`: https://github.com/carissalow/rapids/blob/0c53fd275e628819cf79cf5b87006ce1ad9e597c/config.yaml#L81
.. _`Accelerometer Config Code`: https://github.com/carissalow/rapids/blob/0c53fd275e628819cf79cf5b87006ce1ad9e597c/config.yaml#L115
.. _`Applications Foreground Config Code`: https://github.com/carissalow/rapids/blob/0c53fd275e628819cf79cf5b87006ce1ad9e597c/config.yaml#L125
.. _`Battery Config Code`: https://github.com/carissalow/rapids/blob/0c53fd275e628819cf79cf5b87006ce1ad9e597c/config.yaml#L95
.. _`Activity Recognition Config Code`: https://github.com/carissalow/rapids/blob/0c53fd275e628819cf79cf5b87006ce1ad9e597c/config.yaml#L87
.. _`Light Config Code`: https://github.com/carissalow/rapids/blob/0c53fd275e628819cf79cf5b87006ce1ad9e597c/config.yaml#L109
.. _`Location (Barnett’s) Config Code`: https://github.com/carissalow/rapids/blob/0c53fd275e628819cf79cf5b87006ce1ad9e597c/config.yaml#L71
.. _`Screen Config Code`: https://github.com/carissalow/rapids/blob/0c53fd275e628819cf79cf5b87006ce1ad9e597c/config.yaml#L101
.. _`Fitbit: Sleep Config Code`: https://github.com/carissalow/rapids/blob/0c53fd275e628819cf79cf5b87006ce1ad9e597c/config.yaml#L162
.. _`version 1`: https://dev.fitbit.com/build/reference/web-api/sleep-v1/
.. _`version 1.2`: https://dev.fitbit.com/build/reference/web-api/sleep/
.. _`Conversation Config Code`: https://github.com/carissalow/rapids/blob/0c53fd275e628819cf79cf5b87006ce1ad9e597c/config.yaml#L191
.. _`this Fitbit forum post`: https://community.fitbit.com/t5/Alta/What-does-Restless-mean-in-sleep-tracking/td-p/2989011
.. _shortData: https://dev.fitbit.com/build/reference/web-api/sleep/#interpreting-the-sleep-stage-and-short-data
.. _`Fitbit: Heart Rate Config Code`: https://github.com/carissalow/rapids/blob/0c53fd275e628819cf79cf5b87006ce1ad9e597c/config.yaml#L138
.. _`Fitbit: Steps Config Code`: https://github.com/carissalow/rapids/blob/0c53fd275e628819cf79cf5b87006ce1ad9e597c/config.yaml#L145
.. _`Fitbit documentation`: https://help.fitbit.com/articles/en_US/Help_article/1565
.. _top1global: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L108
.. _`Beiwe Summary Statistics`: http://wiki.beiwe.org/wiki/Summary_Statistics
.. _`Pause-Flight Model`: https://academic.oup.com/biostatistics/advance-article/doi/10.1093/biostatistics/kxy059/5145908
