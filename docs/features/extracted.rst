.. _rapids_features:

RAPIDS Features
===============

Global Parameters
"""""""""""""""""

.. _sensor-list:

- ``SENSORS`` - List of sensors to include in the pipeline that have to match existent tables in your AWARE_ database. See SENSORS_ variable in ``config`` file.  

.. _fitbit-table:

- ``FITBIT_TABLE`` - The name of table in your database that contains Fitbit data. Its ``fitbit_data`` field should contain the data coming from the Fitbit API in JSON format. 

.. _fitbit-sensors:

- ``FITBIT_SENSORS`` - The list of sensors to be parsed from the fitbit table: ``heartrate``, ``steps``, ``sleep``.

.. _pid: 

- ``PID`` - The list of participant ids to be included in the analysis. These should match the names of the files created in the ``data/external`` directory  (:ref:`see more details<db-configuration>`).

.. _day-segments: 

- ``DAY_SEGMENTS`` - The list of day epochs that features can be segmented into: ``daily``, ``morning`` (6am-12pm), ``afternnon`` (12pm-6pm), ``evening`` (6pm-12am) and ``night`` (12am-6am). This list can be modified globally or on a per sensor basis. See DAY_SEGMENTS_ in ``config`` file.

.. _timezone:

- ``TIMEZONE`` - The time zone where data was collected. Use the timezone names from this `List of Timezones`_. Double check your chosen name is correct, for example US Eastern Time is named New America/New_York, not EST.

.. _database_group:

- ``DATABASE_GROUP`` - The name of your database credentials group, it should match the one in ``.env`` (:ref:`see the datbase configuration<db-configuration>`). 

.. _download-dataset:

- ``DOWNLOAD_DATASET``

    - ``GROUP``. Credentials group to connect to the database containing ``SENSORS``. By default it points to ``DATABASE_GROUP``.

.. _readable-datetime:

- ``READABLE_DATETIME`` - Configuration to convert UNIX timestamps into readbale date time strings.

    - ``FIXED_TIMEZONE``. See ``TIMEZONE`` above. This assumes that all data of all participants was collected within one time zone.
    - Support for multiple time zones for each participant coming soon.

.. _phone-valid-sensed-days:

- ``PHONE_VALID_SENSED_DAYS``.
    
    Contains three attributes: ``BIN_SIZE``, ``MIN_VALID_HOURS``, ``MIN_BINS_PER_HOUR``. 

    On any given day, Aware could have sensed data only for a few minutes or for 24 hours. Daily estimates of features should be considered more reliable the more hours Aware was running and logging data (for example, 10 calls logged on a day when only one hour of data was recorded is a less reliable measurement compared to 10 calls on a day when 23 hours of data were recorded. 

    Therefore, we define a valid hour as those that contain at least a certain number of valid bins. In turn, a valid bin are those that contain at least one row of data from any sensor logged within that period. We divide an hour into N bins of size ``BIN_SIZE`` (in minutes) and we mark an hour as valid if contains at least ``MIN_BINS_PER_HOUR`` of valid bins (out of the total possible number of bins that can be captured in an hour i.e. out of 60min/``BIN_SIZE`` bins). Days with valid sensed hours less than ``MIN_VALID_HOURS`` will be excluded form the output of this file. See PHONE_VALID_SENSED_DAYS_ in ``config.yaml``.

    In RAPIDS, you will find that we use ``phone_sensed_bins`` (a list of all valid and invalid bins of all monitored days) to improve the estimation of features that are ratios over time periods like ``episodepersensedminutes`` of :ref:`Screen<screen-sensor-doc>`.


.. _individual-sensor-settings:


.. _sms-sensor-doc:

SMS
"""""

See `SMS Config Code`_

**Available Epochs:**      

- daily 
- morning
- afternoon
- evening
- night

**Available Platforms:**    

- Android

**Snakefile Entry:**

..    - Download raw SMS dataset: ``expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``

..    - Apply readable datetime to SMS dataset: ``expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``

- Extract SMS features:

      | ``expand("data/processed/{pid}/sms_{sms_type}_{day_segment}.csv".``
      |                     ``pid=config["PIDS"],``
      |                     ``sms_type = config["SMS"]["TYPES"],``
      |                     ``day_segment = config["SMS"]["DAY_SEGMENTS"]),``

**Rule Chain:**

- **Rule:** ``rules/preprocessing.snakefile/download_dataset`` - See the download_dataset_ rule.

    - **Script:** ``src/data/download_dataset.R`` - See the download_dataset.R_ script.
    
- **Rule:** ``rules/preprocessing.snakefile/readable_datetime`` - See the readable_datetime_ rule.

    - **Script:** ``src/data/readable_datetime.R`` - See the readable_datetime.R_ script.

- **Rule:** ``rules/features.snakefile/sms_features`` - See the sms_features_ rule.

    - **Script:** ``src/features/sms_features.R`` - See the sms_features.R_ script.


.. _sms-parameters:

**SMS Rule Parameters:**

============    ===================
Name	        Description
============    ===================
sms_type        The particular ``sms_type`` that will be analyzed. The options for this parameter are ``received`` or ``sent``.
day_segment     The particular ``day_segment`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
features        The different measures that can be retrieved from the dataset. These features are available for both ``sent`` and ``received`` SMS messages. See :ref:`Available SMS Features <sms-available-features>` Table below
============    ===================

.. _sms-available-features:

**Available SMS Featues**

The following table shows a list of the available featues for both ``sent`` and ``received`` SMS. 

=========================   =========     =============
Name                        Units         Description
=========================   =========     =============
count                       SMS           Number of SMS of type ``sms_type`` that occurred during a particular ``day_segment``.
distinctcontacts            contacts      Number of distinct contacts that are associated with a particular ``sms_type`` during a particular ``day_segment``.
timefirstsms                minutes       Number of minutes between 12:00am (midnight) and the first ``SMS`` of a particular ``sms_type``.
timelastsms                 minutes       Number of minutes between 12:00am (midnight) and the last ``SMS`` of a particular ``sms_type``.
countmostfrequentcontact    SMS           The count of the number of ``SMS`` messages of a particular ``sms_type`` for the most contacted contact for a particular ``day_segment``.
=========================   =========     =============

**Assumptions/Observations:** 

    #. ``TYPES`` and ``FEATURES`` keys need to match. From example::

        SMS:
            TYPES: [sent]
            FEATURES: 
                sent: [count, distinctcontacts, timefirstsms, timelastsms, countmostfrequentcontact]

In the above config setting code the ``TYPE`` ``sent`` matches the ``FEATURES`` key ``sent``.


.. _call-sensor-doc:

Calls
""""""

See `Call Config Code`_

**Available Epochs:**      

- daily 
- morning
- afternoon
- evening
- night

**Available Platforms:**    

- Android
- iOS

**Snakefile Entry:**

..    - Download raw Calls dataset: ``expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``

..    - Apply readable datetime to Calls dataset: ``expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``
    
- Extract Calls Features
    
      | ``expand("data/processed/{pid}/call_{call_type}_{segment}.csv",``
      |                      ``pid=config["PIDS"],`` 
      |                      ``call_type=config["CALLS"]["TYPES"],``
      |                      ``segment = config["CALLS"]["DAY_SEGMENTS"]),``
    
**Rule Chain:**

- **Rule:** ``rules/preprocessing.snakefile/download_dataset`` - See the download_dataset_ rule.

    - **Script:** ``src/data/download_dataset.R`` - See the download_dataset.R_ script.

- **Rule:** ``rules/preprocessing.snakefile/readable_datetime`` - See the readable_datetime_ rule.

    - **Script:** ``src/data/readable_datetime.R`` - See the readable_datetime.R_ script.

- **Rule:** ``rules/features.snakefile/call_features`` - See the call_features_ rule.

    - **Script:** ``src/features/call_features.R`` - See the call_features.R_ script.

    
.. _calls-parameters:

**Call Rule Parameters:**

============    ===================
Name	        Description
============    ===================
call_type       The particular ``call_type`` that will be analyzed. The options for this parameter are ``incoming``, ``outgoing`` or ``missed``.
day_segment     The particular ``day_segment`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
features         The different measures that can be retrieved from the calls dataset. Note that the same features are available for both ``incoming`` and ``outgoing`` calls, while ``missed`` calls has its own set of features. See :ref:`Available Incoming and Outgoing Call Features <available-in-and-out-call-features>` Table and :ref:`Available Missed Call Features <available-missed-call-features>` Table below.
============    ===================

.. _available-in-and-out-call-features:

**Available Incoming and Outgoing Call Features**

The following table shows a list of the available features for ``incoming`` and ``outgoing`` calls. 

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
timefirstcall               hours         The time in hours between 12:00am (midnight) and the first call of ``call_type``.
timelastcall                hours         The time in hours between 12:00am (midnight) and the last call of ``call_type``.
countmostfrequentcontact    calls         The number of calls of a particular ``call_type`` during a particular ``day_segment`` of the most frequent contact throughout the monitored period.
=========================   =========     =============

.. _available-missed-call-features:

**Available Missed Call Features**

The following table shows a list of the available features for ``missed`` calls. 

=========================   =========     =============
Name                        Units         Description
=========================   =========     =============
count                       calls         Number of ``missed`` calls that occurred during a particular ``day_segment``.
distinctcontacts            contacts      Number of distinct contacts that are associated with ``missed`` calls for a particular ``day_segment``
timefirstcall               hours         The time in hours from 12:00am (Midnight) that the first ``missed`` call occurred.
timelastcall                hours         The time in hours from 12:00am (Midnight) that the last ``missed`` call occurred.
countmostfrequentcontact    calls         The number of ``missed`` calls during a particular ``day_segment`` of the most frequent contact throughout the monitored period.
=========================   =========     =============

**Assumptions/Observations:** 

    #. ``TYPES`` and ``FEATURES`` keys need to match. From example::

        CALLS:
            TYPES: [missed]
            FEATURES: 
                missed: [count, distinctcontacts, timefirstcall, timelastcall, countmostfrequentcontact]

In the above config setting code the ``TYPE`` ``missed`` matches the ``FEATURES`` key ``missed``.


.. _bluetooth-sensor-doc:

Bluetooth
""""""""""

See `Bluetooth Config Code`_

**Available Epochs:**      

- daily 
- morning
- afternoon
- evening
- night

**Available Platforms:**    

- Android
- iOS

**Snakefile Entry:**

..    - Download raw Bluetooth dataset: ``expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``

..    - Apply readable datetime to Bluetooth dataset: ``expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``
    
- Extract Bluetooth Features
    
      | ``expand("data/processed/{pid}/bluetooth_{segment}.csv",``
      |          ``pid=config["PIDS"],`` 
      |          ``segment = config["BLUETOOTH"]["DAY_SEGMENTS"]),``
    
**Rule Chain:**

- **Rule:** ``rules/preprocessing.snakefile/download_dataset`` - See the download_dataset_ rule.

    - **Script:** ``src/data/download_dataset.R`` See the download_dataset.R_ script.

- **Rule:** ``rules/preprocessing.snakefile/readable_datetime`` - See the readable_datetime_ rule.

    - **Script:** ``src/data/readable_datetime.R`` See the readable_datetime.R_ script.

- **Rule:** ``rules/features.snakefile/bluetooth_features`` - See the bluetooth_feature_ rule.

    - **Script:** ``src/features/bluetooth_features.R`` - See the bluetooth_features.R_ script.

    
.. _bluetooth-parameters:

**Bluetooth Rule Parameters:**

============    ===================
Name	        Description
============    ===================
day_segment     The particular ``day_segment`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
features        The different measures that can be retrieved from the Bluetooth dataset. See :ref:`Available Bluetooth Features <bluetooth-available-features>` Table below
============    ===================

.. _bluetooth-available-features:

**Available Bluetooth Features**

The following table shows a list of the available features for Bluetooth. 

===========================   =========     =============
Name                          Units         Description
===========================   =========     =============
countscans                    devices       Number of scanned devices during a ``day_segment``, a device can be detected multiple times over time and these appearances are counted separately
uniquedevices                 devices       Number of unique devices during a ``day_segment`` as identified by their hardware address
countscansmostuniquedevice    scans         Number of scans of the most scanned device during a ``day_segment`` across the whole monitoring period
===========================   =========     =============

**Assumptions/Observations:** N/A 



.. _accelerometer-sensor-doc:

Accelerometer
""""""""""""""

See `Accelerometer Config Code`_

**Available epochs:**      

- daily 
- morning
- afternoon
- evening
- night

**Available platforms:**    

- Android
- iOS

**Snakefile entry:**

..  - Download raw Accelerometer dataset: ``expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``

..  - Apply readable datetime to Accelerometer dataset: ``expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``

- Extract Accelerometer Features

    | ``expand("data/processed/{pid}/accelerometer_{day_segment}.csv",``
    |                      ``pid=config["PIDS"],`` 
    |                      ``day_segment = config["ACCELEROMETER"]["DAY_SEGMENTS"]),``

**Rule chain:**

- **Rule:** ``rules/preprocessing.snakefile/download_dataset`` - See the download_dataset_ rule.

    - **Script:** ``src/data/download_dataset.R`` - See the download_dataset.R_ script.

- **Rule:** ``rules/preprocessing.snakefile/readable_datetime`` - See the readable_datetime_ rule.

    - **Script:** ``src/data/readable_datetime.R`` - See the readable_datetime.R_ script.

- **Rule:** ``rules/features.snakefile/accelerometer_features`` - See the accelerometer_features_ rule.

    - **Script:** ``src/features/accelerometer_features.py`` - See the accelerometer_features.py_ script.

    
.. _Accelerometer-parameters:

**Accelerometer Rule Parameters:**

============    ===================
Name	        Description
============    ===================
day_segment     The particular ``day_segment`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
features        The different measures that can be retrieved from the dataset. See :ref:`Available Accelerometer Features <accelerometer-available-features>` Table below
============    ===================

.. _accelerometer-available-features:

**Available Accelerometer Features**

The following table shows a list of the available features the accelerometer sensor data for a particular ``day_segment``. 

====================================   ==============    =============
Name                                   Units             Description
====================================   ==============    =============
maxmagnitude                           m/s\ :sup:`2`     The maximum magnitude of acceleration (:math:`\|acceleration\| = \sqrt{x^2 + y^2 + z^2}`).
minmagnitude                           m/s\ :sup:`2`     The minimum magnitude of acceleration.
avgmagnitude                           m/s\ :sup:`2`     The average magnitude of acceleration.
medianmagnitude                        m/s\ :sup:`2`     The median magnitude of acceleration.
stdmagnitude                           m/s\ :sup:`2`     The standard deviation of acceleration.
ratioexertionalactivityepisodes                          The ratio of exertional activity time periods to total time periods.
sumexertionalactivityepisodes          minutes           Total duration of all exertional activity episodes during ``day_segment``.
longestexertionalactivityepisode       minutes           Duration of the longest exertional activity episode during ``day_segment``.
longestnonexertionalactivityepisode    minutes           Duration of the longest non-exertional activity episode during ``day_segment``.
countexertionalactivityepisodes        episodes          Number of the exertional activity episodes during ``day_segment``.
countnonexertionalactivityepisodes     episodes          Number of the non-exertional activity episodes during ``day_segment``.
====================================   ==============    =============

**Assumptions/Observations:** N/A



.. _applications-foreground-sensor-doc:

Applications Foreground
""""""""""""""""""""""""

See `Applications Foreground Config Code`_

**Available Epochs:**      

- daily 
- morning
- afternoon
- evening
- night

**Available Platforms:**    

- Android

**Snakefile entry:**

..  - Download raw Applications Foreground dataset: ``expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``

..  - Apply readable dateime Applications Foreground dataset: ``expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``
    
..  - Genre categorization of Applications Foreground dataset: ``expand("data/interim/{pid}/applications_foreground_with_datetime_with_genre.csv", pid=config["PIDS"]),``

- Extract Applications Foreground Features:

    | ``expand("data/processed/{pid}/applications_foreground_{day_segment}.csv",``
    |                      ``pid=config["PIDS"],`` 
    |                      ``day_segment = config["APPLICATIONS_FOREGROUND"]["DAY_SEGMENTS"]),``

**Rule Chain:**

- **Rule:** ``rules/preprocessing.snakefile/download_dataset`` - See the download_dataset_ rule.

        - **Script:** ``src/data/download_dataset.R`` - See the download_dataset.R_ script.

- **Rule:** ``rules/preprocessing.snakefile/readable_datetime`` - See the readable_datetime_ rule.

    - **Script:** ``src/data/readable_datetime.R`` - See the readable_datetime.R_ script.

- **Rule:** ``rules/preprocessing.snakefile/application_genres`` - See the application_genres_ rule

    - **Script:** ``../src/data/application_genres.R`` - See the application_genres.R_ script

- **Rule:** ``rules/features.snakefile/applications_foreground_features`` - See the applications_foreground_features_ rule.

    - **Script:** ``src/features/applications_foreground_features.py`` - See the applications_foreground_features.py_ script.
   
.. _applications-foreground-parameters:

**Applications Foreground Rule Parameters:**

====================    ===================
Name	                Description
====================    ===================
day_segment             The particular ``day_segment`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
single_categories       App categories to be included in the feature extraction computation. See ``APPLICATION_GENRES`` in this file to add new categories or use the catalogue we provide and read :ref:`Assumtions and Observations <applications-foreground-observations>` for more information.
multiple_categories     You can group multiple categories into meta categories, for example ``social: ["socialnetworks", "socialmediatools"]``.
single_apps             Apps to be included in the feature extraction computation. Use their package name, for example, ``com.google.android.youtube`` or the reserved word ``top1global`` (the most used app by a participant over the whole monitoring study).
excluded_categories     App categories to be excluded in the feature extraction computation. See ``APPLICATION_GENRES`` in this file to add new categories or use the catalogue we provide and read :ref:`Assumtions and Observations <applications-foreground-observations>` for more information.
excluded_apps           Apps to be excluded in the feature extraction computation. Use their package name, for example: ``com.google.android.youtube``
features                The features to be extracted. See :ref:`Available Applications Foreground Features <applications-foreground-available-features>` Table below
====================    ===================

.. _applications-foreground-available-features:

**Available Applications Foreground Features**

The following table shows a list of the available features for the Applications Foreground dataset 

==================   =========   =============
Name                 Units       Description
==================   =========   =============
count                apps        Number of times a single app or apps within a category were used (i.e. they were brought to the foreground either by tapping their icon or switching to it from another app).
timeoffirstuse       contacts    The time in minutes between 12:00am (midnight) and the first use of a single app or apps within a category during a ``day_segment``.
timeoflastuse        minutes     The time in minutes between 12:00am (midnight) and the last use of a single app or apps within a category during a ``day_segment``.
frequencyentropy     nats        The entropy of the used apps within a category during a ``day_segment`` (each app is seen as a unique event, the more apps were used, the higher the entropy). This is especially relevant when computed over all apps. Entropy cannot be obtained for a single app.
==================   =========   =============

.. _applications-foreground-observations:

**Assumptions/Observations:** 

Features can be computed by app, by apps grouped under a single category (genre) and by multiple categories grouped together (meta categories). For example, we can get features for Facebook, for Social Network Apps (including Facebook and others) or for a meta category called Social formed by Social Network and Social Media Tools categories. 

We provide three ways of classifying and app within a category (genre): a) by automatically scraping its official category from the Google Play Store, b) by using the catalogue created by Stachl et al. which we provide in RAPIDS (``data/external/``), or c) by manually creating a personalized catalogue.

The way you choose strategy a, b or c is by modifying ``APPLICATION_GENRES`` keys and values. Set ``CATALOGUE_SOURCE`` to ``FILE`` if you want to use a CSV file as catalogue or to ``GOOGLE`` if you want to scrape the genres from the Play Store. By default ``CATALOGUE_FILE`` points to the catalogue created by  Stachl et al. and you can change this path to your own catalogue that follows the same format. In addition, set ``SCRAPE_MISSING_GENRES`` to true if you are using a FILE catalogue and you want to scrape from the Play Store any missing genres and ``UPDATE_CATALOGUE_FILE`` to true if you want to save those scrapped genres back into the FILE.

.. _battery-sensor-doc:

Battery
"""""""""

See `Battery Config Code`_

**Available Epochs:**      

- daily 
- morning
- afternoon
- evening
- night

**Available Platforms:**    

- Android
- iOS

**Snakefile entry:**

..  - Download raw Battery dataset: ``expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``

..  - Apply readable dateime to Battery dataset: ``expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``
    
..  - Extract the deltas in Battery charge : ``expand("data/processed/{pid}/battery_deltas.csv", pid=config["PIDS"]),``

- Extract Battery Features:

    | ``expand("data/processed/{pid}/battery_{day_segment}.csv",``
    |                      ``pid=config["PIDS"],`` 
    |                      ``day_segment = config["BATTERY"]["DAY_SEGMENTS"]),``
    
**Rule Chain:**

- **Rule:** ``rules/preprocessing.snakefile/download_dataset`` - See the download_dataset_ rule.

        - **Script:** ``src/data/download_dataset.R`` - See the download_dataset.R_ script.

- **Rule:** ``rules/preprocessing.snakefile/readable_datetime`` - See the readable_datetime_ rule.

    - **Script:** ``src/data/readable_datetime.R`` - See the readable_datetime.R_ script.

- **Rule:** ``rules/features.snakefile/battery_deltas`` - See the battery_deltas_ rule.

    - **Script:** ``src/features/battery_deltas.R`` - See the battery_deltas.R_ script.
    
- **Rule:** ``rules/features.snakefile/battery_features`` - See the battery_features_ rule

    - **Script:** ``src/features/battery_features.py`` - See the battery_features.py_ script.
    
.. _battery-parameters:

**Battery Rule Parameters:**

============    ===================
Name	        Description
============    ===================
day_segment     The particular ``day_segment`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
features        The different measures that can be retrieved from the Battery dataset. See :ref:`Available Battery Features <battery-available-features>` Table below
============    ===================

.. _battery-available-features:

**Available Battery Features**

The following table shows a list of the available features for Battery data. 

=====================   ===============   =============
Name                    Units             Description
=====================   ===============   =============
countdischarge          episodes          Number of discharging episodes.
sumdurationdischarge    hours             The total duration of all discharging episodes.
countcharge             episodes          Number of battery charging episodes.
sumdurationcharge       hours             The total duration of all charging episodes.
avgconsumptionrate      episodes/hours    The average of all episodes’ consumption rates. An episode’s consumption rate is defined as the ratio between its battery delta and duration
maxconsumptionrate      episodes/hours    The highest of all episodes’ consumption rates. An episode’s consumption rate is defined as the ratio between its battery delta and duration
=====================   ===============   =============

**Assumptions/Observations:** 


.. _activity-recognition-sensor-doc:

Activity Recognition
""""""""""""""""""""""""""""

**Available Epochs:** daily, morning, afternoon, evening, night

**Available Platforms:** Android and iOS

**Snakefile entry to compute these features:**

    | expand("data/processed/{pid}/activity_recognition_{segment}.csv",pid=config["PIDS"], 
    |                        segment = config["ACTIVITY_RECOGNITION"]["DAY_SEGMENTS"]),
    
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

**Available Epochs:**      

    - daily 
    - morning
    - afternoon
    - evening
    - night

**Available Platforms:**    

    - Android

**Snakefile entry:**

..    - Download raw Sensor dataset: ``expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``

..    - Apply readable dateime to Sensor dataset: ``expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``
    
- Extract Light Features:

    | ``expand("data/processed/{pid}/light_{day_segment}.csv",``
    |                      ``pid=config["PIDS"],`` 
    |                      ``day_segment = config["LIGHT"]["DAY_SEGMENTS"]),``
    
**Rule Chain:**

- **Rule:** ``rules/preprocessing.snakefile/download_dataset`` - See the download_dataset_ rule.

    - **Script:** ``src/data/download_dataset.R`` - See the download_dataset.R_ script.

- **Rule:** ``rules/preprocessing.snakefile/readable_datetime`` - See the readable_datetime_ rule.

    - **Script:** ``src/data/readable_datetime.R`` - See the readable_datetime.R_ script.

- **Rule:** ``rules/features.snakefile/light_features`` - See the light_features_ rule.

    - **Script:** ``src/features/light_features.py`` - See the light_features.py_ script.

.. _light-parameters:

**Light Rule Parameters:**

============    ===================
Name	        Description
============    ===================
day_segment     The particular ``day_segment`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
features        The different measures that can be retrieved from the Light dataset. See :ref:`Available Light Features <light-available-features>` Table below
============    ===================

.. _light-available-features:

**Available Light Features**

The following table shows a list of the available features for the Light dataset. 

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
are computed (https://arxiv.org/abs/1606.06328)

See `Location (Barnett’s) Config Code`_

**Available Epochs:**      

    - daily 

**Available Platforms:**    

    - Android
    - iOS

**Snakefile entry:**

..    - Download raw Sensor dataset: ``expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``

..    - Apply readable dateime to Sensor dataset: ``expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``

- Extract Sensor Features: ``expand("data/processed/{pid}/location_barnett.csv", pid=config["PIDS"]),``
    
**Rule Chain:**

- **Rule:** ``rules/preprocessing.snakefile/download_dataset`` - See the download_dataset_ rule.

    - **Script:** ``src/data/download_dataset.R`` - See the download_dataset.R_ script.

- **Rule:** ``rules/preprocessing.snakefile/readable_datetime`` - See the readable_datetime_ rule.

    - **Script:** ``src/data/readable_datetime.R`` - See the readable_datetime.R_ script.

- **Rule:** ``rules/preprocessing.snakefile/phone_sensed_bins`` - See the phone_sensed_bins_ rule.

    - **Script:** ``src/data/phone_sensed_bins.R`` - See the phone_sensed_bins.R_ script.

- **Rule:** ``rules/preprocessing.snakefile/resample_fused_location`` - See the resample_fused_location_ rule.

    - **Script:** ``src/data/resample_fused_location.R`` - See the resample_fused_location.R_ script.

- **Rule:** ``rules/features.snakefile/location_barnett_features`` - See the location_barnett_features_ rule.

    - **Script:** ``src/features/location_barnett_features.R`` - See the location_barnett_features.R_ script.

    
.. _location-parameters:

**Location Rule Parameters:**

=================    ===================
Name	             Description
=================    ===================
location_to_use      The specifies which of the location data will be use in the analysis. Possible options are ``ALL``, ``ALL_EXCEPT_FUSED`` OR ``RESAMPLE_FUSED``
accuracy_limit       This is in meters. The sensor drops location coordinates with an accuracy higher than this. This number means there's a 68% probability the true location is within this radius specified.
timezone             The timezone used to calculate location. 
features             The different measures that can be retrieved from the Location dataset. See :ref:`Available Location Features <location-available-features>` Table below
=================    ===================

.. _location-available-features:

**Available Location Features**

The following table shows a list of the available features for Location dataset. 

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
avgflightdur       meters        Mean duration of all flights.
stdflightdur       meters        The standard deviation of the duration of all flights.
probpause                        The fraction of a day spent in a pause (as opposed to a flight)
siglocentropy      nats          Shannon’s entropy measurement based on the proportion of time spent at each significant location visited during a day.
circdnrtn           	         A continuous metric quantifying a person’s circadian routine that can take any value between 0 and 1, where 0 represents a daily routine completely different from any other sensed days and 1 a routine the same as every other sensed day.
wkenddayrtn        Weekend       Same as circdnrtn but computed separately for weekends and weekdays.
================   =========     =============

**Assumptions/Observations:** 

*Significant Locations Identified*

(i.e. The clustering method used)
Significant locations are determined using K-means clustering on locations that a patient visit over the course of the period of data collection. By setting K=K+1 and repeat clustering until two significant locations are within 100 meters of one another, the results from the previous step (K-1) can   be used as the total number of significant locations. See `Beiwe Summary Statistics`_. 

*Definition of Stationarity*

(i.e., The length of time a person have to be not moving to qualify)
This is based on a Pause-Flight model, The parameters used is a minimum pause duration of 300sec and a minimum pause distance of 60m. See the `Pause-Flight Model`_.

*The Circadian Calculation*

For a detailed description of how this measure is calculated, see Canzian and Musolesi's 2015 paper in the Proceedings of the 2015 ACM International Joint Conference on Pervasive and Ubiquitous Computing, titled "Trajectories of depression: unobtrusive monitoring of depressive states by means of smartphone mobility traces analysis." Their procedure was followed using 30-min increments as a bin size. See `Beiwe Summary Statistics`_.

.. _screen-sensor-doc:

Screen
""""""""

See `Screen Config Code`_

**Available Epochs:**      

    - daily 
    - morning
    - afternoon
    - evening
    - night

**Available Platforms:**    

    - Android
    - iOS

**Snakefile entry:**

..    - Download raw Screen dataset: ``expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``
      - Apply readable dateime to Screen dataset: ``expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``
      - Extract the deltas from the Screen dataset: expand("data/processed/{pid}/screen_deltas.csv", pid=config["PIDS"]),
    
- Extract Screen Features:
    
      | ``expand("data/processed/{pid}/screen_{day_segment}.csv",``
      |                      ``pid=config["PIDS"],`` 
      |                      ``day_segment = config["SCREEN"]["DAY_SEGMENTS"]),``
    
**Rule Chain:**

- **Rule:** ``rules/preprocessing.snakefile/download_dataset`` - See the download_dataset_ rule.

    - **Script:** ``src/data/download_dataset.R`` - See the download_dataset.R_ script.

- **Rule:** ``rules/preprocessing.snakefile/readable_datetime`` - See the readable_datetime_ rule.

    - **Script:** ``src/data/readable_datetime.R`` - See the readable_datetime.R_ script.

- **Rule:** ``rules/features.snakefile/screen_deltas`` - See the screen_deltas_ rule.

    - **Script:** ``src/features/screen_deltas.R`` - See the screen_deltas.R_ script.

- **Rule:** ``rules/features.snakefile/screen_features`` - See the screen_features_ rule.

    - **Script:** ``src/features/screen_features.py`` - See the screen_features.py_ script.

.. _screen-parameters:

**Screen Rule Parameters:**

=========================    ===================
Name	                     Description
=========================    ===================
day_segment                  The particular ``day_segments`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
reference_hour_first_use     The reference point from which ``firstuseafter`` is to be computed, default is midnight
features_deltas              The different measures that can be retrieved from the episodes extracted from the Screen dataset. See :ref:`Available Screen Episodes Features <screen-episodes-available-features>` Table below
episode_types                The action that defines an episode
=========================    ===================

.. _screen-events-available-features:

.. 
    **Available Screen Events Features**
    The following table shows a list of the available features for Screen Events. 
        =================   ==============    =============
        Name                Units             Description
        =================   ==============    =============
        counton             `ON` events       Count on: A count of screen `ON` events (only available for Android)
        countunlock         Unlock events     Count unlock: A count of screen unlock events.
        unlocksperminute    Unlock events     Unlock events per minute: The average of the number of unlock events that occur in a minute 
        =================   ==============    =============

.. _screen-episodes-available-features:

**Available Screen Episodes Features**

The following table shows a list of the available features for Screen Episodes. 

=========================   =================   =============
Name                        Units               Description
=========================   =================   =============
sumduration                 seconds             Total duration of all unlock episodes.
maxduration                 seconds             Longest duration of any unlock episode.
minduration                 seconds             Shortest duration of any unlock episode.
avgduration                 seconds             Average duration of all unlock episodes.
stdduration                 seconds             Standard deviation duration of all unlock episodes.
countepisode                episodes            Number of all unlock episodes
episodepersensedminutes     episodes/minute     The ratio between the total number of episodes in an epoch divided by the total time (minutes) the phone was sensing data.
firstuseafter               seconds             Seconds until the first unlock episode.
=========================   =================   =============

**Assumptions/Observations:** 

An ``unlock`` episode is considered as the time between an ``unlock`` event and a ``lock`` event. iOS recorded these episodes reliable (albeit some duplicated ``lock`` events within milliseconds from each other). However, in Android there are some events unrelated to the screen state because of multiple consecutive ``unlock``/``lock`` events, so we keep the closest pair. In the experiments these are less than 10% of the screen events collected. This happens because ``ACTION_SCREEN_OFF`` and ``ON`` are "sent when the device becomes non-interactive which may have nothing to do with the screen turning off". Additionally in Android it is possible to measure the time spent on the ``lock`` screen onto the ``unlock`` event and the total screen time (i.e. ``ON`` to ``OFF``) events but we are only keeping ``unlock`` episodes (``unlock`` to ``OFF``) to be consistent with iOS. 

.. ------------------------------- Begin Fitbit Section ----------------------------------- ..

.. _fitbit-sleep-sensor-doc:

Fitbit: Sleep
"""""""""""""""""""

See `Fitbit: Sleep Config Code`_

**Available Epochs:**      

    - daily 

**Available Platforms:**    

    - Fitbit

**Snakefile entry:**

- Extract Sensor Features:

    | ``expand("data/processed/{pid}/fitbit_sleep_{day_segment}.csv",``
    |                      ``pid = config["PIDS"],``
    |                      ``day_segment = config["SLEEP"]["DAY_SEGMENTS"]),``

    
**Rule Chain:**

- **Rule:** ``rules/preprocessing.snakefile/download_dataset`` - See the download_dataset_ rule.

    - **Script:** ``src/data/download_dataset.R`` - See the download_dataset.R_ script.

- **Rule:** ``rules/preprocessing.snakefile/fitbit_with_datetime`` - See the fitbit_with_datetime_ rule.

    - **Script:** ``src/data/fitbit_readable_datetime.py`` - See the fitbit_readable_datetime.py_ script.

- **Rule:** ``rules/features.snakefile/fitbit_sleep_features`` - See the fitbit_sleep_features_ rule.

    - **Script:** ``src/features/fitbit_sleep_features.py`` - See the fitbit_sleep_features.py_ script.

    
.. _fitbit-sleep-parameters:

**Fitbit: Sleep Rule Parameters:**

==================================    ===================
Name	                              Description
==================================    ===================
day_segment                           The particular ``day_segment`` that will be analyzed. For this sensor only ``daily`` is used.
sleep_types                           The types of sleep provided by Fitbit: ``main``, ``nap``, ``all``.
daily_features_from_summary_data      The sleep features that can be computed based on Fitbit's summary data. See :ref:`Available Fitbit: Sleep Features <fitbit-sleep-available-features>` Table below
==================================    ===================

.. _fitbit-sleep-available-features:

**Available Fitbit: Sleep Features**

The following table shows a list of the available features for the Fitbit: Sleep dataset. 

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

N/A



.. _fitbit-heart-rate-sensor-doc:

Fitbit: Heart Rate
"""""""""""""""""""

See `Fitbit: Heart Rate Config Code`_

**Available Epochs:**      

    - daily 
    - morning
    - afternoon
    - evening
    - night

**Available Platforms:**    

    - Fitbit

**Snakefile entry:**

..    - Download raw Fitbit: Heart Rate dataset: ``expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["FITBIT_TABLE"]),``

..    - Apply readable datetime to Fitbit: Heart Rate dataset: 

..    
      | ``expand("data/raw/{pid}/fitbit_{fitbit_sensor}_with_datetime.csv",``
      |                      ``pid=config["PIDS"],``
      |                     ``fitbit_sensor=config["FITBIT_SENSORS"]),``
      
- Extract Sensor Features:

    | ``expand("data/processed/{pid}/fitbit_heartrate_{day_segment}.csv",``
    |                      ``pid=config["PIDS"],`` 
    |                      ``day_segment = config["HEARTRATE"]["DAY_SEGMENTS"]),``
    
**Rule Chain:**

- **Rule:** ``rules/preprocessing.snakefile/download_dataset`` - See the download_dataset_ rule.

    - **Script:** ``src/data/download_dataset.R`` - See the download_dataset.R_ script.

- **Rule:** ``rules/preprocessing.snakefile/fitbit_with_datetime`` - See the fitbit_with_datetime_ rule.

    - **Script:** ``src/data/fitbit_readable_datetime.py`` - See the fitbit_readable_datetime.py_ script.

- **Rule:** ``rules/features.snakefile/fitbit_heartrate_features`` - See the fitbit_heartrate_features_ rule.

    - **Script:** ``src/features/fitbit_heartrate_features.py`` - See the fitbit_heartrate_features.py_ script.

    
.. _fitbit-heart-rate-parameters:

**Fitbit: Heart Rate Rule Parameters:**

============    ===================
Name	        Description
============    ===================
day_segment     The particular ``day_segment`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
features        The heartrate features that can be computed. See :ref:`Available Fitbit: Heart Rate Features <fitbit-heart-rate-available-features>` Table below
============    ===================

.. _fitbit-heart-rate-available-features:

**Available Fitbit: Heart Rate Features**

The following table shows a list of the available features for the Fitbit: Heart Rate dataset. 

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
lengthZONE           minutes        Number of minutes the user's heartrate fell within each ``heartrate_zone`` during ``day_segment`` epoch.
==================   ===========    =============

**Assumptions/Observations:** 

There are four heart rate zones: ``out_of_range``, ``fat_burn``, ``cardio``, and ``peak``. Please refer to `Fitbit documentation`_ for more information about the way they are computed.

Calories' accuracy depends on the users’ Fitbit profile (weight, height, etc.).

.. _fitbit-steps-sensor-doc:

Fitbit: Steps
"""""""""""""""

See `Fitbit: Steps Config Code`_

**Available Epochs:**      

    - daily 
    - morning
    - afternoon
    - evening
    - night

**Available Platforms:**    

    - Fitbit

**Snakefile entry:**

..    - Download raw Fitbit: Steps dataset: ``expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["FITBIT_TABLE"]),``

.. 
    - Apply readable datetime to Fitbit: Steps dataset: 
      | ``expand("data/raw/{pid}/fitbit_{fitbit_sensor}_with_datetime.csv",``
      |                      ``pid=config["PIDS"],``
      |                     ``fitbit_sensor=config["FITBIT_SENSORS"]),``
 
- Extract Fitbit: Steps Features:

    | ``expand("data/processed/{pid}/fitbit_step_{day_segment}.csv",``
    |                      ``pid=config["PIDS"],`` 
    |                      ``day_segment = config["STEP"]["DAY_SEGMENTS"]),``
    
**Rule Chain:**

- **Rule:** ``rules/preprocessing.snakefile/download_dataset`` - See the download_dataset_ rule.

    - **Script:** ``src/data/download_dataset.R`` - See the download_dataset.R_ script.

- **Rule:** ``rules/preprocessing.snakefile/fitbit_with_datetime`` - See the fitbit_with_datetime_ rule.

    - **Script:** ``src/data/fitbit_readable_datetime.py`` - See the fitbit_readable_datetime.py_ script.

- **Rule:** ``rules/features.snakefile/fitbit_step_features`` - See the fitbit_step_features.py_ rule.

    - **Script:** ``src/features/fitbit_step_features.py`` - See the fitbit_step_features.py_ script.

    
.. _fitbit-steps-parameters:

**Fitbit: Steps Rule Parameters:**

=======================    ===================
Name	                   Description
=======================    ===================
day_segment                The particular ``day_segment`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
features                   The features that can be computed. See :ref:`Available Fitbit: Steps Features <fitbit-steps-available-features>` Table below
threshold_active_bout      Every minute with Fitbit step data wil be labelled as ``sedentary`` if its step count is below this threshold, otherwise, ``active``. 
=======================    ===================

.. _fitbit-steps-available-features:

**Available Fitbit: Steps Features**

The following table shows a list of the available features for the Fitbit: Steps dataset. 

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

.. _SENSORS: https://github.com/carissalow/rapids/blob/f22d1834ee24ab3bcbf051bc3cc663903d822084/config.yaml#L2
.. _`SMS Config Code`: https://github.com/carissalow/rapids/blob/f22d1834ee24ab3bcbf051bc3cc663903d822084/config.yaml#L38
.. _AWARE: https://awareframework.com/what-is-aware/
.. _`List of Timezones`: https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
.. _sms_features: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L1
.. _sms_features.R: https://github.com/carissalow/rapids/blob/master/src/features/sms_featues.R
.. _download_dataset: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/preprocessing.snakefile#L9
.. _download_dataset.R: https://github.com/carissalow/rapids/blob/master/src/data/download_dataset.R
.. _readable_datetime: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/preprocessing.snakefile#L21
.. _readable_datetime.R: https://github.com/carissalow/rapids/blob/master/src/data/readable_datetime.R
.. _DAY_SEGMENTS: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L13
.. _PHONE_VALID_SENSED_DAYS: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L60
.. _`Call Config Code`: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L46
.. _call_features: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L13
.. _call_features.R: https://github.com/carissalow/rapids/blob/master/src/features/call_features.R
.. _`Bluetooth Config Code`: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L76
.. _bluetooth_feature: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L63
.. _bluetooth_features.R: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/src/features/bluetooth_features.R
.. _`Accelerometer Config Code`: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L98
.. _accelerometer_features: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L124
.. _accelerometer_features.py: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/src/features/accelerometer_featues.py
.. _`Applications Foreground Config Code`: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L102
.. _`Application Genres Config`: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L54
.. _application_genres: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/preprocessing.snakefile#L81
.. _application_genres.R: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/src/data/application_genres.R
.. _applications_foreground_features: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L135
.. _applications_foreground_features.py: https://github.com/carissalow/rapids/blob/master/src/features/accelerometer_features.py
.. _`Battery Config Code`: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L84
.. _battery_deltas: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L25
.. _battery_deltas.R: https://github.com/carissalow/rapids/blob/master/src/features/battery_deltas.R
.. _battery_features: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L86
.. _battery_features.py : https://github.com/carissalow/rapids/blob/master/src/features/battery_features.py
.. _`Google Activity Recognition Config Code`: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L80
.. _google_activity_recognition_deltas: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L41
.. _google_activity_recognition_deltas.R: https://github.com/carissalow/rapids/blob/master/src/features/google_activity_recognition_deltas.R
.. _activity_features: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L74
.. _google_activity_recognition.py: https://github.com/carissalow/rapids/blob/master/src/features/google_activity_recognition.py
.. _`Light Config Code`: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L94
.. _light_features: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L113
.. _light_features.py: https://github.com/carissalow/rapids/blob/master/src/features/light_features.py
.. _`Location (Barnett’s) Config Code`: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L70
.. _phone_sensed_bins: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/preprocessing.snakefile#L46
.. _phone_sensed_bins.R: https://github.com/carissalow/rapids/blob/master/src/data/phone_sensed_bins.R
.. _resample_fused_location: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/preprocessing.snakefile#L67
.. _resample_fused_location.R: https://github.com/carissalow/rapids/blob/master/src/data/resample_fused_location.R
.. _location_barnett_features: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L49
.. _location_barnett_features.R: https://github.com/carissalow/rapids/blob/master/src/features/location_barnett_features.R
.. _`Screen Config Code`: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L88
.. _screen_deltas: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L33
.. _screen_deltas.R: https://github.com/carissalow/rapids/blob/master/src/features/screen_deltas.R
.. _screen_features: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L97
.. _screen_features.py: https://github.com/carissalow/rapids/blob/master/src/features/screen_features.py
.. _fitbit_with_datetime: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/preprocessing.snakefile#L94
.. _fitbit_readable_datetime.py: https://github.com/carissalow/rapids/blob/master/src/data/fitbit_readable_datetime.py
.. _`Fitbit: Sleep Config Code`: https://github.com/carissalow/rapids/blob/e952e27350c7ae02703bd444e8f92979e37d9ba6/config.yaml#L129
.. _fitbit_sleep_features: https://github.com/carissalow/rapids/blob/e952e27350c7ae02703bd444e8f92979e37d9ba6/rules/features.snakefile#L209
.. _fitbit_sleep_features.py: https://github.com/carissalow/rapids/blob/master/src/features/fitbit_sleep_features.py
.. _`Fitbit: Heart Rate Config Code`: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L113
.. _fitbit_heartrate_features: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L151
.. _fitbit_heartrate_features.py: https://github.com/carissalow/rapids/blob/master/src/features/fitbit_heartrate_features.py
.. _`Fitbit: Steps Config Code`: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L117
.. _fitbit_step_features: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L162
.. _fitbit_step_features.py: https://github.com/carissalow/rapids/blob/master/src/features/fitbit_step_features.py
.. _`Fitbit documentation`: https://help.fitbit.com/articles/en_US/Help_article/1565
.. _`Custom Catalogue File`: https://github.com/carissalow/rapids/blob/master/data/external/stachl_application_genre_catalogue.csv
.. _top1global: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L108
.. _`Beiwe Summary Statistics`: http://wiki.beiwe.org/wiki/Summary_Statistics
.. _`Pause-Flight Model`: https://academic.oup.com/biostatistics/advance-article/doi/10.1093/biostatistics/kxy059/5145908
