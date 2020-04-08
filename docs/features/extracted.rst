.. _rapids_metrics:

RAPIDS Metrics
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

- ``DAY_SEGMENTS`` - The list of day epochs that metrics can be segmented into: ``daily``, ``morning`` (6am-12pm), ``afternnon`` (12pm-6pm), ``evening`` (6pm-12am) and ``night`` (12am-6am). This list can be modified globally or on a per sensor basis. See DAY_SEGMENTS_ in ``config`` file.

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

    On any given day, Aware could have sensed data only for a few minutes or for 24 hours. Daily estimates of metrics should be considered more reliable the more hours Aware was running and logging data (for example, 10 calls logged on a day when only one hour of data was recorded is a less reliable measurement compared to 10 calls on a day when 23 hours of data were recorded. 

    Therefore, we define a valid hour as those that contain at least a certain number of valid bins. In turn, a valid bin are those that contain at least one row of data from any sensor logged within that period. We divide an hour into N bins of size ``BIN_SIZE`` (in minutes) and we mark an hour as valid if contains at least ``MIN_BINS_PER_HOUR`` of valid bins (out of the total possible number of bins that can be captured in an hour i.e. out of 60min/``BIN_SIZE`` bins). Days with valid sensed hours less than ``MIN_VALID_HOURS`` will be excluded form the output of this file. See PHONE_VALID_SENSED_DAYS_ in ``config.yaml``.

    In RAPIDS, you will find that we use ``phone_sensed_bins`` (a list of all valid and invalid bins of all monitored days) to improve the estimation of metrics that are ratios over time periods like ``episodepersensedminutes`` of :ref:`Screen<screen-sensor-doc>`.


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

- Extract SMS metrics:

      | ``expand("data/processed/{pid}/sms_{sms_type}_{day_segment}.csv".``
      |                     ``pid=config["PIDS"],``
      |                     ``sms_type = config["SMS"]["TYPES"],``
      |                     ``day_segment = config["SMS"]["DAY_SEGMENTS"]),``

**Rule Chain:**

- **Rule:** ``rules/preprocessing.snakefile/download_dataset`` - See the download_dataset_ rule.

    - **Script:** ``src/data/download_dataset.R`` - See the download_dataset.R_ script.
    
- **Rule:** ``rules/preprocessing.snakefile/readable_datetime`` - See the readable_datetime_ rule.

    - **Script:** ``src/data/readable_datetime.R`` - See the readable_datetime.R_ script.

- **Rule:** ``rules/features.snakefile/sms_metrics`` - See the sms_metric_ rule.

    - **Script:** ``src/features/sms_metrics.R`` - See the sms_metrics.R_ script.


.. _sms-parameters:

**SMS Rule Parameters:**

============    ===================
Name	        Description
============    ===================
sms_type        The particular ``sms_type`` that will be analyzed. The options for this parameter are ``received`` or ``sent``.
day_segment     The particular ``day_segments`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
metrics         The different measures that can be retrieved from the dataset. These metrics are available for both ``sent`` and ``received`` SMS messages. See :ref:`Available SMS Metrics <sms-available-metrics>` Table below
============    ===================

.. _sms-available-metrics:

**Available SMS Metrics**

The following table shows a list of the available metrics for both ``sent`` and ``received`` SMS. 

=========================   =========     =============
Name                        Units         Description
=========================   =========     =============
count                       SMS           A count of the number of times that particular ``sms_type`` occurred for a particular ``day_segment``.
distinctcontacts            contacts      A count of distinct contacts that were communicated for a particular ``sms_type`` for a particular ``day_segment``.
timefirstsms                minutes       The time in minutes from 12:00am (Midnight) that the first of a particular ``sms_type`` occurred.
timelastsms                 minutes       The time in minutes from 12:00am (Midnight) that the last of a particular ``sms_type`` occurred.
countmostfrequentcontact    SMS           The count of the number of sms messages of a particular``sms_type`` for the most contacted contact for a particular ``day_segment``.
=========================   =========     =============

**Assumptions/Observations:** 

    #. ``TYPES`` and ``METRICS`` keys need to match. From example::

        SMS:
            TYPES: [sent]
            METRICS: 
                sent: [count, distinctcontacts, timefirstsms, timelastsms, countmostfrequentcontact]

In the above config setting code the ``TYPE`` ``sent`` matches the ``METRICS`` key ``sent``.


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
day_segment     The particular ``day_segments`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
features         The different measures that can be retrieved from the calls dataset. Note that the same features are available for both ``incoming`` and ``outgoing`` calls, while ``missed`` calls has its own set of features. See :ref:`Available Incoming and Outgoing Call Features <available-in-and-out-call-features>` Table and :ref:`Available Missed Call Features <available-missed-call-features>` Table below.
============    ===================

.. _available-in-and-out-call-features:

**Available Incoming and Outgoing Call Features**

The following table shows a list of the available features for ``incoming`` and ``outgoing`` calls. 

=========================   =========     =============
Name                        Units         Description
=========================   =========     =============
count                       calls         A count of the number of times that a particular ``call_type`` occurred for a particular ``day_segment``.
distinctcontacts            contacts      A count of distinct contacts that were communicated with for a particular ``call_type`` for a particular ``day_segment`` 
meanduration                minutes       The mean duration of all calls for a particular ``call_type`` and ``day_segment``.
sumduration                 minutes       The sum of the duration of all calls for a particular ``call_type`` and ``day_segment``.
minduration                 minutes       The duration of the shortest call for a particular ``call_type`` and ``day_segment``.
maxduration                 minutes       The duration of the longest call for a particular ``call_type`` and ``day_segment``.
stdduration                 minutes       The standard deviation of all the calls for a particular ``call_type`` and ``day_segment``.
modeduration                minutes       The mode duration of all the calls for a particular ``call_type`` and ``day_segment``.
hubermduration                            The generalized Huber M-estimator of location of the MAD for the durations of all the calls for a particular ``call_type`` and ``day_segment``.
varqnduration                             The Location-Free Scale Estimator Qn of the durations of all the calls for a particular ``call_type`` and ``day_segment``.
entropyduration                           The estimate of the Shannon entropy H of the durations of all the calls for a particular ``call_type`` and ``day_segment``.
timefirstcall               minutes       The time in minutes from 12:00am (Midnight) that the first of ``call_type`` occurred.
timelastcall                minutes       The time in minutes from 12:00am (Midnight) that the last of ``call_type`` occurred.
countmostfrequentcontact    calls         The count of the number of calls of a particular ``call_type`` and ``day_segment`` for the most contacted contact.
=========================   =========     =============

.. _available-missed-call-features:

**Available Missed Call Features**

The following table shows a list of the available features for ``missed`` calls. 

=========================   =========     =============
Name                        Units         Description
=========================   =========     =============
count                       calls         A count of the number of times a ``missed`` call occurred for a particular ``day_segment``.
distinctcontacts            contacts      A count of distinct contacts whose calls were ``missed``.
timefirstcall               minutes       The time in minutes from 12:00am (Midnight) that the first ``missed`` call occurred.
timelastcall                minutes       The time in minutes from 12:00am (Midnight) that the last ``missed`` call occurred.
countmostfrequentcontact    CALLS           The count of the number of ``missed`` calls for the contact with the most ``missed`` calls.
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
day_segment     The particular ``day_segments`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
features         The different measures that can be retrieved from the Bluetooth dataset. See :ref:`Available Bluetooth Features <bluetooth-available-features>` Table below
============    ===================

.. _bluetooth-available-features:

**Available Bluetooth Features**

The following table shows a list of the available features for Bluetooth. 

===========================   =========     =============
Name                          Units         Description
===========================   =========     =============
countscans                    scans         Count of scans (a scan is a row containing a single Bluetooth device detected by Aware)
uniquedevices                 devices       Unique devices (number of unique devices identified by their hardware address -bt_address field)
countscansmostuniquedevice    scans         Count of scans of the most unique device across each participant’s dataset
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
day_segment     The particular ``day_segments`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
features         The different measures that can be retrieved from the dataset. See :ref:`Available Accelerometer Features <accelerometer-available-features>` Table below
============    ===================

.. _accelerometer-available-features:

**Available Accelerometer Features**

The following table shows a list of the available features the accelerometer sensor data for a particular ``day_segment``. 

====================================   ==============    =============
Name                                   Units             Description
====================================   ==============    =============
maxmagnitude                           m/s\ :sup:`2`      The maximum magnitude of acceleration (:math:`\|acceleration\| = \sqrt{x^2 + y^2 + z^2}`).
minmagnitude                           m/s\ :sup:`2`     The minimum magnitude of acceleration.
avgmagnitude                           m/s\ :sup:`2`     The average magnitude of acceleration.
medianmagnitude                        m/s\ :sup:`2`     The median magnitude of acceleration.
stdmagnitude                           m/s\ :sup:`2`     The standard deviation of acceleration.
ratioexertionalactivityepisodes                          The ratio of exertional activity time periods to total time periods.
sumexertionalactivityepisodes          minutes           The total minutes of performing exertional activity during the epoch
longestexertionalactivityepisode       minutes           The longest episode of performing exertional activity
longestnonexertionalactivityepisode    minutes           The longest episode of performing non-exertional activity
countexertionalactivityepisodes        episodes          The count of the episodes of performing exertional activity
countnonexertionalactivityepisodes     episodes          The count of the episodes of performing non-exertional activity
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
- iOS

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
day_segment             The particular ``day_segments`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
single_categories       A single category of apps that will be included  for the data collection. The available categories can be defined in the ``APPLICATION_GENRES`` in the ``config`` file. See :ref:`Assumtions and Observations <applications-foreground-observations>`.
multiple_categories     Categories of apps that will be included  for the data collection. The available categories can be defined in the ``APPLICATION_GENRES`` in the ``config`` file. See :ref:`Assumtions and Observations <applications-foreground-observations>`. 
single_apps             Any Android app can be included in the list of apps used to collect data by adding the package name to this list. (E.g. Youtube)
excluded_categories     Categories of apps that will be excluded for the data collection. The available categories can be defined in the ``APPLICATION_GENRES`` in the ``config`` file. See :ref:`Assumtions and Observations <applications-foreground-observations>`. 
excluded_apps           Any Android app can be excluded from the list of apps used to collect data by adding the package name to this list.
features                 The different measures that can be retrieved from the dataset. See :ref:`Available Applications Foreground Features <applications-foreground-available-features>` Table below
====================    ===================

.. _applications-foreground-available-features:

**Available Applications Foreground Features**

The following table shows a list of the available features for the Applications Foreground dataset 

==================   =========   =============
Name                 Units       Description
==================   =========   =============
count                apps        A count number of times using ``all_apps``, ``single_app``, ``single_category`` apps or ``multiple_category`` apps.
timeoffirstuse       contacts    The time in minutes from 12:00am (Midnight) to first use of any app (i.e. ``all_apps``), ``single_app``, ``single_category`` apps or ``multiple_category`` apps.
timeoflastuse        minutes     The time in minutes from 12:00am (Midnight) to the last of use of any app (i.e. ``all_apps``), ``single_app``, ``single_category`` apps or ``multiple_category`` apps.
frequencyentropy     shannons    The entropy of the apps frequency for ``all_apps``, ``single_category`` apps or ``multiple_category`` apps. There is no entropy for ``single_app`` apos.
==================   =========   =============

.. _applications-foreground-observations:

**Assumptions/Observations:** 

The ``APPLICATION_GENRES`` configuration (See `Application Genres Config`_ setting defines that catalogue of categories of apps that available for the pipeline. The ``CATALOGUE_SOURCE`` defines the source of the catalogue which can be ``FILE`` i.e. a custom file like the file provided with this project (See `Custom Catalogue File`_) or ``GOOGLE`` which is category classifications provided by Google. The ``CATALOGUE_FILE`` variable defines the path to the location of the custom file that contains the custom app catalogue. If ``CATALOGUE_SOURCE`` is equal to ``FILE``, the ``UPDATE_CATALOGUE_FILE`` variable specifies (``TRUE`` or ``FALSE``) whether or not to update ``CATALOGUE_FILE``, if ``CATALOGUE_SOURCE`` is equal to ``GOOGLE`` all scraped genres will be saved to ``CATALOGUE_FILE``. The ``SCRAPE_MISSING_GENRES`` is a ``TRUE`` or ``FALSE`` variable that specifies whether or not to scrape missing genres, only effective if ``CATALOGUE_SOURCE`` is equal to ``FILE``. If ``CATALOGUE_SOURCE`` is equal to ``GOOGLE``, all genres are scraped anyway. It should be noted that the ``top1global`` option finds and uses the most used app for that participant for the study. 



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

- Extract Battery Metrics:

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
    
- **Rule:** ``rules/features.snakefile/battery_metrics`` - See the battery_metrics_ rule

    - **Script:** ``src/features/battery_metrics.py`` - See the battery_metrics.py_ script.
    
.. _battery-parameters:

**Battery Rule Parameters:**

============    ===================
Name	        Description
============    ===================
day_segment     The particular ``day_segments`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
metrics         The different measures that can be retrieved from the Battery dataset. See :ref:`Available Battery Metrics <battery-available-metrics>` Table below
============    ===================

.. _battery-available-metrics:

**Available Battery Metrics**

The following table shows a list of the available metrics for Battery data. 

=====================   ===============   =============
Name                    Units             Description
=====================   ===============   =============
countdischarge          episodes          A count of the number of battery discharging episodes
sumdurationdischarge    hours             The total duration of all discharging episodes (time the phone was discharging)
countcharge             episodes          A count of the number of battery charging episodes
sumdurationcharge       hours             The total duration of all charging episodes (time the phone was charging)
avgconsumptionrate      episodes/hours    The average of the ratios between discharging episodes’ battery delta and duration
maxconsumptionrate      episodes/hours    The maximum of the ratios between discharging episodes’ battery delta and duration
=====================   ===============   =============

**Assumptions/Observations:** 


.. _google-activity-recognition-sensor-doc:

Google Activity Recognition
""""""""""""""""""""""""""""

See `Google Activity Recognition Config Code`_

**Available Epochs:**      

- daily 
- morning
- afternoon
- evening
- night

**Available Platforms:**    

- Android

**Snakefile entry:**

..  - Download raw Google Activity Recognition dataset: ``expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``

..  - Apply readable dateime to Google Activity Recognition dataset: ``expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``
    
..  - Extract the deltas in Google Activity Recognition dataset: ``expand("data/processed/{pid}/plugin_google_activity_recognition_deltas.csv", pid=config["PIDS"]),``
    
- Extract Sensor Metrics:

    | ``expand("data/processed/{pid}/google_activity_recognition_{segment}.csv",pid=config["PIDS"],``
    |                ``segment = config["GOOGLE_ACTIVITY_RECOGNITION"]["DAY_SEGMENTS"]),``
    
**Rule Chain:**

- **Rule:** ``rules/preprocessing.snakefile/download_dataset`` - See the download_dataset_ rule.

    - **Script:** ``src/data/download_dataset.R`` - See the download_dataset.R_ script.

- **Rule:** ``rules/preprocessing.snakefile/readable_datetime`` - See the readable_datetime_ rule.

    - **Script:** ``src/data/readable_datetime.R`` - See the readable_datetime.R_ script.

- **Rule:** ``rules/features.snakefile/google_activity_recognition_deltas`` - See the google_activity_recognition_deltas_ rule.

    - **Script:** ``src/features/google_activity_recognition_deltas.R`` - See the google_activity_recognition_deltas.R_ script.

- **Rule:** ``rules/features.snakefile/activity_metrics`` - See the activity_metrics_ rule.

    - **Script:** ``ssrc/features/google_activity_recognition.py`` - See the google_activity_recognition.py_ script.
    
.. _google-activity-recognition-parameters:

**Google Activity Recognition Rule Parameters:**

============    ===================
Name	        Description
============    ===================
day_segment     The particular ``day_segments`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
metrics         The different measures that can be retrieved from the Google Activity Recognition dataset. See :ref:`Available Google Activity Recognition Metrics <google-activity-recognition-available-metrics>` Table below
============    ===================

.. _google-activity-recognition-available-metrics:

**Available Google Activity Recognition Metrics**

The following table shows a list of the available metrics for the Google Activity Recognition dataset. 

======================   ============    =============
Name                     Units           Description
======================   ============    =============
count                    rows            A count of the number of rows of registered activities.
mostcommonactivity                       The most common activity.
countuniqueactivities    activities       A count of the number of unique activities.
activitychangecount      transitions     A count of any transition between two different activities, sitting to running for example.
sumstationary            minutes         The total duration of episodes of still and tilting (phone) activities.
summobile                minutes         The total duration of episodes of on foot, running, and on bicycle activities
sumvehicle               minutes         The total duration of episodes of on vehicle activity
======================   ============    =============

**Assumptions/Observations:** N/A

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
    
- Extract Light Metrics:

    | ``expand("data/processed/{pid}/light_{day_segment}.csv",``
    |                      ``pid=config["PIDS"],`` 
    |                      ``day_segment = config["LIGHT"]["DAY_SEGMENTS"]),``
    
**Rule Chain:**

- **Rule:** ``rules/preprocessing.snakefile/download_dataset`` - See the download_dataset_ rule.

    - **Script:** ``src/data/download_dataset.R`` - See the download_dataset.R_ script.

- **Rule:** ``rules/preprocessing.snakefile/readable_datetime`` - See the readable_datetime_ rule.

    - **Script:** ``src/data/readable_datetime.R`` - See the readable_datetime.R_ script.

- **Rule:** ``rules/features.snakefile/light_metrics`` - See the light_metrics_ rule.

    - **Script:** ``src/features/light_metrics.py`` - See the light_metrics.py_ script.

.. _light-parameters:

**Light Rule Parameters:**

============    ===================
Name	        Description
============    ===================
day_segment     The particular ``day_segments`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
metrics         The different measures that can be retrieved from the Light dataset. See :ref:`Available Light Metrics <light-available-metrics>` Table below
============    ===================

.. _light-available-metrics:

**Available Light Metrics**

The following table shows a list of the available metrics for the Light dataset. 

===========   =========     =============
Name          Units         Description
===========   =========     =============
count         rows          A count of the number of rows that light sensor recorded.
maxlux        lux           The maximum ambient luminance in lux units
minlux        lux           The minimum ambient luminance in lux units
avglux        lux           The average ambient luminance in lux units
medianlux     lux           The median ambient luminance in lux units
stdlux        lux           The standard deviation of ambient luminance in lux units
===========   =========     =============

**Assumptions/Observations:** N/A


.. _location-sensor-doc:

Location (Barnett’s) Features
""""""""""""""""""""""""""""""
Barnett’s location features are based on the concept of flights and pauses. GPS coordinates are converted into a 
sequence of flights (straight line movements) and pauses (time spent stationary). Data is imputed before metrics 
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

- Extract Sensor Metrics: ``expand("data/processed/{pid}/location_barnett.csv", pid=config["PIDS"]),``
    
**Rule Chain:**

- **Rule:** ``rules/preprocessing.snakefile/download_dataset`` - See the download_dataset_ rule.

    - **Script:** ``src/data/download_dataset.R`` - See the download_dataset.R_ script.

- **Rule:** ``rules/preprocessing.snakefile/readable_datetime`` - See the readable_datetime_ rule.

    - **Script:** ``src/data/readable_datetime.R`` - See the readable_datetime.R_ script.

- **Rule:** ``rules/preprocessing.snakefile/phone_sensed_bins`` - See the phone_sensed_bins_ rule.

    - **Script:** ``src/data/phone_sensed_bins.R`` - See the phone_sensed_bins.R_ script.

- **Rule:** ``rules/preprocessing.snakefile/resample_fused_location`` - See the resample_fused_location_ rule.

    - **Script:** ``src/data/resample_fused_location.R`` - See the resample_fused_location.R_ script.

- **Rule:** ``rules/features.snakefile/location_barnett_metrics`` - See the location_barnett_metrics_ rule.

    - **Script:** ``src/features/location_barnett_metrics.R`` - See the location_barnett_metrics.R_ script.

    
.. _location-parameters:

**Location Rule Parameters:**

=================    ===================
Name	             Description
=================    ===================
location_to_use      The specifies which of the location data will be use in the analysis. Possible options are ``ALL``, ``ALL_EXCEPT_FUSED`` OR ``RESAMPLE_FUSED``
accuracy_limit       This is in meters. The sensor drops location coordinates with an accuracy higher than this. This number means there's a 68% probability the true location is within this radius specified.
timezone             The timezone used to calculate location. 
metrics              The different measures that can be retrieved from the Location dataset. See :ref:`Available Location Metrics <location-available-metrics>` Table below
=================    ===================

.. _location-available-metrics:

**Available Location Metrics**

The following table shows a list of the available metrics for Location dataset. 

================   =========     =============
Name               Units         Description
================   =========     =============
hometime           minutes       Time at home. Time spent at home in minutes. Home is the most visited significant location between 8 pm and 8 am including any pauses within a 200-meter radius.
disttravelled      meters        Distance travelled. This is total distance travelled over a day.
rog                meters        The Radius of Gyration (RoG). It is a measure in meters of the area covered by a person over a day. A centroid is calculated for all the places (pauses) visited during a day and a weighted distance between all the places and the centroid is computed. The weights are proportional to the time spent in each place.
maxdiam            meters        The Maximum diameter. The largest distance in meters between any two pauses.
maxhomedist        meters        Max home distance. The maximum distance from home in meters.
siglocsvisited     locations     Significant locations. The number of significant locations visited during the day. Significant locations are computed using k-means clustering over pauses found in the whole monitoring period. The number of clusters is found iterating from 1 to 200 stopping until the centroids of two significant locations are within 400 meters of one another.
avgflightlen       meters        Avg flight length. Mean length of all flights
stdflightlen       meters        Std flight length. The standard deviation of the length of all flights.
avgflightdur       meters        Avg flight duration. Mean duration of all flights.
stdflightdur       meters        Std flight duration. The standard deviation of the duration of all flights.
probpause                        Pause probability. The fraction of a day spent in a pause (as opposed to a flight)
siglocentropy                    Significant location entropy. Entropy measurement based on the proportion of time spent at each significant location visited during a day.
minsmissing                            
circdnrtn           	         Circadian routine. A continuous metric that can take any value between 0 and 1, where 0 represents a daily routine completely different from any other sensed days and 1 a routine the same as every other sensed day.
wkenddayrtn        Weekend       circadian routine. Same as Circadian routine but computed separately for weekends and weekdays.
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

===============    ===================
Name	           Description
===============    ===================
day_segment        The particular ``day_segments`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
features_events     The different measures that can be retrieved from the events in the Screen dataset. See :ref:`Available Screen Events Features <screen-events-available-features>` Table below
features_deltas     The different measures that can be retrieved from the episodes extracted from the Screen dataset. See :ref:`Available Screen Episodes Features <screen-episodes-available-features>` Table below
episodes           The action that defines an episode
===============    ===================

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

=============   =========    =============
Name            Units        Description
=============   =========    =============
sumduration     seconds      Sum duration unlock: The sum duration of unlock episodes 
maxduration     seconds      Max duration unlock: The maximum duration of unlock episodes
minduration     seconds      Min duration unlock: The minimum duration of unlock episodes
avgduration     seconds      Average duration unlock: The average duration of unlock episodes
stdduration     seconds      Std duration unlock: The standard deviation of the duration of unlock episodes
=============   =========    =============

**Assumptions/Observations:** 

An ``unlock`` episode is considered as the time between an ``unlock`` event and a ``lock`` event. iOS recorded these episodes reliable (albeit some duplicated ``lock`` events within milliseconds from each other). However, in Android there are some events unrelated to the screen state because of multiple consecutive ``unlock``/``lock`` events, so we keep the closest pair. In the experiments these are less than 10% of the screen events collected. This happens because ``ACTION_SCREEN_OFF`` and ``ON`` are "sent when the device becomes non-interactive which may have nothing to do with the screen turning off". Additionally in Android it is possible to measure the time spent on the ``lock`` screen onto the ``unlock`` event and the total screen time (i.e. ``ON`` to ``OFF``) events but we are only keeping ``unlock`` episodes (``unlock`` to ``OFF``) to be consistent with iOS. 

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
      
- Extract Sensor Metrics:

    | ``expand("data/processed/{pid}/fitbit_heartrate_{day_segment}.csv",``
    |                      ``pid=config["PIDS"],`` 
    |                      ``day_segment = config["HEARTRATE"]["DAY_SEGMENTS"]),``
    
**Rule Chain:**

- **Rule:** ``rules/preprocessing.snakefile/download_dataset`` - See the download_dataset_ rule.

    - **Script:** ``src/data/download_dataset.R`` - See the download_dataset.R_ script.

- **Rule:** ``rules/preprocessing.snakefile/fitbit_with_datetime`` - See the fitbit_with_datetime_ rule.

    - **Script:** ``src/data/fitbit_readable_datetime.py`` - See the fitbit_readable_datetime.py_ script.

- **Rule:** ``rules/features.snakefile/fitbit_heartrate_metrics`` - See the fitbit_heartrate_metrics_ rule.

    - **Script:** ``src/features/fitbit_heartrate_metrics.py`` - See the fitbit_heartrate_metrics.py_ script.

    
.. _fitbit-heart-rate-parameters:

**Fitbit: Heart Rate Rule Parameters:**

============    ===================
Name	        Description
============    ===================
day_segment     The particular ``day_segments`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
metrics         The different measures that can be retrieved from the Fitbit: Heart Rate dataset. 
                See :ref:`Available Fitbit: Heart Rate Metrics <fitbit-heart-rate-available-metrics>` Table below
============    ===================

.. _fitbit-heart-rate-available-metrics:

**Available Fitbit: Heart Rate Metrics**

The following table shows a list of the available metrics for the Fitbit: Heart Rate dataset. 

==================   ===========    =============
Name                 Units          Description
==================   ===========    =============
maxhr                beats/mins     The maximum heart rate.
minhr                beats/mins     The minimum heart rate.
avghr                beats/mins     The average heart rate.
medianhr             beats/mins     The median heart rate.
modehr               beats/mins     The mode heart rate.
stdhr                beats/mins     The standard deviation of heart rate.
diffmaxmodehr        beats/mins     Diff max mode heart rate: The maximum heart rate minus mode heart rate.
diffminmodehr        beats/mins     Diff min mode heart rate: The mode heart rate minus minimum heart rate.
entropyhr                           Entropy heart rate: The entropy of heart rate.
lengthoutofrange     minutes        Length out of range: The duration of time the heart rate is in the ``out_of_range`` zone in minute.
lengthfatburn        minutes        Length fat burn: The duration of time the heart rate is in the ``fat_burn`` zone in minute.
lengthcardio         minutes        Length cardio: The duration of time the heart rate is in the ``cardio`` zone in minute.
lengthpeak           minutes        Length peak: The duration of time the heart rate is in the ``peak`` zone in minute
==================   ===========    =============

**Assumptions/Observations:** Heart rate zones contain 4 zones: ``out_of_range`` zone, ``fat_burn`` zone, ``cardio`` zone, and ``peak`` zone. Please refer to the `Fitbit documentation`_ for detailed information of how to define those zones.

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
day_segment                The particular ``day_segments`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, ``evening``, ``night``
metrics                    The different measures that can be retrieved from the dataset. See :ref:`Available Fitbit: Steps Metrics <fitbit-steps-available-metrics>` Table below
threshold_active_bout      The maximum number of steps per minute necessary for a bout to be ``sedentary``. That is, if the step count per minute is greater than this value the bout has a status of ``active``. 
=======================    ===================

.. _fitbit-steps-available-metrics:

**Available Fitbit: Steps Metrics**

The following table shows a list of the available metrics for the Fitbit: Steps dataset. 

=========================   =========     =============
Name                        Units         Description
=========================   =========     =============
sumallsteps                 steps         Sum all steps: The total step count.
maxallsteps                 steps         Max all steps: The maximum step count
minallsteps                 steps         Min all steps: The minimum step count
avgallsteps                 steps         Avg all steps: The average step count
stdallsteps                 steps         Std all steps: The standard deviation of step count
countsedentarybout          bouts         Count sedentary bout: A count of sedentary bouts
maxdurationsedentarybout    minutes       Max duration sedentary bout: The maximum duration of sedentary bouts
mindurationsedentarybout    minutes       Min duration sedentary bout: The minimum duration of sedentary bouts
avgdurationsedentarybout    minutes       Avg duration sedentary bout: The average duration of sedentary bouts
stddurationsedentarybout    minutes       Std duration sedentary bout: The standard deviation of the duration of sedentary bouts
countactivebout             bouts         Count active bout: A count of active bouts
maxdurationactivebout       minutes       Max duration active bout: The maximum duration of active bouts
mindurationactivebout       minutes       Min duration active bout: The minimum duration of active bouts
avgdurationactivebout       minutes       Avg duration active bout: The average duration of active bouts
stddurationactivebout       minutes       Std duration active bout: The standard deviation of the duration of active bouts
=========================   =========     =============

**Assumptions/Observations:** If the step count per minute smaller than the ``THRESHOLD_ACTIVE_BOUT`` (default value is 10), it is defined as sedentary status. Otherwise, it is defined as active status. One active/sedentary bout is a period during with the user is under ``active``/``sedentary`` status.
	

.. -------------------------Links ------------------------------------ ..

.. _SENSORS: https://github.com/carissalow/rapids/blob/f22d1834ee24ab3bcbf051bc3cc663903d822084/config.yaml#L2
.. _`SMS Config Code`: https://github.com/carissalow/rapids/blob/f22d1834ee24ab3bcbf051bc3cc663903d822084/config.yaml#L38
.. _AWARE: https://awareframework.com/what-is-aware/
.. _`List of Timezones`: https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
.. _sms_metric: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L1
.. _sms_metrics.R: https://github.com/carissalow/rapids/blob/master/src/features/sms_metrics.R
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
.. _accelerometer_metrics: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L124
.. _accelerometer_metrics.py: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/src/features/accelerometer_metrics.py
.. _`Applications Foreground Config Code`: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L102
.. _`Application Genres Config`: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L54
.. _application_genres: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/preprocessing.snakefile#L81
.. _application_genres.R: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/src/data/application_genres.R
.. _applications_foreground_features: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L135
.. _applications_foreground_features.py: https://github.com/carissalow/rapids/blob/master/src/features/accelerometer_features.py
.. _`Battery Config Code`: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L84
.. _battery_deltas: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L25
.. _battery_deltas.R: https://github.com/carissalow/rapids/blob/master/src/features/battery_deltas.R
.. _battery_metrics: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L86
.. _battery_metrics.py : https://github.com/carissalow/rapids/blob/master/src/features/battery_metrics.py
.. _`Google Activity Recognition Config Code`: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L80
.. _google_activity_recognition_deltas: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L41
.. _google_activity_recognition_deltas.R: https://github.com/carissalow/rapids/blob/master/src/features/google_activity_recognition_deltas.R
.. _activity_metrics: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L74
.. _google_activity_recognition.py: https://github.com/carissalow/rapids/blob/master/src/features/google_activity_recognition.py
.. _`Light Config Code`: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L94
.. _light_metrics: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L113
.. _light_metrics.py: https://github.com/carissalow/rapids/blob/master/src/features/light_metrics.py
.. _`Location (Barnett’s) Config Code`: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L70
.. _phone_sensed_bins: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/preprocessing.snakefile#L46
.. _phone_sensed_bins.R: https://github.com/carissalow/rapids/blob/master/src/data/phone_sensed_bins.R
.. _resample_fused_location: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/preprocessing.snakefile#L67
.. _resample_fused_location.R: https://github.com/carissalow/rapids/blob/master/src/data/resample_fused_location.R
.. _location_barnett_metrics: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L49
.. _location_barnett_metrics.R: https://github.com/carissalow/rapids/blob/master/src/features/location_barnett_metrics.R
.. _`Screen Config Code`: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L88
.. _screen_deltas: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L33
.. _screen_deltas.R: https://github.com/carissalow/rapids/blob/master/src/features/screen_deltas.R
.. _screen_features: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L97
.. _screen_features.py: https://github.com/carissalow/rapids/blob/master/src/features/screen_features.py
.. _`Fitbit: Heart Rate Config Code`: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L113
.. _fitbit_with_datetime: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/preprocessing.snakefile#L94
.. _fitbit_readable_datetime.py: https://github.com/carissalow/rapids/blob/master/src/data/fitbit_readable_datetime.py
.. _fitbit_heartrate_metrics: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L151
.. _fitbit_heartrate_metrics.py: https://github.com/carissalow/rapids/blob/master/src/features/fitbit_heartrate_metrics.py
.. _`Fitbit: Steps Config Code`: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L117
.. _fitbit_step_features: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L162
.. _fitbit_step_features.py: https://github.com/carissalow/rapids/blob/master/src/features/fitbit_step_features.py
.. _`Fitbit documentation`: https://help.fitbit.com/articles/en_US/Help_article/1565
.. _`Custom Catalogue File`: https://github.com/carissalow/rapids/blob/master/data/external/stachl_application_genre_catalogue.csv
.. _top1global: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L108
.. _`Beiwe Summary Statistics`: http://wiki.beiwe.org/wiki/Summary_Statistics
.. _`Pause-Flight Model`: https://academic.oup.com/biostatistics/advance-article/doi/10.1093/biostatistics/kxy059/5145908