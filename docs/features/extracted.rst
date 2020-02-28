.. _rapids_metrics:

RAPIDS Metrics
===============

This following is documentation of on the RAPIDS metrics settings in the configuation file. 

.. _sensor-list:

    - ``SENSORS`` - This varable stores a list of the names of the sensor data that are being pulled from the AWARE_ database. These names are the actual names of the  tables that the data is found in the database. See SENSORS_ variable in ``config`` file.  

.. _fitbit-table:

    - ``FITBIT_TABLE`` - The name of the fitbit database 

.. _fitbit-sensors:

    - ``FITBIT_SENSORS`` - The list of sensors that to be pulled from the fitbit database

.. _pid: 

    - ``PID`` - The list of participant ids included in the analysis. Remember that you must create a file named ``pXXX`` for each participant in the ``data/external`` directory containing there device_id. (Remember installation :ref:`step 8 <install-step-8>`)

.. _day-segments: 

    - ``DAY_SEGMENTS`` - The list of common day segments (time frequency/checkpoints) that data would be analyzed. See DAY_SEGMENTS_ in ``config`` file.

.. _timezone:

    - ``TIMEZONE`` - The timezone of the server. Use the timezone names from this `List of Timezones`_. Double check your code, for example EST is not US Eastern Time.

.. _database_group:

    - ``DATABASE_GROUP`` - The name of the research project database. 

.. _download-dataset:

    - ``DOWNLOAD_DATASET`` - The name of the dataset for the research project. 

.. _readable-datetime:

    - ``READABLE_DATETIME`` - Readable datetime configuration. Defines the format that the readable date and time should be. 

.. _phone-valid-sensed-days:

    - ``PHONE_VALID_SENSED_DAYS`` - Specifies the ``BIN_SIZE``, ``MIN_VALID_HOURS``, ``MIN_BINS_PER_HOUR``. ``BIN_SIZE`` is the time that the data is aggregated. ``MIN_VALID_HOURS`` is the minimum numbers of hours data will be gathered within a 24 hour period (a day). Finally ``MIN_BINS_PER_HOUR`` specifies minimum number of bins that are captured per hour. This is out of the total possible number of bins that can be captured in an hour i.e. out of 60min/``BIN_SIZE`` bins. See PHONE_VALID_SENSED_DAYS_ in ``config`` file.


.. _individual-sensor-settings:

List of Indvidual Sensors and There Settings
---------------------------------------------

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

    - Download raw SMS dataset: ``expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``

    - Download raw SMS dataset with readable: ``expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``

    - Extract SMS metrics

      | ``expand("data/processed/{pid}/sms_{sms_type}_{day_segment}.csv".``
      |                     ``pid=config["PIDS"],``
      |                     ``sms_type = config["SMS"]["TYPES"],``
      |                     ``day_segment = config["SMS"]["DAY_SEGMENTS"]),``

**Rule Chain:**

    - ``rules/preprocessing.snakefile/download_dataset`` - See the download_dataset_ rule.

        - ``src/data/download_dataset.R`` - See the download_dataset.R_ script.
    
    - ``rules/preprocessing.snakefile/readable_datetime`` - See the readable_datetime_ rule.

        - ``src/data/readable_datetime.R`` - See the readable_datetime.R_ script.

    - ``rules/features.snakefile/sms_metrics`` - See the sms_metric_ rule.

        - ``src/features/sms_metrics.R`` - See the sms_metrics.R_ script.

.. _sms-parameters:

**SMS Rule Parameters:**

============    ===================
Name	        Description
============    ===================
sms_type        The particular ``sms_type`` that will be analyzed. The options for this parameter is ``received`` or ``sent``.
day_segment     The particular ``day_segments`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, 
                ``evening``, ``night``
metrics         The different measures that can be retrieved from the dataset. These metrics are available for both ``sent`` and ``received``
                SMS messages. See :ref:`Available SMS Metrices <sms-available-metrics>` Table below
============    ===================

.. _sms-available-metrics:

**Available SMS Metrics**

The following table shows a list of the available metrics for both ``sent`` and ``received`` SMS. 

=========================   =========     =============
Name                        Units         Description
=========================   =========     =============
count                       SMS           A count of the number of times that particular ``sms_type`` occured for a particular ``day_segment``.
distinctcontacts            contacts      A count of distinct contacts that were comunicated for a particular ``sms_type`` for a particular 
                                          ``day_segment``.
timefirstsms                minutes       The time in minutes from 12:00am (Midnight) that the first of a particular ``sms_type`` occured.
timelastsms                 minutes       The time in minutes from 12:00am (Midnight) that the last of a particular ``sms_type`` occured.
countmostfrequentcontact    SMS           The count of the number of sms meassages of a particular``sms_type`` for the most contacted contact for 
                                          a particular ``day_segment``.
=========================   =========     =============

Assumptions/Observations: 

    #. ``TYPES`` and ``METRICS`` keys need to match. From example::

        SMS:
            TYPES : [sent]
            METRICS: 
                sent: [count, distinctcontacts, timefirstsms, timelastsms, countmostfrequentcontact]

In the above config setting code the ``TYPE``  ``sent`` matches the ``METRICS`` key ``sent``.


.. _call-sensor-doc:

Calls
"""""""""""""

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

    - Download raw Calls dataset: ``expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``

    - Download raw Calls dataset with readable: ``expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``
    - Extract Calls Metrics
    
      | ``expand("data/processed/{pid}/call_{call_type}_{segment}.csv",``
      |                      ``pid=config["PIDS"],`` 
      |                      ``call_type=config["CALLS"]["TYPES"],``
      |                      ``segment = config["CALLS"]["DAY_SEGMENTS"]),``
    
**Rule Chain:**

    - ``rules/preprocessing.snakefile/download_dataset`` - See the download_dataset_ rule.

         - ``src/data/download_dataset.R`` - See the download_dataset.R_ script.
    
    - ``rules/preprocessing.snakefile/readable_datetime`` - See the readable_datetime_ rule.

        - ``src/data/readable_datetime.R`` - See the readable_datetime.R_ script.

    - ``rules/features.snakefile/call_metrics`` - See the call_metrics_ rule.

        - ``src/features/call_metrics.R`` - See the call_metrics.R_ script.

    
.. _calls-parameters:

**Sensor Rule Parameters:**

============    ===================
Name	        Description
============    ===================
call_type       The particular ``call_type`` that will be analyzed. The options for this parameter are ``incoming``, ``outgoing`` or ``missed``.
day_segment     The particular ``day_segment`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, 
                ``evening``, ``night``
metrics         The different measures that can be retrieved from the calls dataset. Note that the same metrics are available for both 
                ``incoming`` and ``outgoing`` calls,  while ``missed`` calls has its own set of metrics. See :ref:`Available Incoming and Outgoing Call Metrices <available-in-and-out-call-metrics>` Table and :ref:`Available Missed Call Metrices <available-missed-call-metrics>` Table below.
============    ===================

.. _available-in-and-out-call-metrics:

**Available Incoming and Outgoing Call Metrices**

The following table shows a list of the available metrics for ``incoming`` and ``outgoing`` calls. 

=========================   =========     =============
Name                        Units         Description
=========================   =========     =============
count                       calls         A count of the number of times that a particular ``call_type`` occured for a particular ``day_segment``.
distinctcontacts            contacts      A count of distinct contacts that were comunicated with for a particular ``call_type`` for a particular 
                                          ``day_segment`` 
meanduration                minutes       The mean duration of all calls for a particular ``call_type`` and ``day_segment``.
sumduration                 minutes       The sum of the duration of all calls for a particular ``call_type`` and ``day_segment``.
minduration                 minutes       The duration of the shortest call for a particular ``call_type`` and ``day_segment``.
maxduration                 minutes       The duration of the longest call for a particular ``call_type`` and ``day_segment``.
stdduration                 minutes       The standard deviation of all the calls for a particular ``call_type`` and ``day_segment``.
modeduration                minutes       The mode duration of all the calls for a particular ``call_type`` and ``day_segment``.
hubermduration                            The generalized Huber M-estimator of location of the MAD for the durations of all the calls for a 
                                          particular ``call_type`` and ``day_segment``.
varqnduration                             The location-Free Scale Estimator Qn of the durations of all the calls for a particular ``call_type`` 
                                          and ``day_segment``.
entropyduration                           The estimates the Shannon entropy H of the durations of all the calls for a particular ``call_type`` 
                                          and ``day_segment``.
timefirstcall               minutes       The time in minutes from 12:00am (Midnight) that the first of ``call_type`` occured.
timelastcall                minutes       The time in minutes from 12:00am (Midnight) that the last of ``call_type`` occured.
countmostfrequentcontact    calls         The count of the number of calls of a particular ``call_type`` and ``day_segment`` for the most contacted contact.
=========================   =========     =============

.. _available-missed-call-metrics:

**Available Missed Call Metrices**

The following table shows a list of the available metrics for ``missed`` calls. 

=========================   =========     =============
Name                        Units         Description
=========================   =========     =============
count                       calls         A count of the number of times a ``missed`` call occured for a particular ``day_segment``.
distinctcontacts            contacts      A count of distinct contacts whose calls were ``missed``.
timefirstcall               minutes       The time in minutes from 12:00am (Midnight) that the first ``missed`` call occured.
timelastcall                minutes       The time in minutes from 12:00am (Midnight) that the last ``missed`` call occured.
countmostfrequentcontact    SMS           The count of the number of ``missed`` calls for the contact with the most ``missed`` calls.
=========================   =========     =============

Assumptions/Observations: 

    #. ``TYPES`` and ``METRICS`` keys need to match. From example::

        SMS:
            TYPES : [missed]
            METRICS: 
                missed: [count, distinctcontacts, timefirstsms, timelastsms, countmostfrequentcontact]

In the above config setting code the ``TYPE``  ``missed`` matches the ``METRICS`` key ``missed``.


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

    - Download raw Bluetooth dataset: ``expand("data/raw/{pid}/{sensor}_raw.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``

    - Download raw Bluetooth dataset with readable: ``expand("data/raw/{pid}/{sensor}_with_datetime.csv", pid=config["PIDS"], sensor=config["SENSORS"]),``
    - Extract Bluetooth Metrics
    
      | ``expand("data/processed/{pid}/bluetooth_{segment}.csv",``
      |          ``pid=config["PIDS"],`` 
      |          ``segment = config["BLUETOOTH"]["DAY_SEGMENTS"]),``
    
**Rule Chain:**

    - ``rules/preprocessing.snakefile/download_dataset`` - See the download_dataset_ rule.

        - ``src/data/download_dataset.R`` See the download_dataset.R_ script.
    
    - ``rules/preprocessing.snakefile/readable_datetime`` - See the readable_datetime_ rule.

        - ``src/data/readable_datetime.R`` See the readable_datetime.R_ script.

    - ``rules/features.snakefile/bluetooth_metrics`` - See the bluetooth_metric_ rule.

        - ``src/features/bluetooth_metrics.R`` - See the bluetooth_metrics.R_ script.

    
.. _bluetooth-parameters:

**Bluetooth Rule Parameters:**

============    ===================
Name	        Description
============    ===================
day_segment     The particular ``day_segment`` that will be analyzed. The available options are ``daily``, ``morning``, ``afternoon``, 
                ``evening``, ``night``
metrics         The different measures that can be retrieved from the Bluetooth dataset. See :ref:`Available Bluetooth Metrices <bluetooth-available-metrics>` Table below
============    ===================

.. _bluetooth-available-metrics:

**Available Bluetooth Metrics**

The following table shows a list of the available metrics for Bluetooth. 

===========================   =========     =============
Name                          Units         Description
===========================   =========     =============
countscans                    scans         Count of scans (a scan is a row containing a single Bluetooth device detected by Aware)
uniquedevices                 devices       Unique devices (number of unique devices identified by their hardware address -bt_address field)
countscansmostuniquedevice    scans         Count of scans of the most unique device across each participant’s dataset
===========================   =========     =============

Assumptions/Observations: N/A 



.. _accelerometer:

Accelerometer
--------------

Available epochs: daily, morning, afternoon, evening, and night

- Max magnitude: maximum magnitude of acceleration (:math:`\|acceleration\| = \sqrt{x^2 + y^2 + z^2}`)
- Min magnitude: minimum magnitude of acceleration
- Avg magnitude: average magnitude of acceleration
- Median magnitude: median magnitude of acceleration
- Std magnitude: standard deviation of acceleration
- Ratio exertional activity episodes: ratio of exertional activity time periods to total time periods
- Sum exertional activity episodes: total minutes of performing exertional activity during the epoch
- Longest exertional activity episode: longest episode of performing exertional activity
- Longest non-exertional activity episode: longest episode of performing non-exertional activity
- Count exertional activity episodes: count of the episods of performing exertional activity
- Count non-exertional activity episodes: count of the episodes of performing non-exertional activity

.. _applications_foreground:

Applications_foreground
-------------------------

Available epochs: daily, morning, afternoon, evening, and night

- Count: number of times using all_apps/single_app/single_category/multiple_category
- Time of first use: time of first use all_apps/single_app/single_category/multiple_category in minutes
- Time of last use: time of last use all_apps/single_app/single_category/multiple_category in minutes
- Frenquency entropy: the entropy of the apps frequency for all_apps/single_app/single_category/multiple_category. There is no entropy for single_app.

.. _battery:

Battery
--------

Available epochs: daily, morning, afternoon, evening, and night

-	Count discharge: number of battery discharging episodes
-	Sum duration discharge: total duration of all discharging episodes (time the phone was discharging)
-	Average consumption rate: average of the ratios between discharging episodes’ battery delta and duration
-	Max consumption rate: max of the ratios between discharging episodes’ battery delta and duration
-	Count charge: number of battery charging episodes
-	Sum duration charge: total duration of all charging episodes (time the phone was charging)

.. _google-activity-recognition:

Google Activity Recognition
---------------------------

Available epochs: daily, morning, afternoon, evening, and night

-	Count (number of rows)
-	Most common activity
-	Number of unique activities
-	Activity change count (any transition between two different activities, sitting to running for example)
-	Sum stationary: total duration of episodes of still and tilting (phone) activities
-	Sum mobile: total duration of episodes of on foot, running, and on bicycle activities
-	Sum vehicle: total duration of episodes of on vehicle activity

.. _light:

Light
-----

Available epochs: daily, morning, afternoon, evening, and night

- Count (number of rows)
- Max lux: maximum ambient luminance in lux units
- Min lux: minimum ambient luminance in lux units
- Avg lux: average ambient luminance in lux units
- median lux: median ambient luminance in lux units
- Std lux: standard deviation of ambient luminance in lux units

.. _location-features:

Location (Barnett’s) Features
-----------------------------

Available epochs: daily

Barnett’s location features are based on the concept of flights and pauses. GPS coordinates are converted into a sequence of flights (straight line movements) and pauses (time spent stationary). Data is imputed before metrics are computed (https://arxiv.org/abs/1606.06328)

-	Time at home. Time spent at home in minutes. Home is the most visited significant location between 8 pm and 8 am including any pauses within a 200-meter radius.
-	Max home distance. Maximum distance from home in meters.
-	Pause probability. The fraction of a day spent in a pause (as opposed to a flight)
-	Circadian routine. A continuous metric that can take any value between 0 and 1, where 0 represents a daily routine completely different from any other sensed days and 1 a routine the same as every other sensed day.
-	Wkn circadian routine. Same as Circadian routine but computed separately for weekends and weekdays.
-	Distance travelled. Total distance travelled over a day.
-	Radius of Gyration (RoG). It is a measure in meters of the area covered by a person over a day. A centroid is calculated for all the places (pauses) visited during a day and a weighted distance between all places and the centroid is computed. The weights are proportional to the time spent in each place.
-	Maximum diameter. Largest distance in meters between any two pauses.
-	Avg flight duration. Mean duration of all flights.
-	Avg flight length. Mean length of all flights
-	Std flight duration. The standard deviation of the duration of all flights.
-	Std flight length. The standard deviation of the length of all flights.
-	Significant locations. The number of significant locations visited during the day. Significant locations are computed using k-means clustering over pauses found in the whole monitoring period. The number of clusters is found iterating from 1 to 200 stopping until the centroids of two significant locations are within 400 meters of one another.
-	Significant location entropy. Entropy measurement based on the proportion of time spent at each significant location visited during a day.

.. _screen:

Screen
------

Available epochs: daily, morning, afternoon, evening, and night

Notes. An unlock episode is considered as the time between an unlock event and a lock event. iOS recorded these episodes reliable (albeit duplicated lock events within milliseconds from each other). However, in Android there are multiple consecutive unlock/lock events so we keep the closest pair. This happens because ACTION_SCREEN_OFF and ON are "sent when the device becomes non-interactive which may have nothing to do with the screen turning off" see this link

-	Count on: count of screen on events (only available for Android)
-	Count unlock: count of screen unlock events
-	Diff count on off: For debug purposes, on and off events should come in pairs, difference should be close to zero then.
-	Diff count unlock lock, For debug purposes, unlock and lock events should come in pairs, difference should be close to zero then.
-	Sum duration unlock: sum duration of unlock episodes 
-	Max duration unlock: maximum duration of unlock episodes
-	Min duration unlock: minimum duration of unlock episodes
-	Average duration unlock: average duration of unlock episodes
-	Std duration unlock: standard deviation of the duration of unlock episodes

.. _fitbit-heart-rate:

Fitbit: heart rate
------------------

Available epochs: daily, morning, afternoon, evening, and night

Notes. eart rate zones contain 4 zones: out_of_range zone, fat_burn zone, cardio zone, and peak zone. Please refer to the [Fitbit documentation](https://help.fitbit.com/articles/en_US/Help_article/1565) for the detailed informations of how to define those zones.

- Max hr: maximum heart rate
- Min hr: minimum heart rate
- Avg hr: average heart rate
- Median hr: median heart rate
- Mode hr: mode heart rate
- Std hr: standard deviation of heart rate
- Diff max mode hr: maximum heart rate minus mode heart rate
- Diff min mode hr: mode heart rate minus minimum heart rate
- Entropy hr: entropy of heart rate
- Length out of range: duration of heart rate in out_of_range zone in minute
- Length fat burn: duration of heart rate in fat_burn zone in minute
- Length cardio: duration of heart rate in cardio zone in minute
- Length peak: duration of heart rate in peak zone in minute

.. _fitbit-steps:

Fitbit: steps
-------------

Available epochs: daily, morning, afternoon, evening, and night

Notes. If the step count per minute smaller than the THRESHOLD_ACTIVE_BOUT (default value is 10), it is defined as sedentary status. Otherwise, it is defined as active status. One active/sedentary bout is a period during with the user is under active/sedentary status.

- Sum all steps: total step count
- Max all steps: maximum step count
- Min all steps: minimum step count
- Avg all steps: average step count
- Std all steps: standard deviation of step count
- Count sedentary bout: count of sedentary bouts
- Max duration sedentary bout: maximum duration of sedentary bouts
- Min duration sedentary bout: minimum duration of sedentary bouts
- Avg duration sedentary bout: average duration of sedentary bouts
- Std duration sedentary bout: standard deviation of the duration of sedentary bouts
- Count active bout: count of active bouts
- Max duration active bout: maximum duration of active bouts
- Min duration active bout: minimum duration of active bouts
- Avg duration active bout: average duration of active bouts
- Std duration active bout: standard deviation of the duration of active bouts



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
.. _call_metrics: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L13
.. _call_metrics.R: https://github.com/carissalow/rapids/blob/master/src/features/call_metrics.R
.. _`Bluetooth Config Code`: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/config.yaml#L76
.. _bluetooth_metric: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/rules/features.snakefile#L63
.. _bluetooth_metrics.R: https://github.com/carissalow/rapids/blob/765bb462636d5029a05f54d4c558487e3786b90b/src/features/bluetooth_metrics.R
