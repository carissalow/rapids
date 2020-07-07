.. _test-cases:

Test Cases
-----------

Along with the continued development and the addition of new sensors and features to the RAPIDS pipeline, tests for the currently available sensors and features are being implemented. Since this is a Work In Progress this page will be updated with the list of sensors and features for which testing is available. For each of the sensors listed a description of the data used for testing (test cases) are outline. Currently for all intent and testing purposes the ``tests/data/raw/test01/`` contains all the test data files for testing android data formats and ``tests/data/raw/test02/`` contains all the test data files for testing iOS data formats. It follows that the expected (verified output) are contained in the ``tests/data/processed/test01/`` and ``tests/data/processed/test02/`` for Android and iOS respectively. ``tests/data/raw/test03/`` and ``tests/data/raw/test04/`` contain data files for testing empty raw data files for android and iOS respectively. 

List of Sensor with Tests
^^^^^^^^^^^^^^^^^^^^^^^^^^
The following is a list of the sesors that testing is currently available. 


Messages (SMS)
"""""""""""""""

    - The raw message data file contains data for 2 separate days. 
    - The data for the first day contains records 5 records for every ``epoch``.
    - The second day's data contains 6 records for each of only 2 ``epoch`` (currently ``morning`` and ``evening``)
    - The raw message data contains records for both ``message_types`` (i.e. ``recieved`` and ``sent``) in both days in all epochs. The number records with each ``message_types`` per epoch is randomly distributed There is at least one records with each ``message_types`` per epoch.
    - There is one raw message data file each, as described above, for testing both iOS and Android data. 
    - There is also an additional empty data file for both android and iOS for testing empty data files

Calls
"""""""

    Due to the difference in the format of the raw call data for iOS and Android (see the **Assumptions/Observations** section of :ref:`Calls<call-sensor-doc>`) the following is the expected results the ``calls_with_datetime_unified.csv``. This would give a better idea of the use cases being tested since the ``calls_with_datetime_unified.csv`` would make both the iOS and Android data comparable. 

    - The call data would contain data for 2 days. 
    - The data for the first day contains 6 records for every ``epoch``. 
    - The second day's data contains 6 records for each of only 2 ``epoch`` (currently ``morning`` and ``evening``)
    - The call data contains records for all ``call_types`` (i.e. ``incoming``, ``outgoing`` and ``missed``) in both days in all epochs. The number records with each of the ``call_types`` per epoch is randomly distributed. There is at least one records with each ``call_types`` per epoch.
    - There is one call data file each, as described above, for testing both iOS and Android data. 
    - There is also an additional empty data file for both android and iOS for testing empty data files

Screen
""""""""

    Due to the difference in the format of the raw screen data for iOS and Android (see the **Assumptions/Observations** section of :ref:`Screen<screen-sensor-doc>`) the following is the expected results the ``screen_deltas.csv``. This would give a better idea of the use cases being tested since the ``screen_deltas.csv`` would make both the iOS and Android data comparable. These files are used to calculate the features for the screen sensor. 

    - The screen delta data file contains data for 1 day. 
    - The screen delta data contains 1 record to represent an ``unlock`` episode that falls within an ``epoch`` for every ``epoch``. 
    - The screen delta data contains 1 record to represent an ``unlock`` episode that falls across the boundary of 2 epochs. Namely the ``unlock`` episode starts in one epoch and ends in the next, thus there is a record for ``unlock`` episodes that fall across ``night`` to ``morning``, ``morning`` to ``afternoon`` and finally ``afternoon`` to ``night``
    - The testing is done for ``unlock`` episode_type.
    - There is one screen data file each for testing both iOS and Android data formats.
    - There is also an additional empty data file for both android and iOS for testing empty data files

Battery
"""""""""

    Due to the difference in the format of the raw battery data for iOS and Android as well as versions of iOS (see the **Assumptions/Observations** section of :ref:`Battery<battery-sensor-doc>`) the following is the expected results the ``battery_deltas.csv``. This would give a better idea of the use cases being tested since the ``battery_deltas.csv`` would make both the iOS and Android data comparable. These files are used to calculate the features for the battery sensor. 

    - The battery delta data file contains data for 1 day. 
    - The battery delta data contains 1 record each for a ``charging`` and ``discharging`` episode that falls within an ``epoch`` for every ``epoch``. Thus for the ``daily`` epoch there would be multiple ``charging`` and ``discharging`` episodes
    - Since either a ``charging`` episode or a ``discharging`` episode and not both can occur across epochs, in order to test epsiodes that occur across epochs alternating episodes of ``charging`` and ``discharging`` episodes that fall across ``night`` to ``morning``, ``morning`` to ``afternoon`` and finally ``afternoon`` to ``night`` are present in the battery delta data. This starts with a ``discharging`` episode that begins in ``night`` and end in ``morning``.
    - There is one battery data file each, for testing both iOS and Android data formats.
    - There is also an additional empty data file for both android and iOS for testing empty data files

Bluetooth
""""""""""

    - The raw bluetooth data file contains data for 1 day. 
    - The raw bluetooth data contains at least 2 records for each ``epoch``. Each ``epoch`` has a record with a ``timestamp`` for the beginning boundary for that ``epoch`` and a record with a ``timestamp`` for the ending boundary for that ``epoch``. (e.g. For the ``morning`` epoch there is a record with a ``timestamp`` for ``6:00AM`` and another record with a ``timestamp`` for ``11:59:59AM``. These are to test edge cases) 
    - An option of 5 bluetooth devices are randomly distributed throughout the data records.
    - There is one raw bluetooth data file each, for testing both iOS and Android data formats.
    - There is also an additional empty data file for both android and iOS for testing empty data files.

WIFI
"""""

    - The raw WIFI data file contains data for 1 day. 
    - The raw WIFI data contains at least 2 records for each ``epoch``. Each ``epoch`` has a record with a ``timestamp`` for the beginning boundary for that ``epoch`` and a record with a ``timestamp`` for the ending boundary for that ``epoch``. (e.g. For the ``morning`` epoch there is a record with a ``timestamp`` for ``6:00AM`` and another record with a ``timestamp`` for ``11:59:59AM``. These are to test edge cases) 
    - An option of 5 access point devices is randomly distributed throughout the data records.
    - There is one raw WIFI data file each, for testing both iOS and Android data formats.
    - There is also an additional empty data file for both android and iOS for testing empty data files.
