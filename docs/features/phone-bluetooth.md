# Phone Bluetooth

Sensor parameters description for `[PHONE_BLUETOOTH]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[CONTAINER]`| Data stream [container](../../datastreams/data-streams-introduction/) (database table, CSV file, etc.) where the bluetooth data is stored

## RAPIDS provider

!!! warning
    The features of this provider are deprecated in favor of `DORYAB` provider (see below).

!!! info "Available time segments and platforms"
    - Available for all time segments
    - Available for Android only

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/phone_bluetooth_raw.csv
    - data/raw/{pid}/phone_bluetooth_with_datetime.csv
    - data/interim/{pid}/phone_bluetooth_features/phone_bluetooth_{language}_{provider_key}.csv
    - data/processed/features/{pid}/phone_bluetooth.csv"
    ```


Parameters description for `[PHONE_BLUETOOTH][PROVIDERS][RAPIDS]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`| Set to `True` to extract `PHONE_BLUETOOTH` features from the `RAPIDS` provider|
|`[FEATURES]` |         Features to be computed, see table below


Features description for `[PHONE_BLUETOOTH][PROVIDERS][RAPIDS]`:

|Feature                    |Units      |Description|
|-------------------------- |---------- |---------------------------|
| {--countscans--}                 | devices | Number of scanned devices during a time segment, a device can be detected multiple times over time and these appearances are counted separately |
| {--uniquedevices--}              | devices | Number of unique devices during a time segment as identified by their hardware (`bt_address`) address                                                          |
| {--countscansmostuniquedevice--} | scans   | Number of scans of the most sensed device within each time segment instance                                              |

!!! note "Assumptions/Observations"
    - From `v0.2.0` `countscans`, `uniquedevices`, `countscansmostuniquedevice` were deprecated because they overlap with the respective features for `ALL` devices of the `PHONE_BLUETOOTH` `DORYAB` provider

## DORYAB provider
This provider is adapted from the work by [Doryab et al](../../citation#doryab-bluetooth). 

!!! info "Available time segments and platforms"
    - Available for all time segments
    - Available for Android only

!!! info "File Sequence"
    ```bash
    - data/raw/{pid}/phone_bluetooth_raw.csv
    - data/raw/{pid}/phone_bluetooth_with_datetime.csv
    - data/interim/{pid}/phone_bluetooth_features/phone_bluetooth_{language}_{provider_key}.csv
    - data/processed/features/{pid}/phone_bluetooth.csv"
    ```


Parameters description for `[PHONE_BLUETOOTH][PROVIDERS][DORYAB]`:

|Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------
|`[COMPUTE]`| Set to `True` to extract `PHONE_BLUETOOTH` features from the `DORYAB` provider|
|`[FEATURES]` |         Features to be computed, see table below. These features are computed for three device categories: `all` devices, `own` devices and `other` devices.


Features description for `[PHONE_BLUETOOTH][PROVIDERS][DORYAB]`:

|Feature&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                     |Units      |Description|
|-------------------------- |---------- |---------------------------|
| countscans                 | scans | Number of scans (rows) from the devices sensed during a time segment instance. The more scans a bluetooth device has the longer it remained within range of the participant's phone |
| uniquedevices              | devices | Number of unique bluetooth devices sensed during a time segment instance as identified by their hardware addresses (`bt_address`) |
| meanscans | scans| Mean of the scans of every sensed device within each time segment instance|
| stdscans | scans| Standard deviation of the scans of every sensed device within each time segment instance|
| countscans{==most==}frequentdevice{==within==}segments | scans   | Number of scans of the **most** sensed device **within** each time segment instance|
| countscans{==least==}frequentdevice{==within==}segments| scans| Number of scans of the **least** sensed device **within** each time segment instance |
| countscans{==most==}frequentdevice{==across==}segments | scans   | Number of scans of the **most** sensed device **across** time segment instances of the same type|
| countscans{==least==}frequentdevice{==across==}segments| scans| Number of scans of the **least** sensed device **across** time segment instances of the same type per device|
| countscans{==most==}frequentdevice{==acrossdataset==} | scans   | Number of scans of the **most** sensed device **across** the entire dataset of every participant|
| countscans{==least==}frequentdevice{==acrossdataset==}| scans| Number of scans of the **least** sensed device **across** the entire dataset of every participant |


!!! note "Assumptions/Observations"
    - Devices are classified as belonging to the participant (`own`) or to other people (`others`) using k-means based on the number of times and the number of days each device was detected across each participant's dataset. See [Doryab et al](../../citation#doryab-bluetooth) for more details.
    - If ownership cannot be computed because all devices were detected on only one day, they are all considered as `other`. Thus `all` and `other` features will be equal. The likelihood of this scenario decreases the more days of data you have.
    - When searching for the most frequent device across 30-minute segments, the search range is equivalent to the sum of all segments of the same time period. For instance, the `countscansmostfrequentdeviceacrosssegments` for the time segment (`Fri 00:00:00, Fri 00:29:59`) will get the count in that segment of the most frequent device found within all (`00:00:00, 00:29:59`) time segments. To find `countscansmostfrequentdeviceacrosssegments` for `other` devices, the search range needs to filter out all `own` devices. But no need to do so for `countscansmostfrequentdeviceacrosssedataset`. The most frequent device across the dataset stays the same for `countscansmostfrequentdeviceacrossdatasetall`, `countscansmostfrequentdeviceacrossdatasetown` and `countscansmostfrequentdeviceacrossdatasetother`. Same rule applies to the least frequent device across the dataset. 
    - The most and least frequent devices will be the same across time segment instances and across the entire dataset when every time segment instance covers every hour of a dataset. For example, daily segments (00:00 to 23:59) fall in this category but morning segments (06:00am to 11:59am) or periodic 30-minute segments don't.

    ??? info "Example"
        
        ??? example "Simplified raw bluetooth data"
            The following is a simplified example with bluetooth data from three days and two time segments: morning and afternoon. There are two `own` devices: `5C836F5-487E-405F-8E28-21DBD40FA4FF` detected seven times across two days and `499A1EAF-DDF1-4657-986C-EA5032104448` detected eight times on a single day.
            ```csv
            local_date	segment	    bt_address                              own_device
            2016-11-29	morning	    55C836F5-487E-405F-8E28-21DBD40FA4FF              1
            2016-11-29	morning	    55C836F5-487E-405F-8E28-21DBD40FA4FF              1
            2016-11-29	morning	    55C836F5-487E-405F-8E28-21DBD40FA4FF              1
            2016-11-29	morning	    55C836F5-487E-405F-8E28-21DBD40FA4FF              1
            2016-11-29	morning	    48872A52-68DE-420D-98DA-73339A1C4685              0
            2016-11-29	afternoon	55C836F5-487E-405F-8E28-21DBD40FA4FF              1
            2016-11-29	afternoon	48872A52-68DE-420D-98DA-73339A1C4685              0
            2016-11-30	morning	    55C836F5-487E-405F-8E28-21DBD40FA4FF              1
            2016-11-30	morning	    48872A52-68DE-420D-98DA-73339A1C4685              0
            2016-11-30	morning	    25262DC7-780C-4AD5-AD3A-D9776AEF7FC1              0
            2016-11-30	morning	    5B1E6981-2E50-4D9A-99D8-67AED430C5A8              0
            2016-11-30	morning	    5B1E6981-2E50-4D9A-99D8-67AED430C5A8              0
            2016-11-30	afternoon	55C836F5-487E-405F-8E28-21DBD40FA4FF              1
            2017-05-07	morning	    5C5A9C41-2F68-4CEB-96D0-77DE3729B729              0
            2017-05-07	morning	    25262DC7-780C-4AD5-AD3A-D9776AEF7FC1              0
            2017-05-07	morning	    5B1E6981-2E50-4D9A-99D8-67AED430C5A8              0
            2017-05-07	morning	    6C444841-FE64-4375-BC3F-FA410CDC0AC7              0
            2017-05-07	morning	    4DC7A22D-9F1F-4DEF-8576-086910AABCB5              0
            2017-05-07	afternoon	5B1E6981-2E50-4D9A-99D8-67AED430C5A8              0
            2017-05-07  afternoon   499A1EAF-DDF1-4657-986C-EA5032104448              1
            2017-05-07  afternoon   499A1EAF-DDF1-4657-986C-EA5032104448              1
            2017-05-07  afternoon   499A1EAF-DDF1-4657-986C-EA5032104448              1
            2017-05-07  afternoon   499A1EAF-DDF1-4657-986C-EA5032104448              1
            2017-05-07  afternoon   499A1EAF-DDF1-4657-986C-EA5032104448              1
            2017-05-07  afternoon   499A1EAF-DDF1-4657-986C-EA5032104448              1
            2017-05-07  afternoon   499A1EAF-DDF1-4657-986C-EA5032104448              1
            2017-05-07  afternoon   499A1EAF-DDF1-4657-986C-EA5032104448              1
            ```
        

        

        ??? example "The most and least frequent `OTHER` devices (`own_device == 0`) during morning segments"
            The most and least frequent `ALL`|`OWN`|`OTHER` devices are computed within each time segment instance, across time segment instances of the same type and across the entire dataset of each person. These are the most and least frequent devices for `OTHER` devices during morning segments.
            ```csv
            most frequent device across 2016-11-29 morning:   '48872A52-68DE-420D-98DA-73339A1C4685'  (this device is the only one in this instance)
            least frequent device across 2016-11-29 morning:  '48872A52-68DE-420D-98DA-73339A1C4685'  (this device is the only one in this instance)
            most frequent device across 2016-11-30 morning:   '5B1E6981-2E50-4D9A-99D8-67AED430C5A8'
            least frequent device across 2016-11-30 morning:  '25262DC7-780C-4AD5-AD3A-D9776AEF7FC1'  (when tied, the first occurance is chosen)
            most frequent device across 2017-05-07 morning:   '25262DC7-780C-4AD5-AD3A-D9776AEF7FC1'  (when tied, the first occurance is chosen)
            least frequent device across 2017-05-07 morning:  '25262DC7-780C-4AD5-AD3A-D9776AEF7FC1'  (when tied, the first occurance is chosen)
            
            most frequent across morning segments:            '5B1E6981-2E50-4D9A-99D8-67AED430C5A8'
            least frequent across morning segments:           '6C444841-FE64-4375-BC3F-FA410CDC0AC7' (when tied, the first occurance is chosen)
            
            most frequent across dataset:                     '499A1EAF-DDF1-4657-986C-EA5032104448' (only taking into account "morning" segments)
            least frequent across dataset:                    '4DC7A22D-9F1F-4DEF-8576-086910AABCB5' (when tied, the first occurance is chosen)
            ```

        ??? example "Bluetooth features for  `OTHER` devices and morning segments"
            For brevity we only show the following features for morning segments:
            ```yaml
            OTHER: 
                DEVICES: ["countscans", "uniquedevices", "meanscans", "stdscans"]
                SCANS_MOST_FREQUENT_DEVICE: ["withinsegments", "acrosssegments", "acrossdataset"]
            ```

            Note that `countscansmostfrequentdeviceacrossdatasetothers` is all `0`s because `499A1EAF-DDF1-4657-986C-EA5032104448` is excluded from the count as is labelled as an `own` device (not `other`).
            ```csv
            local_segment       countscansothers	uniquedevicesothers	meanscansothers	stdscansothers	countscansmostfrequentdevicewithinsegmentsothers	countscansmostfrequentdeviceacrosssegmentsothers	countscansmostfrequentdeviceacrossdatasetothers
            2016-11-29-morning	1	                1	                1.000000	    NaN             1	                                                0.0	                                                0.0
            2016-11-30-morning	4	                3	                1.333333	    0.57735	        2	                                                2.0	                                                2.0
            2017-05-07-morning	5	                5	                1.000000	    0.00000	        1	                                                1.0	                                                1.0
            ```
