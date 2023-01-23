
# Configuration

You need to follow these steps to configure your RAPIDS deployment before you can extract behavioral features.

0. Verify RAPIDS can process your [data streams](#supported-data-streams)
3. Create your [participants files](#participant-files)
4. Select what [time segments](#time-segments) you want to extract features on
2. Select the [timezone of your study](#timezone-of-your-study)
5. Configure your [data streams](#data-stream-configuration)
6. Select what [sensors and features](#sensor-and-features-to-process) you want to process

When you are done with this configuration, go to [executing RAPIDS](../execution).

!!! hint
    Every time you see `config["KEY"]` or `[KEY]` in these docs, we are referring to the corresponding key in the `config.yaml` file.

---

## Supported data streams

A data stream refers to sensor data collected using a specific **device** with a specific **format** and stored in a specific **container**. For example, the `aware_mysql` data stream handles smartphone data (**device**) collected with the [AWARE Framework](https://awareframework.com/) (**format**) stored in a MySQL database (**container**).

Check the table in [introduction to data streams](../../datastreams/data-streams-introduction) to know what data streams we support. If your data stream is supported, continue to the next configuration section, **you will use its label later in this guide** (e.g. `aware_mysql`). If your steam is not supported, but you want to implement it, follow the tutorial to [add support for new data streams](../../datastreams/add-new-data-streams) and [open a new discussion](https://github.com/carissalow/rapids/discussions) in Github with any questions.

---

## Participant files

Participant files link together multiple devices (smartphones and wearables) to specific participants and identify them throughout RAPIDS. You can create these files manually or [automatically](#automatic-creation-of-participant-files). Participant files are stored in `data/external/participant_files/pxx.yaml` and follow a unified [structure](#structure-of-participants-files).

??? important "Remember to modify the `config.yaml` file with your PIDS"
    The list `PIDS` in `config.yaml` needs to have the participant file names of the people you want to process. For example, if you created `p01.yaml`, `p02.yaml` and `p03.yaml` files in `/data/external/participant_files/ `, then `PIDS` should be:
    ```yaml
    PIDS: [p01, p02, p03] 
    ```

??? info "Optional: Migrating participants files with the old format"
    If you were using the pre-release version of RAPIDS with participant files in plain text (as opposed to yaml), you could run the following command, and your old files will be converted into yaml files stored in `data/external/participant_files/`

    ```bash
    python tools/update_format_participant_files.py
    ```

### Structure of participants files

??? example "Example of the structure of a participant file"

    In this example, the participant used an android phone, an ios phone, a Fitbit device, and an Empatica device throughout the study between April 23rd, 2020, and October 28th, 2020

     If your participants didn't use a `[PHONE]`, `[FITBIT]` or `[EMPATICA]` device, it is not necessary to include that section in their participant file. In other words, you can analyze data from 1 or more devices per participant.

    ```yaml
    PHONE:
      DEVICE_IDS: [a748ee1a-1d0b-4ae9-9074-279a2b6ba524, dsadas-2324-fgsf-sdwr-gdfgs4rfsdf43]
      PLATFORMS: [android,ios]
      LABEL: test01
      START_DATE: 2020-04-23
      END_DATE: 2020-10-28
    FITBIT:
      DEVICE_IDS: [fitbit1]
      LABEL: test01
      START_DATE: 2020-04-23
      END_DATE: 2020-10-28
    EMPATICA:
      DEVICE_IDS: [empatica1]
      LABEL: test01
      START_DATE: 2020-04-23
      END_DATE: 2020-10-28
    ```

=== "[PHONE]"

    | Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            | Description                                                                                                                                                                                                                                                                                                                                |
    |-------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
    | `[DEVICE_IDS]` | An array of the strings that uniquely identify each smartphone, you can have more than one for when participants changed phones in the middle of the study.                                                                           |
    | `[PLATFORMS]`  | An array that specifies the OS of each smartphone in  `[DEVICE_IDS]` , use a combination of  `android`  or  `ios`  (we support participants that changed platforms in the middle of your study!). You can set `[PLATFORMS]: [infer]`, and RAPIDS will infer them automatically (each phone data stream infer this differently, e.g., `aware_mysql` uses the `aware_device` table). |
    | `[LABEL]`      | A string that is used in reports and visualizations.        |
    | `[START_DATE]` | A string with format `YYYY-MM-DD` or `YYYY-MM-DD HH:MM:SS`. Only data collected  *after*  this date-time will be included in the analysis. By default, `YYYY-MM-DD` is interpreted as `YYYY-MM-DD 00:00:00`.                               |
    | `[END_DATE]`   | A string with format `YYYY-MM-DD` or `YYYY-MM-DD HH:MM:SS`. Only data collected  *before*  this date-time will be included in the analysis. By default, `YYYY-MM-DD` is interpreted as `YYYY-MM-DD 00:00:00`.                              |

=== "[FITBIT]"

    | Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;              | Description                                                                                               |
    |------------------|-----------------------------------------------------------------------------------------------------------|
    | `[DEVICE_IDS]`   | An array of the strings that uniquely identify each Fitbit, you can have more than one in case the participant changed devices in the middle of the study. |
    | `[LABEL]`        | A string that is used in reports and visualizations.                                              |
    | `[START_DATE]`   | A string with format `YYYY-MM-DD` or `YYYY-MM-DD HH:MM:SS`. Only data collected  *after*  this date-time will be included in the analysis. By default, `YYYY-MM-DD` is interpreted as `YYYY-MM-DD 00:00:00`.                               |
    | `[END_DATE]`     | A string with format `YYYY-MM-DD` or `YYYY-MM-DD HH:MM:SS`. Only data collected  *before*  this date-time will be included in the analysis. By default, `YYYY-MM-DD` is interpreted as `YYYY-MM-DD 00:00:00`.                              |

=== "[EMPATICA]"

    | Key&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;              | Description                                                                                               |
    |------------------|-----------------------------------------------------------------------------------------------------------|
    | `[DEVICE_IDS]`   | An array of the strings that uniquely identify each Empatica device used by this participant. Since the most common use case involves having multiple zip files from a single device for each person, set this device id to an arbitrary string (we usually use their `pid`) |
    | `[LABEL]`        | A string that is used in reports and visualizations.                                                                                              |
    | `[START_DATE]`   | A string with format `YYYY-MM-DD` or `YYYY-MM-DD HH:MM:SS`. Only data collected  *after*  this date-time will be included in the analysis. By default, `YYYY-MM-DD` is interpreted as `YYYY-MM-DD 00:00:00`.                               |
    | `[END_DATE]`     | A string with format `YYYY-MM-DD` or `YYYY-MM-DD HH:MM:SS`. Only data collected  *before*  this date-time will be included in the analysis. By default, `YYYY-MM-DD` is interpreted as `YYYY-MM-DD 00:00:00`.                              |
### Automatic creation of participant files

You can use a CSV file with a row per participant to automatically create participant files. 

??? "`AWARE_DEVICE_TABLE` was deprecated"
    In previous versions of RAPIDS, you could create participant files automatically using the `aware_device` table. We deprecated this option, but you can still achieve the same results if you export the output of the following SQL query as a CSV file and follow the instructions below:
    
    ```sql
    SELECT device_id, device_id as fitbit_id, CONCAT("p", _id) as empatica_id, CONCAT("p", _id) as pid, if(brand = "iPhone", "ios", "android") as platform, CONCAT("p", _id)  as label, DATE_FORMAT(FROM_UNIXTIME((timestamp/1000)- 86400), "%Y-%m-%d") as start_date, CURRENT_DATE as end_date from aware_device order by _id;
    ```

In your `config.yaml`:

1. Set `CSV_FILE_PATH` to a CSV file path that complies with the specs described below
2. Set the devices (`PHONE`, `FITBIT`, `EMPATICA`) `[ADD]` flag to `TRUE` depending on what devices you used in your study.

```yaml
CREATE_PARTICIPANT_FILES:
  CSV_FILE_PATH: "your_path/to_your.csv"
  PHONE_SECTION:
    ADD: TRUE # or FALSE
    IGNORED_DEVICE_IDS: []
  FITBIT_SECTION:
    ADD: TRUE # or FALSE
    IGNORED_DEVICE_IDS: []
  EMPATICA_SECTION:
    ADD: TRUE # or FALSE
    IGNORED_DEVICE_IDS: []
```

Your CSV file (`[CSV_FILE_PATH]`) should have the following columns (headers), but the values within each column can be empty:

| Column           | Description                                                                                               |
|------------------|-----------------------------------------------------------------------------------------------------------|
| device_id        | Phone device id. Separate multiple ids with `;`   |
| fitbit_id        | Fitbit device id. Separate multiple ids with `;`  |
| empatica_id      | Empatica device id. Since the most common use case involves having various zip files from a single device for each person, set this device id to an arbitrary string (we usually use their `pid`)  |
| pid              | Unique identifiers with the format pXXX (your participant files will be named with this string)            |
| platform         | Use `android`, `ios` or `infer` as explained above, separate values with `;`            |
| label            | A human-readable string that is used in reports and visualizations.                                       |
| start_date       | A string with format `YYY-MM-DD` or `YYYY-MM-DD HH:MM:SS`. By default, `YYYY-MM-DD` is interpreted as `YYYY-MM-DD 00:00:00`. |
| end_date         | A string with format `YYY-MM-DD` or `YYYY-MM-DD HH:MM:SS`. By default, `YYYY-MM-DD` is interpreted as `YYYY-MM-DD 00:00:00`. |

!!! example
    We added white spaces to this example to make it easy to read, but you don't have to.

    ```csv
    device_id                                                                ,fitbit_id, empatica_id ,pid ,label ,platform    ,start_date ,end_date
    a748ee1a-1d0b-4ae9-9074-279a2b6ba524;dsadas-2324-fgsf-sdwr-gdfgs4rfsdf43 ,fitbit1  , p01         ,p01 ,julio ,android;ios ,2020-01-01 ,2021-01-01
    4c4cf7a1-0340-44bc-be0f-d5053bf7390c                                     ,fitbit2  , p02         ,p02 ,meng  ,ios         ,2021-01-01 ,2022-01-01
    ```

Then run 

```bash
snakemake -j1 create_participants_files
```

---

## Time Segments

Time segments (or epochs) are the time windows on which you want to extract behavioral features. For example, you might want to process data every day, every morning, or only during weekends. RAPIDS offers three categories of time segments that are flexible enough to cover most use cases: **frequency** (short time windows every day), **periodic** (arbitrary time windows on any day), and **event** (arbitrary time windows around events of interest). See also our [examples](#segment-examples).

=== "Frequency Segments"

    These segments are computed every day, and all have the same duration (for example, 30 minutes). Set the following keys in your `config.yaml`

    ```yaml
    TIME_SEGMENTS: &time_segments
      TYPE: FREQUENCY
      FILE: "data/external/your_frequency_segments.csv"
      INCLUDE_PAST_PERIODIC_SEGMENTS: FALSE
    ```

    The file pointed by `[TIME_SEGMENTS][FILE]` should have the following format and only have 1 row.

    | Column | Description                                                          |
    |--------|----------------------------------------------------------------------|
    | label  | A string that is used as a prefix in the name of your time segments   |
    | length | An integer representing the duration of your time segments in minutes |

    !!! example

        ```csv
        label,length
        thirtyminutes,30
        ```
        
        This configuration will compute 48 time segments for every day when any data from any participant was sensed. For example:

        ```csv
        start_time,length,label
        00:00,30,thirtyminutes0000
        00:30,30,thirtyminutes0001
        01:00,30,thirtyminutes0002
        01:30,30,thirtyminutes0003
        ...
        ```

=== "Periodic Segments"

    These segments can be computed every day or on specific days of the week, month, quarter, and year. Their minimum duration is 1 minute, but they can be as long as you want. Set the following keys in your `config.yaml`.

    ```yaml
    TIME_SEGMENTS: &time_segments
      TYPE: PERIODIC
      FILE: "data/external/your_periodic_segments.csv"
      INCLUDE_PAST_PERIODIC_SEGMENTS: FALSE # or TRUE
    ```

    If `[INCLUDE_PAST_PERIODIC_SEGMENTS]` is set to `TRUE`, RAPIDS will consider instances of your segments back enough in the past to include the first row of data of each participant. For example, if the first row of data from a participant happened on Saturday, March 7th, 2020, and the requested segment duration is 7 days starting on every Sunday, the first segment to be considered would begin on Sunday, March 1st if `[INCLUDE_PAST_PERIODIC_SEGMENTS]` is `TRUE` or on Sunday, March 8th if `FALSE`.

    The file pointed by `[TIME_SEGMENTS][FILE]` should have the following format and can have multiple rows.

    | Column        | Description                                                                                                                                                                                                   |
    |---------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
    | label         | A string that is used as a prefix in the name of your time segments. It has to be **unique** between rows                                                                                                      |
    | start_time    | A string with format `HH:MM:SS` representing the starting time of this segment on any day                                                                                                                                 |
    | length        | A string representing the length of this segment. It can have one or more of the following strings **`XXD XXH XXM XXS`** to represent days, hours, minutes, and seconds. For example, `7D 23H 59M 59S`                        |
    | repeats_on    | One of the following options `every_day`, `wday`, `qday`, `mday`, and `yday`. The last four represent a week, quarter, month, and year day                                                                        |
    | repeats_value | An integer complementing `repeats_on`. If you set `repeats_on` to `every_day`, set this to `0`, otherwise `1-7` represent a `wday` starting from Mondays, `1-31` represent a `mday`, `1-91` represent a `qday`, and `1-366` represent a `yday` |

    !!! example

        ```csv
        label,start_time,length,repeats_on,repeats_value
        daily,00:00:00,23H 59M 59S,every_day,0
        morning,06:00:00,5H 59M 59S,every_day,0
        afternoon,12:00:00,5H 59M 59S,every_day,0
        evening,18:00:00,5H 59M 59S,every_day,0
        night,00:00:00,5H 59M 59S,every_day,0
        ```

        This configuration will create five segment instances (`daily`, `morning`, `afternoon`, `evening`, `night`) on any given day (`every_day` set to 0). The `daily` segment will start at midnight and last `23:59:59`; the other four segments will begin at 6am, 12pm, 6pm, and 12am, respectively, and last for `05:59:59`. 

=== "Event segments"

    These segments can be computed before or after an event of interest (defined as any UNIX timestamp). Their minimum duration is 1 minute, but they can be as long as you want. The start of each segment can be shifted backward or forwards from the specified timestamp. Set the following keys in your `config.yaml`.

    ```yaml
    TIME_SEGMENTS: &time_segments
      TYPE: EVENT
      FILE: "data/external/your_event_segments.csv"
      INCLUDE_PAST_PERIODIC_SEGMENTS: FALSE # or TRUE
    ```

    The file pointed by `[TIME_SEGMENTS][FILE]` should have the following format and can have multiple rows.

    | Column        | Description                                                                                                                                                                                                   |
    |---------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
    | label         | A string that is used as a prefix in the name of your time segments. If labels are unique, every segment is independent; if two or more segments have the same label, their data will be grouped when computing auxiliary data for features like the `most frequent contact` for calls (the most frequent contact will be calculated across all these segments). There cannot be two *overlapping* event segments with the same label (RAPIDS will throw an error)                                                                                                      |
    | event_timestamp    | A UNIX timestamp that represents the moment an event of interest happened (clinical relapse, survey, readmission, etc.). The corresponding time segment will be computed around this moment using `length`, `shift`, and `shift_direction`                                                                                            |
    | length        | A string representing the length of this segment. It can have one or more of the following keys `XXD XXH XXM XXS` to represent days, hours, minutes, and seconds. For example, `7D 23H 59M 59S`                        |
    | shift    | A string representing the time shift from `event_timestamp`. It can have one or more of the following keys `XXD XXH XXM XXS` to represent days, hours, minutes, and seconds. For example, `7D 23H 59M 59S`. Use this value to change the start of a segment with respect to its `event_timestamp`. For example, set this variable to `1H` to create a segment that starts 1 hour from an event of interest (`shift_direction` determines if it's before or after).                                                                        |
    | shift_direction | An integer representing whether the `shift` is before (`-1`) or after (`1`) an `event_timestamp` |
    |device_id| The device id (smartphone or Fitbit) to whom this segment belongs to. You have to create a line in this event segment file for each event of a participant that you want to analyze. If you have participants with multiple device ids, you can choose any of them|

    !!! example
        ```csv
        label,event_timestamp,length,shift,shift_direction,device_id
        stress1,1587661220000,1H,5M,1,a748ee1a-1d0b-4ae9-9074-279a2b6ba524
        stress2,1587747620000,4H,4H,-1,a748ee1a-1d0b-4ae9-9074-279a2b6ba524
        stress3,1587906020000,3H,5M,1,a748ee1a-1d0b-4ae9-9074-279a2b6ba524
        stress4,1584291600000,7H,4H,-1,a748ee1a-1d0b-4ae9-9074-279a2b6ba524
        stress5,1588172420000,9H,5M,-1,a748ee1a-1d0b-4ae9-9074-279a2b6ba524
        mood,1587661220000,1H,0,0,a748ee1a-1d0b-4ae9-9074-279a2b6ba524
        mood,1587747620000,1D,0,0,a748ee1a-1d0b-4ae9-9074-279a2b6ba524
        mood,1587906020000,7D,0,0,a748ee1a-1d0b-4ae9-9074-279a2b6ba524
        ```
        
        This example will create eight segments for a single participant (`a748ee1a...`), five independent `stressX` segments with various lengths (1,4,3,7, and 9 hours). Segments `stress1`, `stress3`, and `stress5` are shifted forwards by 5 minutes, and `stress2` and `stress4` are shifted backward by 4 hours (that is, if the `stress4` event happened on March 15th at 1pm EST (`1584291600000`), the time segment will start on that day at 9am and end at 4pm). 
        
        The three `mood` segments are 1 hour, 1 day, and 7 days long and have no shift. In addition, these `mood` segments are grouped together, meaning that although RAPIDS will compute features on each one of them, some information for such computation will be extracted from all three segments, for example, the phone contact that called a participant the most, or the location clusters visited by a participant.

    ??? info "Date time labels of event segments"
        In the final feature file, you will find a row per event segment. The `local_segment` column of each row has a `label`, a start date-time string, and an end date-time string.

        ```bash
        weeklysurvey2060#2020-09-12 01:00:00,2020-09-18 23:59:59
        ```

        All sensor data is always segmented based on timestamps, and the date-time strings are attached for informative purposes. For example, you can plot your features based on these strings. 

        When you configure RAPIDS to work with a single time zone, such time zone code will be used to convert start/end timestamps (the ones you typed in the event segments file) into start/end date-time strings. However, when you configure RAPIDS to work with multiple time zones, RAPIDS will use the most common time zone across all devices of every participant to do the conversion. The most common time zone is the one in which a participant spent the most time.

        In practical terms, this means that the date-time strings of event segments that happened in uncommon time zones will have shifted start/end date-time labels. However, the data within each segment was correctly filtered based on timestamps.

### Segment Examples

=== "5-minutes"
    Use the following `Frequency` segment file to create 288 (12 * 60 * 24) 5-minute segments starting from midnight of every day in your study
    ```csv
    label,length
    fiveminutes,5
    ```
=== "Daily"
    Use the following `Periodic` segment file to create daily segments starting from midnight of every day in your study
    ```csv
    label,start_time,length,repeats_on,repeats_value
    daily,00:00:00,23H 59M 59S,every_day,0
    ```
=== "Morning"
    Use the following `Periodic` segment file to create morning segments starting at 06:00:00 and ending at 11:59:59 of every day in your study
    ```csv
    label,start_time,length,repeats_on,repeats_value
    morning,06:00:00,5H 59M 59S,every_day,0
    ```
=== "Overnight"
    Use the following `Periodic` segment file to create overnight segments starting at 20:00:00 and ending at 07:59:59 (next day) of every day in your study
    ```csv
    label,start_time,length,repeats_on,repeats_value
    morning,20:00:00,11H 59M 59S,every_day,0
    ```
=== "Weekly"
    Use the following `Periodic` segment file to create **non-overlapping** weekly segments starting at midnight of every **Monday** in your study
    ```csv
    label,start_time,length,repeats_on,repeats_value
    weekly,00:00:00,6D 23H 59M 59S,wday,1
    ```
    Use the following `Periodic` segment file to create **overlapping** weekly segments starting at midnight of **every day** in your study
    ```csv
    label,start_time,length,repeats_on,repeats_value
    weekly,00:00:00,6D 23H 59M 59S,every_day,0
    ```
=== "Week-ends"
    Use the following `Periodic` segment file to create week-end segments starting at midnight of every **Saturday** in your study
    ```csv
    label,start_time,length,repeats_on,repeats_value
    weekend,00:00:00,1D 23H 59M 59S,wday,6
    ```
=== "Around surveys"
    Use the following `Event` segment file to create two 2-hour segments that start 1 hour before surveys answered by 3 participants
    ```csv
    label,event_timestamp,length,shift,shift_direction,device_id
    survey1,1587661220000,2H,1H,-1,a748ee1a-1d0b-4ae9-9074-279a2b6ba524
    survey2,1587747620000,2H,1H,-1,a748ee1a-1d0b-4ae9-9074-279a2b6ba524
    survey1,1587906020000,2H,1H,-1,rqtertsd-43ff-34fr-3eeg-efe4fergregr
    survey2,1584291600000,2H,1H,-1,rqtertsd-43ff-34fr-3eeg-efe4fergregr
    survey1,1588172420000,2H,1H,-1,klj34oi2-8frk-2343-21kk-324ljklewlr3
    survey2,1584291600000,2H,1H,-1,klj34oi2-8frk-2343-21kk-324ljklewlr3
    ```
--- 

## Timezone of your study

### Single timezone

If your study only happened in a single time zone or you want to ignore short trips of your participants to different time zones, select the appropriate code from this [list](https://en.wikipedia.org/wiki/List_of_tz_database_time_zones) and change the following config key. Double-check your timezone code pick; for example, US Eastern Time is `America/New_York`, not `EST`.

``` yaml
TIMEZONE: 
    TYPE: SINGLE
    TZCODE: America/New_York
```

### Multiple timezones

If your participants lived in different time zones or they traveled across time zones, and you know when participants' devices were in a specific time zone, RAPIDS can use this data to process your data streams with the correct date-time. You need to provide RAPIDS with the time zone data in a CSV file (`[TZCODES_FILE]`) in the format described below.

``` yaml
TIMEZONE: 
    TYPE: MULTIPLE
    SINGLE:
      TZCODE: America/New_York
    MULTIPLE:
      TZCODES_FILE: path_to/time_zones_csv.file
      IF_MISSING_TZCODE: STOP
      DEFAULT_TZCODE: America/New_York
      FITBIT: 
        ALLOW_MULTIPLE_TZ_PER_DEVICE: False
        INFER_FROM_SMARTPHONE_TZ: False
```

Parameters for `[TIMEZONE]`

|Parameter &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Description |
|--|--|
|`[TYPE]`| Either `SINGLE` or `MULTIPLE` as explained above |
|`[SINGLE][TZCODE]`| The time zone code from this [list](https://en.wikipedia.org/wiki/List_of_tz_database_time_zones) to be used across all devices |
|`[MULTIPLE][TZCODES_FILE]`| A CSV file containing the time zones in which participants' devices sensed data (see the required format below). Multiple devices can be linked to the same person. Read more in [Participants Files](#participant-files) |
|`[MULTIPLE][IF_MISSING_TZCODE]`| When a device is missing from `[TZCODES_FILE]` Set this flag to `STOP` to stop RAPIDS execution and show an error, or to `USE_DEFAULT` to assign the time zone specified in `[DEFAULT_TZCODE]` to any such devices  |
|`[MULTIPLE][FITBIT][ALLOW_MULTIPLE_TZ_PER_DEVICE]`| You only need to care about this flag if one or more Fitbit devices sensed data in one or more time zones, and you want RAPIDS to take into account this in its feature computation. Read more in  "How does RAPIDS handle Fitbit devices?" below. |
|`[MULTIPLE][FITBIT][INFER_FROM_SMARTPHONE_TZ]`| You only need to care about this flag if one or more Fitbit devices sensed data in one or more time zones, and you want RAPIDS to take into account this in its feature computation. Read more in  "How does RAPIDS handle Fitbit devices?" below. |

??? info "Format of `TZCODES_FILE`"
    `TZCODES_FILE` has three columns and a row for each time zone a device visited (a device can be a smartphone or wearable (Fitbit/Empatica)):

    |Column | Description |
    |--|--|
    |`device_id`|A string that uniquely identifies a smartphone or wearable|
    |`tzcode`| A string with the appropriate code from this [list](https://en.wikipedia.org/wiki/List_of_tz_database_time_zones) that represents the time zone where the `device` sensed data|
    |`timestamp`| A UNIX timestamp indicating when was the first time this `device_id` sensed data in `tzcode`|

    ```csv
    device_id,                            tzcode,              timestamp
    13dbc8a3-dae3-4834-823a-4bc96a7d459d, America/New_York,     1587500000000
    13dbc8a3-dae3-4834-823a-4bc96a7d459d, America/Mexico_City,  1587600000000
    13dbc8a3-dae3-4834-823a-4bc96a7d459d, America/Los_Angeles,  1587700000000
    65sa66a5-2d2d-4524-946v-44ascbv4sad7, Europe/Amsterdam,     1587100000000
    65sa66a5-2d2d-4524-946v-44ascbv4sad7, Europe/Berlin,        1587200000000
    65sa66a5-2d2d-4524-946v-44ascbv4sad7, Europe/Amsterdam,     1587300000000
    ```

    Using this file, RAPDIS will create time zone intervals per device, for example for `13dbc8a3-dae3-4834-823a-4bc96a7d459d`:

    -  Interval 1 `[1587500000000, 1587599999999]` for `America/New_York`
    -  Interval 2 `[1587600000000, 1587699999999]` for `America/Mexico_City`
    -  Interval 3 `[1587700000000, now]` for `America/Los_Angeles`

    Any sensor data row from a device will be assigned a timezone if it falls within that interval, for example:

    - A screen row sensed at `1587533333333` will be assigned to `America/New_York` because it falls within Interval 1
    - A screen row sensed at `1587400000000` will be discarded because it was logged outside any interval.
    - **Important** any sensor data that cannot be assigned to a row will be discarded and you will see a warning. To avoid this, add enough time zone intervals. Remember you can set **timestamp** to 0 to avoid having to set an initial timestamp for every device id; see below `What happens if participant X lives in Los Angeles but participant Y lives in Amsterdam and they both stayed there during my study?`

??? note "Can I get the `TZCODES_FILE` from the time zone table collected automatically by the AWARE app?"
    Sure. You can put your timezone table (`timezone.csv`) collected by the AWARE app under `data/external` folder and run:
    ```bash
    python tools/create_multi_timezones_file.py
    ```
    The `TZCODES_FILE` will be saved as `data/external/multiple_timezones.csv`.

??? note "What happens if participant X lives in Los Angeles but participant Y lives in Amsterdam and they both stayed there during my study?"
    Add a row per participant and set timestamp to `0`:
    ```csv
    device_id,                            tzcode,              timestamp
    13dbc8a3-dae3-4834-823a-4bc96a7d459d, America/Los_Angeles,  0
    65sa66a5-2d2d-4524-946v-44ascbv4sad7, Europe/Amsterdam,     0
    ```

??? note "What happens if I forget to add a timezone for one or more devices?"
    It depends on `[IF_MISSING_TZCODE]`. 
    
    If `[IF_MISSING_TZCODE]` is set to `STOP`, RAPIDS will stop its execution and show you an error message.

    If `[IF_MISSING_TZCODE]` is set to `USE_DEFAULT`, it will assign the time zone specified in `[DEFAULT_TZCODE]` to any devices with missing time zone information in `[TZCODES_FILE]`. This is helpful if only a few of your participants had multiple timezones, and you don't want to specify the same time zone for the rest.

??? note "How does RAPIDS handle Fitbit devices?"
    Fitbit devices are not time zone aware, and they always log data with a local date-time string. 

    - When none of the Fitbit devices in your study changed time zones (e.g., `p01` was always in New York and `p02` was always in Amsterdam), you can set a single time zone per Fitbit device id along with a timestamp of 0 (you can still assign multiple time zones to smartphone device ids)
    ```csv
    device_id, tzcode,              timestamp
    fitbit123, America/New_York,     0
    fitbit999, Europe/Amsterdam,     0
    ```

    - On the other hand, when at least one of your Fitbit devices changed time zones **AND** you want RAPIDS to take into account these changes, you need to set `[ALLOW_MULTIPLE_TZ_PER_DEVICE]` to `True`. **You have to manually allow this option because you need to be aware it can produce inaccurate features around the times when time zones changed**. This is because we cannot know precisely when the Fitbit device detected and processed the time zone change.

        If you want to `ALLOW_MULTIPLE_TZ_PER_DEVICE`, you will need to add any time zone changes per device in the `TZCODES_FILE` as explained above. You could obtain this data by hand, but if your participants also used a smartphone during your study, you can use their time zone logs. Recall that in RAPIDS, every participant is represented with a participant file `pXX.yaml`, this file links together multiple devices, and we will use it to know what smartphone time zone data should be applied to Fitbit devices. Thus set `INFER_FROM_SMARTPHONE_TZ` to `TRUE`, if you have included smartphone time zone data in your `TZCODE_FILE` and want to make a participant's Fitbit data time zone aware with their respective smartphone data.

---
## Data Stream Configuration

Modify the following keys in your `config.yaml` depending on the [data stream](../../datastreams/data-streams-introduction) you want to process.

=== "Phone"

    Set `[PHONE_DATA_STREAMS][TYPE]` to the smartphone data stream you want to process (e.g. `aware_mysql`) and configure its parameters (e.g. `[DATABASE_GROUP]`). Ignore the parameters of streams you are not using (e.g. `[FOLDER]` of `aware_csv`).

    ```yaml
    PHONE_DATA_STREAMS:
      USE: aware_mysql

      # AVAILABLE:
      aware_mysql:
        DATABASE_GROUP: MY_GROUP

      aware_csv:
        FOLDER: data/external/aware_csv
    ```

    === "aware_mysql"

        | Key                  | Description                                                                                                                |
        |---------------------|----------------------------------------------------------------------------------------------------------------------------|
        | `[DATABASE_GROUP]`   | A database credentials group. Read the instructions below to set it up    |

        --8<---- "docs/snippets/database.md"

    === "aware_csv"

        | Key                  | Description                                                                                                                |
        |---------------------|----------------------------------------------------------------------------------------------------------------------------|
        | `[FOLDER]`   | Folder where you have to place a CSV file **per** phone sensor. Each file has to contain all the data from every participant you want to process.     |




=== "Fitbit"

    Set `[FITBIT_DATA_STREAMS][TYPE]` to the Fitbit data stream you want to process (e.g. `fitbitjson_mysql`) and configure its parameters (e.g. `[DATABASE_GROUP]`). Ignore the parameters of the other streams you are not using (e.g. `[FOLDER]` of `aware_csv`).

    !!! warning
        You will probably have to tell RAPIDS the name of the columns where you stored your Fitbit data. To do this, modify your chosen stream's `format.yaml` column mappings to match your raw data column names.


    ```yaml
    FITBIT_DATA_STREAMS:
      USE: fitbitjson_mysql

      # AVAILABLE:
      fitbitjson_mysql:
        DATABASE_GROUP: MY_GROUP
        SLEEP_SUMMARY_LAST_NIGHT_END: 660

      fitbitjson_csv:
        FOLDER: data/external/fitbit_csv
        SLEEP_SUMMARY_LAST_NIGHT_END: 660

      fitbitparsed_mysql:
        DATABASE_GROUP: MY_GROUP
        SLEEP_SUMMARY_LAST_NIGHT_END: 660
        
      fitbitparsed_csv:
        FOLDER: data/external/fitbit_csv
        SLEEP_SUMMARY_LAST_NIGHT_END: 660

    ```

    === "fitbitjson_mysql"

        This data stream processes Fitbit data inside a JSON column obtained from the Fitbit API and stored in a MySQL database. Read more about its column mappings and mutations in [`fitbitjson_mysql`](../../datastreams/fitbitjson-mysql#format).


        | Key                  | Description                                                                                                                |
        |---------------------|----------------------------------------------------------------------------------------------------------------------------|
        | `[DATABASE_GROUP]`   | A database credentials group. Read the instructions below to set it up    |
        | `[SLEEP_SUMMARY_LAST_NIGHT_END]`   | Segments are assigned based on this parameter. Any sleep episodes that start between today's SLEEP_SUMMARY_LAST_NIGHT_END (LNE) and tomorrow's LNE are regarded as today's sleep episodes. While today's bedtime is based on today's sleep episodes, today's wake time is based on yesterday's sleep episodes.  |

        --8<---- "docs/snippets/database.md"

    === "fitbitjson_csv"

        This data stream processes Fitbit data inside a JSON column obtained from the Fitbit API and stored in a CSV file. Read more about its column mappings and mutations in [`fitbitjson_csv`](../../datastreams/fitbitjson-csv#format).

        | Key                  | Description                                                                                                                |
        |---------------------|----------------------------------------------------------------------------------------------------------------------------|
        | `[FOLDER]`   | Folder where you have to place a CSV file **per** Fitbit sensor. Each file has to contain all the data from every participant you want to process.     |
        | `[SLEEP_SUMMARY_LAST_NIGHT_END]`   | Segments are assigned based on this parameter. Any sleep episodes that start between today's SLEEP_SUMMARY_LAST_NIGHT_END (LNE) and tomorrow's LNE are regarded as today's sleep episodes. While today's bedtime is based on today's sleep episodes, today's wake time is based on yesterday's sleep episodes.  |


    === "fitbitparsed_mysql"

        This data stream process Fitbit data stored in multiple columns after being parsed from the JSON column returned by Fitbit API and stored in a MySQL database. Read more about its column mappings and mutations in [`fitbitparsed_mysql`](../../datastreams/fitbitparsed-mysql#format).
        

        | Key                  | Description                                                                                                                |
        |---------------------|----------------------------------------------------------------------------------------------------------------------------|
        | `[DATABASE_GROUP]`   | A database credentials group. Read the instructions below to set it up    |
        | `[SLEEP_SUMMARY_LAST_NIGHT_END]`   | Segments are assigned based on this parameter. Any sleep episodes that start between today's SLEEP_SUMMARY_LAST_NIGHT_END (LNE) and tomorrow's LNE are regarded as today's sleep episodes. While today's bedtime is based on today's sleep episodes, today's wake time is based on yesterday's sleep episodes.  |

        --8<---- "docs/snippets/database.md"
        
    === "fitbitparsed_csv"

        This data stream process Fitbit data stored in multiple columns (plain text) after being parsed from the JSON column returned by Fitbit API and stored in a CSV file. Read more about its column mappings and mutations in [`fitbitparsed_csv`](../../datastreams/fitbitparsed-csv#format).

        | Key                  | Description                                                                                                                |
        |---------------------|----------------------------------------------------------------------------------------------------------------------------|
        | `[FOLDER]`   | Folder where you have to place a CSV file **per** Fitbit sensor. Each file has to contain all the data from every participant you want to process.     |
        | `[SLEEP_SUMMARY_LAST_NIGHT_END]`   | Segments are assigned based on this parameter. Any sleep episodes that start between today's SLEEP_SUMMARY_LAST_NIGHT_END (LNE) and tomorrow's LNE are regarded as today's sleep episodes. While today's bedtime is based on today's sleep episodes, today's wake time is based on yesterday's sleep episodes.  |

=== "Empatica"

    Set `[USE]` to the Empatica data stream you want to use; see the table in [introduction to data streams](../../datastreams/data-streams-introduction). Configure any parameters as indicated below.

    ```yaml
    EMPATICA_DATA_STREAMS:
      USE: empatica_zip
      
      # AVAILABLE:
      empatica_zip: 
        FOLDER: data/external/empatica

    ```

    === "empatica_zip"

        | Key             | Description                                                                                                                |
        |---------------------|----------------------------------------------------------------------------------------------------------------------------|
        | `[FOLDER]` | The relative path to a folder containing one subfolder per participant. The name of a participant folder should match their device_id assigned in their participant file. Each participant folder can have one or more zip files with any name; in other words, the sensor data in those zip files belong to a single participant. The zip files are [automatically](https://support.empatica.com/hc/en-us/articles/201608896-Data-export-and-formatting-from-E4-connect-) generated by Empatica and have a CSV file per sensor (`ACC`, `HR`, `TEMP`, `EDA`, `BVP`, `TAGS`). All CSV files of the same type contained in one or more zip files are uncompressed, parsed, sorted by timestamp, and joined together.|

        ??? example "Example of an EMPATICA FOLDER"
            In the file tree below, we want to process three participants' data: `p01`, `p02`, and `p03`. `p01` has two zip files, `p02` has only one zip file, and `p03` has three zip files. Each zip has a CSV file per sensor that is joined together and processed by RAPIDS.

            ```bash
            data/ # this folder exists in the root RAPIDS folder
              external/
                empatica/
                  p01/
                    file1.zip
                    file2.zip
                  p02/
                    aaaa.zip
                  p03/
                    t1.zip
                    t2.zip
                    t3.zip
            ```

---

## Sensor and Features to Process

Finally, you need to modify the `config.yaml` section of the sensors you want to extract behavioral features from. All sensors follow the same naming nomenclature (`DEVICE_SENSOR`) and parameter structure which we explain in the [Behavioral Features Introduction](../../features/feature-introduction/). 

!!! done
    Head over to [Execution](../execution/) to learn how to execute RAPIDS.