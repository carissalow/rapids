# Add New Data Streams

A data stream is a set of sensor data collected using a specific type of **device** with a specific **format** and stored in a specific **container**. RAPIDS is agnostic to data streams' formats and containers (see the [Data Streams Introduction](../data-streams-introduction) for a list of supported data streams).

In short, a format describes how raw data logged by a device maps to the data expected by RAPIDS, and a container is a script that connects to the database or file where that data is stored. 

The most common cases when you would want to implement a new data stream are:

- You collected data with a mobile sensing app RAPIDS does not support yet. For example, [Beiwe](https://www.beiwe.org/) data stored in MySQL. You will need to define a new format and a new container.
- You collected data with a mobile sensing app RAPIDS supports but this data is stored in a container that RAPIDS can't connect to yet. For example, AWARE data stored in PostgreSQL. In this case, you can reuse the format of the`aware_mysql` stream but you will need to implement a new container script.

!!! hint
    RAPIDS supports smartphones, Fitbit, and Empatica devices, you can add a new data stream for the first two.

## Formats and Containers in RAPIDS

**CONTAINER**. The container of a data stream is queried using a `container.R` script. This script implements functions that will pull data from a database, file, etc.

**FORMAT**. The format of a data stream is described using a `format.yaml` file. A format file describes the mapping between your stream's raw data and the data that RAPIDS needs.

Both the `container.R` and the `format.yaml` are saved under `src/data/streams/[stream_name]` where `[stream_name]` can be 
`aware_mysql` for example.

## Implement a Container

The `container.R` script of a data stream should be implemented in R. This script must have two functions if you are implementing a stream for phone data, or one function otherwise.

=== "download_data"

    This function returns the data columns for a specific sensor and participant. It has the following parameters:

    | Param              | Description                                                                                           |   
    |--------------------|-------------------------------------------------------------------------------------------------------|
    | data_configuration | Any parameters (keys/values) set by the user in any `[DEVICE_DATA_STREAMS][stream_name]` key of `config.yaml`. For example, `[DATABASE_GROUP]` inside `[FITBIT_DATA_STREAMS][fitbitjson_mysql]` | 
    | sensor_container   | The value set by the user in any `[DEVICE_SENSOR][CONTAINER]` key of `config.yaml`. It can be a table, file path, or whatever data source you want to support that contains the **data from a single sensor for all participants**. For example, `[PHONE_ACCELEROMETER][CONTAINER]`|
    | device             | The device id that you need to get the data for (this is set by the user in the [participant files](../../setup/configuration/#participant-files)). For example, in AWARE this device is a uuid|
    | columns            | A list of the columns that you need to get from `sensor_container`. You specify these columns in your stream's `format.yaml`|


    !!! example
        This is the `download_data` function we implemented for `aware_mysql`. Note that we can `message`, `warn` or `stop` the user during execution.

        ```r
        download_data <- function(data_configuration, device, sensor_container, columns){
            group <- data_configuration$SOURCE$DATABASE_GROUP
            dbEngine <- dbConnect(MariaDB(), default.file = "./.env", group = group)
            
            
            query <- paste0("SELECT ", paste(columns, collapse = ",")," FROM ", sensor_container, " WHERE device_id = '", device,"'")
            # Letting the user know what we are doing
            message(paste0("Executing the following query to download data: ", query)) 
            sensor_data <- dbGetQuery(dbEngine, query)
            
            dbDisconnect(dbEngine)
            
            if(nrow(sensor_data) == 0)
                warning(paste("The device '", device,"' did not have data in ", sensor_container))

            return(sensor_data)
        }
        ```

=== "infer_device_os"

    !!! warning
        This function is only necessary for phone data streams. 
    
    RAPIDS allows users to use the keyword `infer` (previously `multiple`) to [automatically infer](../../setup/configuration/#structure-of-participants-files) the mobile Operative System a device (phone) was running. 
    
    If you have a way to infer the OS of a device id, implement this function. For example, for AWARE data we use the `aware_device` table.
 
    If you don't have a way to infer the OS, call `stop("Error Message")` so other users know they can't use `infer` or the inference failed, and they have to assign the OS manually in the participant file.
    
    This function returns the operative system (`android` or `ios`) for a specific device. It has the following parameters:

    | Param              | Description                                                                                           |   
    |--------------------|-------------------------------------------------------------------------------------------------------|
    | data_configuration | Any parameters (keys/values) set by the user in any `[DEVICE_DATA_STREAMS][stream_name]` key of `config.yaml`. For example, `[DATABASE_GROUP]` inside `[FITBIT_DATA_STREAMS][fitbitjson_mysql]` | 
    | device             | The device id that you need to infer the OS for (this is set by the user in the [participant files](../../setup/configuration/#participant-files)). For example, in AWARE this device is a uuid|


    !!! example
        This is the `infer_device_os` function we implemented for `aware_mysql`. Note that we can `message`, `warn` or `stop` the user during execution.

        ```r
        infer_device_os <- function(data_configuration, device){
            group <- data_configuration$SOURCE$DATABASE_GROUP # specified DB credentials group in config.yaml
            
            dbEngine <- dbConnect(MariaDB(), default.file = "./.env", group = group)
            query <- paste0("SELECT device_id,brand FROM aware_device WHERE device_id = '", device, "'")
            message(paste0("Executing the following query to infer phone OS: ", query)) 
            os <- dbGetQuery(dbEngine, query)
            dbDisconnect(dbEngine)
            
            if(nrow(os) > 0)
                return(os %>% mutate(os = ifelse(brand == "iPhone", "ios", "android")) %>% pull(os))
            else
                stop(paste("We cannot infer the OS of the following device id because it does not exist in the aware_device table:", device))
            
            return(os)
        }
        ```

## Implement a Format

A format describes the mapping between your stream's raw data and the data that RAPIDS needs. This file has a section per sensor (e.g. `PHONE_ACCELEROMETER`), and each section has two keys (attributes):

1. `COLUMN_MAPPINGS` is a mapping between the columns RAPIDS needs and the columns your raw data has. 
2. `MUTATION_SCRIPTS` are a collection of R or Python scripts that transform your raw data into the format RAPIDS needs. 

Let's explain these keys with examples.

### Name mapping

The mapping for some sensors is straightforward. For example, accelerometer data most of the time has a timestamp, three axis (x,y,z) and a device id that produced it. It is likely that AWARE and a different sensing app like Beiwe logged accelerometer data in the same way but with different columns names. In this case we only need to match Beiwe data columns to RAPIDS columns one-to-one:

```yaml hl_lines="4 5 6 7 8"
PHONE_ACCELEROMETER:
  ANDROID:
    COLUMN_MAPPINGS:
      TIMESTAMP: beiwe_timestamp
      DEVICE_ID: beiwe_deviceID
      DOUBLE_VALUES_0: beiwe_x
      DOUBLE_VALUES_1: beiwe_y
      DOUBLE_VALUES_2: beiwe_z
    MUTATION_SCRIPTS: # it's ok if this is empty
```

### Value mapping
For some sensors we need to map column names and values. For example, screen data has ON and OFF events, let's suppose Beiwe represents an ON event with the number `1` but RAPIDS identifies ON events with the number `2`. In this case we need to mutate the raw data coming from Beiwe and replace all `1`s with `2`s.

We do this by listing one or more R or Python scripts in `MUTATION_SCRIPTS` that will be executed in order:

```yaml hl_lines="8"
PHONE_SCREEN:
  ANDROID:
    COLUMN_MAPPINGS:
      TIMESTAMP: beiwe_timestamp
      DEVICE_ID: beiwe_deviceID
      EVENT: beiwe_event
    MUTATION_SCRIPTS:
        - src/data/streams/mutations/phone/beiwe/beiwe_screen_map.py
```

Every `MUTATION_SCRIPT` has a `main` function that receives a data frame with your raw sensor data and should return the mutated data. We usually store all mutation scripts under `src/data/streams/mutations/[device]/[platform]/` and they can be reused across data streams.

!!! hint
    This `MUTATION_SCRIPT` can also be used to clean/preprocess your data before extracting behavioral features.

=== "python"
    Example of a python mutation script
    ```python
    import pandas as pd

    def main(data):
        # mutate data
        return(data)
    ```
=== "R"
    Example of a R mutation script
    ```r
    source("renv/activate.R") # needed to use RAPIDS renv environment
    library(dplyr)

    main <- function(data){
        # mutate data
        return(data)
    }
    ```

### Complex mapping
Sometimes, your raw data doesn't even have the same columns RAPIDS expects for a sensor. For example, let's pretend Beiwe stores `PHONE_ACCELEROMETER` axis data in a single column called `acc_col` instead of three: `x-y-z`. You need to create a `MUTATION_SCRIPT` to split `acc_col` into three columns `x`, `y`, and `z`. 

For this, you mark the missing `COLUMN_MAPPINGS` with the word `FLAG_TO_MUTATE`, map `acc_col` to `FLAG_AS_EXTRA`, and list a Python script under `MUTATION_SCRIPT` with the code to split `acc_col`.

Every column mapped with `FLAG_AS_EXTRA` will be included in the data frame you receive in your mutation script and we recommend deleting them from the returned data frame after they are not needed anymore.

!!! hint
    Note that although `COLUMN_MAPPINGS` keys are in capital letters for readability (e.g. `DOUBLE_VALUES_0`), the names of the final columns you mutate in your scripts should be lower case.

```yaml hl_lines="6 7 8 9 11"
PHONE_ACCELEROMETER:
  ANDROID:
    COLUMN_MAPPINGS:
      TIMESTAMP: beiwe_timestamp
      DEVICE_ID: beiwe_deviceID
      DOUBLE_VALUES_0: FLAG_TO_MUTATE
      DOUBLE_VALUES_1: FLAG_TO_MUTATE
      DOUBLE_VALUES_2: FLAG_TO_MUTATE
      FLAG_AS_EXTRA: acc_col
    MUTATION_SCRIPTS:
        - src/data/streams/mutations/phone/beiwe/beiwe_split_acc.py
```

This is a draft of `beiwe_split_acc.py` `MUTATION_SCRIPT`:
```python
import pandas as pd

def main(data):
    # data has the acc_col
    # split acc_col into three columns: double_values_0, double_values_1, double_values_2 to match RAPIDS format
    # remove acc_col since we don't need it anymore
    return(data)
```

### OS complex mapping
There is a special case for a complex mapping scenario for smartphone data streams. The Android and iOS sensor APIs return data in different formats for certain sensors (like screen, activity recognition, battery, among others). 

In case you didn't notice, the examples we have used so far are grouped under an `ANDROID` key, which means they will be applied to data collected by Android phones. Additionally, each sensor has an `IOS` key for a similar purpose. We use the complex mapping described above to transform iOS data into an Android format (it's always iOS to Android and any new phone data stream must do the same).

For example, this is the `format.yaml` key for `PHONE_ACTVITY_RECOGNITION`. Note that the `ANDROID` mapping is simple (one-to-one) but the `IOS` mapping is complex with two `FLAG_TO_MUTATE` columns, one `FLAG_AS_EXTRA` column, and one `MUTATION_SCRIPT`.

```yaml hl_lines="14 15 17 19"
PHONE_ACTIVITY_RECOGNITION:
  ANDROID:
    COLUMN_MAPPINGS:
      TIMESTAMP: timestamp
      DEVICE_ID: device_id
      ACTIVITY_TYPE: activity_type
      ACTIVITY_NAME: activity_name
      CONFIDENCE: confidence
    MUTATION_SCRIPTS: 
  IOS:
    COLUMN_MAPPINGS:
      TIMESTAMP: timestamp
      DEVICE_ID: device_id
      ACTIVITY_TYPE: FLAG_TO_MUTATE
      ACTIVITY_NAME: FLAG_TO_MUTATE
      CONFIDENCE: confidence
      FLAG_AS_EXTRA: activities
    MUTATION_SCRIPTS:
      - "src/data/streams/mutations/phone/aware/activity_recogniton_ios_unification.R"
```

??? "Example activity_recogniton_ios_unification.R"
    In this `MUTATION_SCRIPT` we create `ACTIVITY_NAME` and `ACTIVITY_TYPE` based on `activities`, and map `confidence` iOS values to Android values.
    ```R
    source("renv/activate.R")
    library("dplyr", warn.conflicts = F)
    library(stringr)

    clean_ios_activity_column <- function(ios_gar){
        ios_gar <- ios_gar %>%
            mutate(activities = str_replace_all(activities, pattern = '("|\\[|\\])', replacement = ""))

        existent_multiple_activities <- ios_gar %>%
            filter(str_detect(activities, ",")) %>% 
            group_by(activities) %>%
            summarise(mutiple_activities = unique(activities), .groups = "drop_last") %>% 
            pull(mutiple_activities)

        known_multiple_activities <- c("stationary,automotive")
        unkown_multiple_actvities <- setdiff(existent_multiple_activities, known_multiple_activities)
        if(length(unkown_multiple_actvities) > 0){
            stop(paste0("There are unkwown combinations of ios activities, you need to implement the decision of the ones to keep: ", unkown_multiple_actvities))
        }

        ios_gar <- ios_gar %>%
            mutate(activities = str_replace_all(activities, pattern = "stationary,automotive", replacement = "automotive"))
        
        return(ios_gar)
    }

    unify_ios_activity_recognition <- function(ios_gar){
        # We only need to unify Google Activity Recognition data for iOS
        # discard rows where activities column is blank
        ios_gar <- ios_gar[-which(ios_gar$activities == ""), ]
        # clean "activities" column of ios_gar
        ios_gar <- clean_ios_activity_column(ios_gar)

        # make it compatible with android version: generate "activity_name" and "activity_type" columns
        ios_gar  <-  ios_gar %>% 
            mutate(activity_name = case_when(activities == "automotive" ~ "in_vehicle",
                                            activities == "cycling" ~ "on_bicycle",
                                            activities == "walking" ~ "walking",
                                            activities == "running" ~ "running",
                                            activities == "stationary" ~ "still"),
                    activity_type = case_when(activities == "automotive" ~ 0,
                                            activities == "cycling" ~ 1,
                                            activities == "walking" ~ 7,
                                            activities == "running" ~ 8,
                                            activities == "stationary" ~ 3,
                                            activities == "unknown" ~ 4),
                    confidence = case_when(confidence == 0 ~ 0,
                                        confidence == 1 ~ 50,
                                        confidence == 2 ~ 100)
                                        ) %>% 
            select(-activities)
        
        return(ios_gar)
    }

    main <- function(data){
        return(unify_ios_activity_recognition(data))
    }
    ```
