library(dplyr)
library(stringr)

unify_ios_screen <- function(ios_screen){
    # In Android we only process UNLOCK to OFF episodes. In iOS we only process UNLOCK to LOCKED episodes,
    # thus, we replace LOCKED with OFF episodes (2 to 0) so we can use Android's code for iOS
    ios_screen <- ios_screen %>% 
        # only keep consecutive pairs of 3,2 events
        filter( (screen_status == 3 & lead(screen_status) == 2) | (screen_status == 2 & lag(screen_status) == 3) ) %>%
        mutate(screen_status = replace(screen_status, screen_status == 2, 0))
    return(ios_screen)
}

unify_ios_battery <- function(ios_battery){
    # We only need to unify battery data for iOS client V1. V2 does it out-of-the-box
    # V1 will not have rows where battery_status is equal to 4
    if(nrow(ios_battery %>% filter(battery_status == 4)) == 0)
        ios_battery <- ios_battery %>%
            mutate(battery_status = replace(battery_status, battery_status == 3, 5),
                battery_status = replace(battery_status, battery_status == 1, 3))
    return(ios_battery)
}

unify_ios_calls <- function(ios_calls){
    # Android’s call types 1=incoming, 2=outgoing, 3=missed
    # iOS' call status 1=incoming, 2=connected, 3=dialing, 4=disconnected
    # iOS' call types based on call status: (1,2,4)=incoming=1, (3,2,4)=outgoing=2, (1,4) or (3,4)=missed=3
    # Sometimes (due to a possible bug in Aware) sequences get logged on the exact same timestamp, thus 3-item sequences can be 2,3,4 or 3,2,4
    # Even tho iOS stores the duration of ringing/dialing for missed calls, we set it to 0 to match Android

    ios_calls <- ios_calls %>%
        arrange(trace, timestamp, call_type) %>% 
        group_by(trace) %>%
                # search for the disconnect event, as it is common to outgoing, received and missed calls
        mutate(completed_call = ifelse(call_type == 4, 2, 0), 
                # assign the same ID to all events before a 4
                completed_call = cumsum(c(1, head(completed_call, -1) != tail(completed_call, -1))), 
                # hack to match ID of last event (4) to that of the previous rows
                completed_call = ifelse(call_type == 4, completed_call - 1, completed_call))

        # We check utc_date_time and local_date_time exist because sometimes we call this function from
        # download_dataset to unify multi-platform participants. At that point such time columns are missing
        if("utc_date_time" %in% colnames(ios_calls) && "local_date_time" %in% colnames(ios_calls)){
            ios_calls <- ios_calls %>% summarise(call_type_sequence = paste(call_type, collapse = ","), # collapse all events before a 4
                        # sanity check, timestamp_diff should be equal or close to duration sum
                        # timestamp_diff = trunc((last(timestamp) - first(timestamp)) / 1000) 
                        # use call_duration = last(call_duration) if you want duration from pick up to hang up
                        # use call_duration = sum(call_duration) if you want duration from dialing/ringing to hang up
                        call_duration = last(call_duration), 
                        timestamp = first(timestamp),
                        utc_date_time = first(utc_date_time),
                        local_date_time = first(local_date_time),
                        local_date = first(local_date),
                        local_time = first(local_time),
                        local_hour = first(local_hour),
                        local_minute = first(local_minute),
                        local_timezone = first(local_timezone),
                        assigned_segments = first(assigned_segments))
        }
        else {
            ios_calls <- ios_calls %>% summarise(call_type_sequence = paste(call_type, collapse = ","), call_duration = sum(call_duration),  timestamp = first(timestamp))
        }
        ios_calls <- ios_calls %>% mutate(call_type = case_when(
            call_type_sequence == "1,2,4" | call_type_sequence == "2,1,4" ~ 1, # incoming
            call_type_sequence == "1,4" ~ 3, # missed
            call_type_sequence == "3,2,4" | call_type_sequence == "2,3,4" ~ 2, # outgoing
            call_type_sequence == "3,4" ~ 4, # outgoing missed, we create this temp missed state to assign a duration of 0 below
            TRUE ~ -1), # other, call sequences without a disconnect (4) event are discarded
            # assign a duration of 0 to incoming and outgoing missed calls
            call_duration = ifelse(call_type == 3 | call_type == 4, 0, call_duration), 
            # get rid of the temp missed call type, set to 2 to match Android. See https://github.com/carissalow/rapids/issues/79
            call_type = ifelse(call_type == 4, 2, call_type) 
        ) %>% 
        # discard sequences without an event 4 (disconnect)
        filter(call_type > 0) %>%
        ungroup() %>%
        arrange(timestamp)

    return(ios_calls)
}

clean_ios_activity_column <- function(ios_gar){
    ios_gar <- ios_gar %>%
        mutate(activities = str_replace_all(activities, pattern = '("|\\[|\\])', replacement = ""))

    existent_multiple_activities <- ios_gar %>%
        filter(str_detect(activities, ",")) %>% 
        group_by(activities) %>%
        summarise(mutiple_activities = unique(activities)) %>% 
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

unify_ios_gar <- function(ios_gar){
    # We only need to unify Google Activity Recognition data for iOS
    # discard rows where activities column is blank
    ios_gar <- ios_gar[-which(ios_gar$activities == ""), ]
    # clean "activities" column of ios_gar
    ios_gar <- clean_ios_activity_column(ios_gar)

    # make it compatible with android version: generate "activity_name" and "activity_type" columns
    ios_gar  <-  ios_gar %>% 
        mutate(activity_name = case_when(activities == "automotive" ~ "in_vehicle",
                                         activities == "cycling" ~ "on_bicycle",
                                         activities == "walking" | activities == "running" ~ "on_foot",
                                         activities == "stationary" ~ "still"),
               activity_type = case_when(activities == "automotive" ~ 0,
                                         activities == "cycling" ~ 1,
                                         activities == "walking" | activities == "running" ~ 2,
                                         activities == "stationary" ~ 3,
                                         activities == "unknown" ~ 4))
    
    return(ios_gar)
}

unify_ios_conversation <- function(conversation){
    if(nrow(conversation) > 0){
        duration_check <- conversation %>% 
            select(double_convo_start, double_convo_end) %>% 
            mutate(start_is_seconds = double_convo_start <= 9999999999,
                end_is_seconds = double_convo_end <= 9999999999) # Values smaller than 9999999999 are in seconds instead of milliseconds
        start_end_in_seconds = sum(duration_check$start_is_seconds) + sum(duration_check$end_is_seconds)

        if(start_end_in_seconds > 0) # convert seconds to milliseconds
            conversation <- conversation %>% mutate(double_convo_start = double_convo_start * 1000, double_convo_end = double_convo_end * 1000)
    }
    return(conversation)
}

# This function is used in download_dataset.R
unify_raw_data <- function(dbEngine, table, start_datetime_utc, end_datetime_utc, aware_multiplatform_tables, unifiable_tables, device_ids, platforms){
  # If platforms is 'multiple', fetch each device_id's platform from aware_device, otherwise, use those given by the user
  if(length(platforms) == 1 && platforms == "multiple")
      devices_platforms <- dbGetQuery(dbEngine, paste0("SELECT device_id,brand FROM aware_device WHERE device_id IN ('", paste0(device_ids, collapse = "','"), "')")) %>% 
        mutate(platform = ifelse(brand == "iPhone", "ios", "android"))
    else
      devices_platforms <- data.frame(device_id = device_ids, platform = platforms)

  # Get existent tables in database
  available_tables_in_db <- dbGetQuery(dbEngine, paste0("SELECT table_name FROM information_schema.tables WHERE table_type = 'base table' AND table_schema='", dbGetInfo(dbEngine)$dbname,"'")) %>% pull(table_name)
  
  # Parse the table names for activity recognition and conversation plugins because they are different between android and ios
  ar_tables <- setNames(aware_multiplatform_tables[1:2], c("android", "ios"))
  conversation_tables <- setNames(aware_multiplatform_tables[3:4], c("android", "ios"))

  participants_sensordata <- list()
  for(i in 1:nrow(devices_platforms)) {
    row <- devices_platforms[i,]
    device_id <- row$device_id
    platform <- row$platform
    
    # Handle special cases when tables for the same sensor have different names for Android and iOS (AR and conversation)
    if(table %in% ar_tables)
      table <- ar_tables[[platform]]
    else if(table %in% conversation_tables)
      table <- conversation_tables[[platform]]

    if(table %in% available_tables_in_db){
      query <- paste0("SELECT * FROM ", table, " WHERE device_id IN ('", device_id, "')")
      if("timestamp" %in% available_columns && !(is.na(start_datetime_utc)) && !(is.na(end_datetime_utc)) && start_datetime_utc < end_datetime_utc){
        query <- paste0(query, "AND timestamp BETWEEN 1000*UNIX_TIMESTAMP('", start_datetime_utc, "') AND 1000*UNIX_TIMESTAMP('", end_datetime_utc, "')")
      }
      sensor_data <- unify_data(dbGetQuery(dbEngine, query), table, platform, unifiable_tables)
      participants_sensordata <- append(participants_sensordata, list(sensor_data))
    }else{
      warning(paste0("Missing ", table, " table. We unified the data from ", paste0(devices_platforms$device_id, collapse = " and "), " but without records from this missing table for ", device_id))
    }
  }
  unified_data <- bind_rows(participants_sensordata)
  return(unified_data)

}

# This function is used in unify_ios_android.R and unify_raw_data function
unify_data <- function(sensor_data, sensor, platform, unifiable_sensors){
    if(sensor == unifiable_sensors$calls){
        if(platform == "ios"){
            sensor_data = unify_ios_calls(sensor_data)
        }
        # android calls remain unchanged
    } else if(sensor == unifiable_sensors$battery){
        if(platform == "ios"){
            sensor_data = unify_ios_battery(sensor_data)
        }
        # android battery remains unchanged
    } else if(sensor == unifiable_sensors$ios_activity_recognition){
        sensor_data = unify_ios_gar(sensor_data)
    } else if(sensor == unifiable_sensors$screen){
        if(platform == "ios"){
            sensor_data = unify_ios_screen(sensor_data)
        }
        # android screen remains unchanged
    } else if(sensor == unifiable_sensors$ios_conversation){
        sensor_data = unify_ios_conversation(sensor_data)
    }
    return(sensor_data)
}