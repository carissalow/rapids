source("renv/activate.R")

library(dplyr)
library(stringr)

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
                completed_call = ifelse(call_type == 4, completed_call - 1, completed_call)) %>% 
        summarise(call_type_sequence = paste(call_type, collapse = ","), # collapse all events before a 4
                    # use this if Android measures calls' duration from pick up to hang up
                    # duration = last(call_duration), 
                    # sanity check, timestamp_diff should be equal or close to duration sum
                    # timestamp_diff = trunc((last(timestamp) - first(timestamp)) / 1000) 
                    # use this if Android measures calls' duration from dialing/ringing to hang up
                    call_duration = sum(call_duration), 

                    timestamp = first(timestamp),
                    utc_date_time = first(utc_date_time),
                    local_date_time = first(local_date_time),
                    local_date = first(local_date),
                    local_time = first(local_time),
                    local_hour = first(local_hour),
                    local_minute = first(local_minute),
                    local_day_segment = first(local_day_segment)
                    ) %>% 
        mutate(call_type = case_when(
            call_type_sequence == "1,2,4" | call_type_sequence == "2,1,4" ~ 1, # incoming
            call_type_sequence == "1,4" ~ 3, # missed
            call_type_sequence == "3,2,4" | call_type_sequence == "2,3,4" ~ 2, # outgoing
            call_type_sequence == "3,4" ~ 4, # outgoing missed, we create this temp missed state to assign a duration of 0 below
            TRUE ~ -1), # other, call sequences without a disconnect (4) event are discarded
            # assign a duration of 0 to incoming and outgoing missed calls
            call_duration = ifelse(call_type == 3 | call_type == 4, 0, call_duration), 
            # get rid of the temp missed call type, set to 3 to match Android
            call_type = ifelse(call_type == 4, 3, call_type) 
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


sensor_data <- read.csv(snakemake@input[["sensor_data"]], stringsAsFactors = FALSE)
participant_info <- snakemake@input[["participant_info"]]
sensor <- snakemake@params[["sensor"]]
platform <- readLines(participant_info, n=2)[[2]]

if(sensor == "calls"){
    if(platform == "ios"){
        sensor_data = unify_ios_calls(sensor_data)
    }
    # android calls remain unchanged
} else if(sensor == "battery"){
    if(platform == "ios"){
        sensor_data = unify_ios_battery(sensor_data)
    }
    # android battery remains unchanged
} else if(sensor == "plugin_ios_activity_recognition"){
    sensor_data = unify_ios_gar(sensor_data)
}
write.csv(sensor_data, snakemake@output[[1]], row.names = FALSE)
