source("renv/activate.R")
library("dplyr", warn.conflicts = F)


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
            ios_calls <- ios_calls %>% summarise(call_type_sequence = paste(call_type, collapse = ","), call_duration = sum(as.numeric(call_duration)),  timestamp = first(timestamp), device_id = first(device_id))
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
    
    ios_calls <- select(ios_calls, -call_type_sequence)

    return(ios_calls)
}

main <- function(data, stream_parameters){
    return(unify_ios_calls(data))
}