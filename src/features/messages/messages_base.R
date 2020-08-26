library('tidyr') 
library('stringr')

base_messages_features <- function(messages, messages_type, day_segment, requested_features){
    # Output dataframe
    features = data.frame(local_segment = character(), stringsAsFactors = FALSE)

    # The name of the features this function can compute
    base_features_names  <- c("countmostfrequentcontact", "count", "distinctcontacts", "timefirstmessage", "timelastmessage")

    # The subset of requested features this function can compute
    features_to_compute  <- intersect(base_features_names, requested_features)

    # Filter the rows that belong to day_segment, and put the segment full name in a new column for grouping
    date_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2}"
    hour_regex = "[0-9]{2}:[0-9]{2}:[0-9]{2}"
    messages <- messages %>% 
        filter(message_type == ifelse(messages_type == "received", "1", ifelse(messages_type == "sent", 2, NA))) %>% 
        filter(grepl(paste0("\\[", day_segment, "#"),assigned_segments)) %>% 
        mutate(local_segment = str_extract(assigned_segments, paste0("\\[", day_segment, "#", date_regex, "#", hour_regex, "#", date_regex, "#", hour_regex, "\\]")),
                local_segment = str_sub(local_segment, 2, -2)) # get rid of first and last character([])
 
    # If there are not features or data to work with, return an empty df with appropiate columns names
    if(length(features_to_compute) == 0)
        return(features)
    if(nrow(messages) < 1)
        return(cbind(features, read.csv(text = paste(paste("messages", messages_type, features_to_compute, sep = "_"), collapse = ","), stringsAsFactors = FALSE)))

    for(feature_name in features_to_compute){
        if(feature_name == "countmostfrequentcontact"){
            # Get the number of messages for the most frequent contact throughout the study
            mostfrequentcontact <- messages %>% 
                group_by(trace) %>% 
                mutate(N=n()) %>% 
                ungroup() %>%
                filter(N == max(N)) %>% 
                head(1) %>% # if there are multiple contacts with the same amount of messages pick the first one only
                pull(trace)
            feature <- messages %>% 
                filter(trace == mostfrequentcontact) %>% 
                group_by(local_segment) %>% 
                summarise(!!paste("messages", messages_type, feature_name, sep = "_") := n())  %>% 
                replace(is.na(.), 0)
            features <- merge(features, feature, by="local_segment", all = TRUE)
        } else {
            feature <- messages %>% 
                group_by(local_segment)
            
            feature <- switch(feature_name,
                    "count" = feature %>% summarise(!!paste("messages", messages_type, feature_name, sep = "_") := n()),
                    "distinctcontacts" = feature %>% summarise(!!paste("messages", messages_type, feature_name, sep = "_") := n_distinct(trace)),
                    "timefirstmessage" = feature %>% summarise(!!paste("messages", messages_type, feature_name, sep = "_") := first(local_hour) * 60 + first(local_minute)),
                    "timelastmessage" = feature %>% summarise(!!paste("messages", messages_type, feature_name, sep = "_") := last(local_hour) * 60 + last(local_minute)))

            features <- merge(features, feature, by="local_segment", all = TRUE)
        }
    }
    features <- features %>% mutate_at(vars(contains("countmostfrequentcontact")), list( ~ replace_na(., 0)))
    return(features)
}