library('tidyr') 
library('stringr')

message_features_of_type <- function(messages, messages_type, time_segment, requested_features){
    # Output dataframe
    features = data.frame(local_segment = character(), stringsAsFactors = FALSE)

    # The name of the features this function can compute
    base_features_names  <- c("countmostfrequentcontact", "count", "distinctcontacts", "timefirstmessage", "timelastmessage")

    # The subset of requested features this function can compute
    features_to_compute  <- intersect(base_features_names, requested_features)
 
    # If there are not features or data to work with, return an empty df with appropiate columns names
    if(length(features_to_compute) == 0)
        return(features)
    if(nrow(messages) < 1)
        return(cbind(features, read.csv(text = paste(paste(messages_type, features_to_compute, sep = "_"), collapse = ","), stringsAsFactors = FALSE)))

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
                group_by(local_segment) %>% 
                summarise(!!paste(messages_type, feature_name, sep = "_") := sum(trace == mostfrequentcontact))
            features <- merge(features, feature, by="local_segment", all = TRUE)
        } else {
            feature <- messages %>% 
                group_by(local_segment)
            
            feature <- switch(feature_name,
                    "count" = feature %>% summarise(!!paste(messages_type, feature_name, sep = "_") := n()),
                    "distinctcontacts" = feature %>% summarise(!!paste(messages_type, feature_name, sep = "_") := n_distinct(trace)),
                    "timefirstmessage" = feature %>% summarise(!!paste(messages_type, feature_name, sep = "_") := first(local_hour) * 60 + first(local_minute)),
                    "timelastmessage" = feature %>% summarise(!!paste(messages_type, feature_name, sep = "_") := last(local_hour) * 60 + last(local_minute)))

            features <- merge(features, feature, by="local_segment", all = TRUE)
        }
    }
    return(features)
}

rapids_features <- function(sensor_data_files, time_segment, provider){
    messages_data <-  read.csv(sensor_data_files[["sensor_data"]], stringsAsFactors = FALSE)
    messages_data <- messages_data %>% filter_data_by_segment(time_segment)
    messages_types = provider[["MESSAGES_TYPES"]]
    messages_features <- setNames(data.frame(matrix(ncol=1, nrow=0)), c("local_segment"))

    for(message_type in messages_types){
        # Filter rows that belong to the message type and time segment of interest
        message_type_label = ifelse(message_type == "received", "1", ifelse(message_type == "sent", "2", NA))
        if(is.na(message_type_label))
            stop(paste("Message type can online be received or sent but instead you typed: ", message_type, " in config[PHONE_MESSAGES][MESSAGES_TYPES]"))

        requested_features <- provider[["FEATURES"]][[message_type]]
        messages_of_type <- messages_data %>% filter(message_type == message_type_label)

        features <- message_features_of_type(messages_of_type, message_type, time_segment, requested_features)
        messages_features <- merge(messages_features, features, all=TRUE)
    }
    messages_features <- messages_features %>% mutate_at(vars(contains("countmostfrequentcontact") | contains("distinctcontacts") | contains("count")), list( ~ replace_na(., 0)))
    return(messages_features)
}