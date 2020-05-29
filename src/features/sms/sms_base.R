filter_by_day_segment <- function(data, day_segment) {
  if(day_segment %in% c("morning", "afternoon", "evening", "night"))
    data <- data %>% filter(local_day_segment == day_segment)
  return(data)
}

base_sms_features <- function(sms, sms_type, day_segment, requested_features){
    # Output dataframe
    features = data.frame(local_date = character(), stringsAsFactors = FALSE)

    # The name of the features this function can compute
    base_features_names  <- c("countmostfrequentcontact", "count", "distinctcontacts", "timefirstsms", "timelastsms")

    # The subset of requested features this function can compute
    features_to_compute  <- intersect(base_features_names, requested_features)

    # Filter rows that belong to the sms type and day segment of interest
    sms <- sms %>% filter(message_type == ifelse(sms_type == "received", "1", "2")) %>% 
        filter_by_day_segment(day_segment)

    # If there are not features or data to work with, return an empty df with appropiate columns names
    if(length(features_to_compute) == 0)
        return(features)
    if(nrow(sms) < 1)
        return(cbind(features, read.csv(text = paste(paste("sms", sms_type, day_segment, features_to_compute, sep = "_"), collapse = ","), stringsAsFactors = FALSE)))

    for(feature_name in features_to_compute){
        if(feature_name == "countmostfrequentcontact"){
            # Get the number of messages for the most frequent contact throughout the study
            feature <- sms %>% group_by(trace) %>% 
                mutate(N=n()) %>% 
                ungroup() %>%
                filter(N == max(N)) %>% 
                head(1) %>% # if there are multiple contacts with the same amount of messages pick the first one only
                group_by(local_date) %>% 
                summarise(!!paste("sms", sms_type, day_segment, feature_name, sep = "_") := n())  %>% 
                replace(is.na(.), 0)

            features <- merge(features, feature, by="local_date", all = TRUE)
        } else {
            feature <- sms %>% 
                group_by(local_date)
            
            feature <- switch(feature_name,
                    "count" = feature %>% summarise(!!paste("sms", sms_type, day_segment, feature_name, sep = "_") := n()),
                    "distinctcontacts" = feature %>% summarise(!!paste("sms", sms_type, day_segment, feature_name, sep = "_") := n_distinct(trace)),
                    "timefirstsms" = feature %>% summarise(!!paste("sms", sms_type, day_segment, feature_name, sep = "_") := first(local_hour) * 60 + first(local_minute)),
                    "timelastsms" = feature %>% summarise(!!paste("sms", sms_type, day_segment, feature_name, sep = "_") := last(local_hour) * 60 + last(local_minute)))

            features <- merge(features, feature, by="local_date", all = TRUE)
        }
    }

    return(features)
}