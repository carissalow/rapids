library('tidyr')

filter_by_day_segment <- function(data, day_segment) {
  if(day_segment %in% c("morning", "afternoon", "evening", "night"))
    data <- data %>% filter(local_day_segment == day_segment)
  else if(day_segment == "daily")
    return(data)
  else 
    return(data %>% head(0))
}

Mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

base_call_features <- function(calls, call_type, day_segment, requested_features){
    # Output dataframe
    features = data.frame(local_date = character(), stringsAsFactors = FALSE)

    # The name of the features this function can compute
    base_features_names  <- c("count", "distinctcontacts", "meanduration", "sumduration", "minduration", "maxduration", "stdduration", "modeduration", "entropyduration", "timefirstcall", "timelastcall", "countmostfrequentcontact")

    # The subset of requested features this function can compute
    features_to_compute  <- intersect(base_features_names, requested_features)

    # Filter rows that belong to the calls type and day segment of interest
    call_type_label = ifelse(call_type == "incoming", "1", ifelse(call_type == "outgoing", "2", ifelse(call_type == "missed", "3", NA)))
    if(is.na(call_type_label))
        stop(paste("Call type can online be incoming, outgoing or missed but instead you typed: ", call_type))
    calls <- calls %>% filter(call_type == call_type_label) %>% filter_by_day_segment(day_segment)

    # If there are not features or data to work with, return an empty df with appropiate columns names
    if(length(features_to_compute) == 0)
        return(features)
    if(nrow(calls) < 1)
        return(cbind(features, read.csv(text = paste(paste("call", call_type, day_segment, features_to_compute, sep = "_"), collapse = ","), stringsAsFactors = FALSE)))

    for(feature_name in features_to_compute){
        if(feature_name == "countmostfrequentcontact"){
            # Get the number of messages for the most frequent contact throughout the study
            mostfrequentcontact <- calls %>% 
                group_by(trace) %>% 
                mutate(N=n()) %>% 
                ungroup() %>%
                filter(N == max(N)) %>% 
                head(1) %>% # if there are multiple contacts with the same amount of messages pick the first one only
                pull(trace)
            feature <- calls %>% 
                filter(trace == mostfrequentcontact) %>% 
                group_by(local_date) %>% 
                summarise(!!paste("call", call_type, day_segment, feature_name, sep = "_") := n())  %>% 
                replace(is.na(.), 0)
            features <- merge(features, feature, by="local_date", all = TRUE)
        } else {
            feature <- calls %>% 
                group_by(local_date)

            feature <- switch(feature_name,
                "count" = feature %>% summarise(!!paste("call", call_type, day_segment, feature_name, sep = "_") := n()),
                "distinctcontacts" = feature %>% summarise(!!paste("call", call_type, day_segment, feature_name, sep = "_") := n_distinct(trace)),
                "meanduration" = feature %>% summarise(!!paste("call", call_type, day_segment, feature_name, sep = "_") := mean(call_duration)),
                "sumduration" = feature %>% summarise(!!paste("call", call_type, day_segment, feature_name, sep = "_") := sum(call_duration)),
                "minduration" = feature %>% summarise(!!paste("call", call_type, day_segment, feature_name, sep = "_") := min(call_duration)),
                "maxduration" = feature %>% summarise(!!paste("call", call_type, day_segment, feature_name, sep = "_") := max(call_duration)),
                "stdduration" = feature %>% summarise(!!paste("call", call_type, day_segment, feature_name, sep = "_") := sd(call_duration)),
                "modeduration" = feature %>% summarise(!!paste("call", call_type, day_segment, feature_name, sep = "_") := Mode(call_duration)),
                "entropyduration" = feature %>% summarise(!!paste("call", call_type, day_segment, feature_name, sep = "_") := entropy.MillerMadow(call_duration)),
                "timefirstcall" = feature %>% summarise(!!paste("call", call_type, day_segment, feature_name, sep = "_") := first(local_hour) * 60 + first(local_minute)),
                "timelastcall" = feature %>% summarise(!!paste("call", call_type, day_segment, feature_name, sep = "_") := last(local_hour) * 60 + last(local_minute)))

            features <- merge(features, feature, by="local_date", all = TRUE)
        }
    }
    features <- features %>% mutate_at(vars(contains("countmostfrequentcontact")), list( ~ replace_na(., 0)))
    return(features)
}