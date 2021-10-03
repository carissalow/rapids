library('tidyr')
library('stringr')
library('entropy')

Mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

call_features_of_type <- function(calls, features_type, call_type, time_segment, requested_features){
    # Output dataframe
    features = data.frame(local_segment = character(), stringsAsFactors = FALSE)

    # The name of the features this function can compute
    base_features_names  <- c("count", "distinctcontacts", "meanduration", "sumduration", "minduration", "maxduration", "stdduration", "modeduration", "entropyduration", "timefirstcall", "timelastcall", "countmostfrequentcontact")
    # The subset of requested features this function can compute
    features_to_compute  <- intersect(base_features_names, requested_features)

    # If there are not features or data to work with, return an empty df with appropiate columns names
    if(length(features_to_compute) == 0)
        return(features)
    if(nrow(calls) < 1)
        return(cbind(features, read.csv(text = paste(paste(call_type, features_to_compute, sep = "_"), collapse = ","), stringsAsFactors = FALSE)))

    if(features_type == "EPISODES"){
        calls <- calls %>% 
            mutate(call_duration = (end_timestamp - start_timestamp) / 1000) %>% 
            separate(local_start_date_time, c("local_date","local_time"), "\\s", remove = FALSE) %>%
            separate(local_time, c("local_hour", "local_minute"), ":", remove = FALSE, extra = "drop") %>%
            mutate(local_hour = as.numeric(local_hour),
                local_minute = as.numeric(local_minute))
    }

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
                group_by(local_segment) %>% 
                summarise(!!paste(call_type, feature_name, sep = "_") := sum(trace == mostfrequentcontact))
            features <- merge(features, feature, by="local_segment", all = TRUE)
        } else {
            feature <- calls %>% 
                group_by(local_segment)

            feature <- switch(feature_name,
                "count" = feature %>% summarise(!!paste(call_type, feature_name, sep = "_") := n()),
                "distinctcontacts" = feature %>% summarise(!!paste(call_type, feature_name, sep = "_") := n_distinct(trace)),
                "meanduration" = feature %>% summarise(!!paste(call_type, feature_name, sep = "_") := mean(call_duration)),
                "sumduration" = feature %>% summarise(!!paste(call_type, feature_name, sep = "_") := sum(call_duration)),
                "minduration" = feature %>% summarise(!!paste(call_type, feature_name, sep = "_") := min(call_duration)),
                "maxduration" = feature %>% summarise(!!paste(call_type, feature_name, sep = "_") := max(call_duration)),
                "stdduration" = feature %>% summarise(!!paste(call_type, feature_name, sep = "_") := sd(call_duration)),
                "modeduration" = feature %>% summarise(!!paste(call_type, feature_name, sep = "_") := Mode(call_duration)),
                "entropyduration" = feature %>% summarise(!!paste(call_type, feature_name, sep = "_") := entropy.MillerMadow(call_duration)),
                "timefirstcall" = feature %>% summarise(!!paste(call_type, feature_name, sep = "_") := first(local_hour) * 60 + first(local_minute)),
                "timelastcall" = feature %>% summarise(!!paste(call_type, feature_name, sep = "_") := last(local_hour) * 60 + last(local_minute)))

            features <- merge(features, feature, by="local_segment", all = TRUE)
        }
    }
    return(features)
}

rapids_features <- function(sensor_data_files, time_segment, provider){
    calls_data <-  read.csv(sensor_data_files[["sensor_data"]], stringsAsFactors = FALSE)
    calls_data <- calls_data %>% filter_data_by_segment(time_segment)

    features_type <- provider[["FEATURES_TYPE"]]
    call_types = provider[["CALL_TYPES"]]
    call_features <- setNames(data.frame(matrix(ncol=1, nrow=0)), c("local_segment"))

    for(call_type in call_types){
        # Filter rows that belong to the calls type and time segment of interest
        call_type_label = ifelse(call_type == "incoming", "1", ifelse(call_type == "outgoing", "2", ifelse(call_type == "missed", "3", NA)))
        if(is.na(call_type_label))
            stop(paste("Call type can online be incoming, outgoing or missed but instead you typed: ", call_type, " in config[CALLS][CALL_TYPES]"))

        requested_features <- provider[["FEATURES"]][[call_type]]
        calls_of_type <- calls_data %>% filter(call_type == call_type_label)

        features <- call_features_of_type(calls_of_type, features_type, call_type, time_segment, requested_features)
        call_features <- merge(call_features, features, all=TRUE)
    }
    call_features <- call_features %>% mutate_at(vars(contains("countmostfrequentcontact") | contains("distinctcontacts") | contains("count") | contains("sumduration") | contains("minduration") | contains("maxduration") | contains("meanduration") | contains("modeduration")), list( ~ replace_na(., 0)))
    return(call_features)
}