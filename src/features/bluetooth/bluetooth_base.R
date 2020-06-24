library(dplyr)
library(tidyr)

filter_by_day_segment <- function(data, day_segment) {
  if(day_segment %in% c("morning", "afternoon", "evening", "night"))
    data <- data %>% filter(local_day_segment == day_segment)

  return(data %>% group_by(local_date))
}

compute_bluetooth_feature <- function(data, feature, day_segment){
  if(feature %in% c("countscans", "uniquedevices")){
    data <- data %>% filter_by_day_segment(day_segment)
    data <- switch(feature,
              "countscans" = data %>% summarise(!!paste("bluetooth", day_segment, feature, sep = "_") := n()),
              "uniquedevices" = data %>% summarise(!!paste("bluetooth", day_segment, feature, sep = "_") := n_distinct(bt_address)))
    return(data)
   } else if(feature == "countscansmostuniquedevice"){
     # Get the most scanned device
    data <- data %>% group_by(bt_address) %>% 
      mutate(N=n()) %>% 
      ungroup() %>%
      filter(N == max(N))
    return(data %>% 
             filter_by_day_segment(day_segment) %>%
             summarise(!!paste("bluetooth", day_segment, feature, sep = "_") := n()))
  }
}

base_bluetooth_features <- function(bluetooth_data, day_segment, requested_features){
    # Output dataframe
    features = data.frame(local_date = character(), stringsAsFactors = FALSE)

    # The name of the features this function can compute
    base_features_names  <- c("countscans", "uniquedevices", "countscansmostuniquedevice")

    # The subset of requested features this function can compute
    features_to_compute  <- intersect(base_features_names, requested_features)

    for(feature_name in features_to_compute){
    feature <- compute_bluetooth_feature(bluetooth_data, feature_name, day_segment)
    features <- merge(features, feature, by="local_date", all = TRUE)
    }

    features <- features %>% mutate_at(vars(contains("countscansmostuniquedevice")), list( ~ replace_na(., 0)))

    return(features)
}