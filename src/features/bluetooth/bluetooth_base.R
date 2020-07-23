library(dplyr)
library(tidyr)

filter_by_day_segment <- function(data, day_segment) {
  if(day_segment != "daily")
    data <- data %>% filter(local_day_segment == day_segment)

  return(data %>% group_by(local_date))
}

compute_bluetooth_feature <- function(data, feature, day_segment){
  data <- data %>% filter_by_day_segment(day_segment)
  if(feature %in% c("countscans", "uniquedevices")){
    data <- switch(feature,
              "countscans" = data %>% summarise(!!paste("bluetooth", day_segment, feature, sep = "_") := n()),
              "uniquedevices" = data %>% summarise(!!paste("bluetooth", day_segment, feature, sep = "_") := n_distinct(bt_address)))
    return(data)
   } else if(feature == "countscansmostuniquedevice"){
     # Get the most scanned device
    mostuniquedevice <- data %>% 
      group_by(bt_address) %>% 
      mutate(N=n()) %>% 
      ungroup() %>%
      filter(N == max(N)) %>% 
      head(1) %>% # if there are multiple device with the same amount of scans pick the first one only
      pull(bt_address)
    return(data %>% 
             filter(bt_address == mostuniquedevice) %>%
             group_by(local_date) %>% 
             summarise(!!paste("bluetooth", day_segment, feature, sep = "_") := n()) %>%
             replace(is.na(.), 0))
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