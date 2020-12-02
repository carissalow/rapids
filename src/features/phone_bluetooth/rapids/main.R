library("dplyr", warn.conflicts = F)
library(tidyr)

compute_bluetooth_feature <- function(data, feature, time_segment){
  data <- data %>% filter_data_by_segment(time_segment)
  if(feature %in% c("countscans", "uniquedevices")){
    data <- data %>% group_by(local_segment)
    data <- switch(feature,
              "countscans" = data %>% summarise(!!feature := n()),
              "uniquedevices" = data %>% summarise(!!feature := n_distinct(bt_address)))
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
    mostuniquedevice
    return(data %>% 
             filter(bt_address == mostuniquedevice) %>%
             group_by(local_segment) %>% 
             summarise(!!feature := n()) %>%
             replace(is.na(.), 0))
  }
}

rapids_features <- function(sensor_data_files, time_segment, provider){
  
  bluetooth_data <-  read.csv(sensor_data_files[["sensor_data"]], stringsAsFactors = FALSE)
  requested_features <- provider[["FEATURES"]]
  
  # Output dataframe
  features = data.frame(local_segment = character(), stringsAsFactors = FALSE)

  # The name of the features this function can compute
  base_features_names  <- c("countscans", "uniquedevices", "countscansmostuniquedevice")

  # The subset of requested features this function can compute
  features_to_compute  <- intersect(base_features_names, requested_features)

  for(feature_name in features_to_compute){
    feature <- compute_bluetooth_feature(bluetooth_data, feature_name, time_segment)
    features <- merge(features, feature, by="local_segment", all = TRUE)
  }

  features <- features %>% mutate_at(vars(contains("countscansmostuniquedevice")), list( ~ replace_na(., 0)))

  return(features)
}