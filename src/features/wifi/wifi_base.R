library(dplyr)

filter_by_day_segment <- function(data, day_segment) {
  if(day_segment %in% c("morning", "afternoon", "evening", "night"))
    data <- data %>% filter(local_day_segment == day_segment)

  return(data %>% group_by(local_date))
}

compute_wifi_feature <- function(data, feature, day_segment){
  if(feature %in% c("countscans", "uniquedevices")){
    data <- data %>% filter_by_day_segment(day_segment)
    data <- switch(feature,
              "countscans" = data %>% summarise(!!paste("wifi", day_segment, feature, sep = "_") := n()),
              "uniquedevices" = data %>% summarise(!!paste("wifi", day_segment, feature, sep = "_") := n_distinct(bssid)))
    return(data)
   } else if(feature == "countscansmostuniquedevice"){
     # Get the most scanned device
    data <- data %>% group_by(bssid) %>% 
      mutate(N=n()) %>% 
      ungroup() %>%
      filter(N == max(N))
    return(data %>% 
             filter_by_day_segment(day_segment) %>%
             summarise(!!paste("wifi", day_segment, feature, sep = "_") := n()))
  }
}

base_wifi_features <- function(wifi_data, day_segment, requested_features){
    # Output dataframe
    features = data.frame(local_date = character(), stringsAsFactors = FALSE)

    # The name of the features this function can compute
    base_features_names  <- c("countscans", "uniquedevices", "countscansmostuniquedevice")

    # The subset of requested features this function can compute
    features_to_compute  <- intersect(base_features_names, requested_features)

    for(feature_name in features_to_compute){
    feature <- compute_wifi_feature(wifi_data, feature_name, day_segment)
    features <- merge(features, feature, by="local_date", all = TRUE)
    }

    return(features)
}
