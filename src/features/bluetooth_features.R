source("packrat/init.R")

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

data <- read.csv(snakemake@input[[1]], stringsAsFactors = FALSE)
day_segment <- snakemake@params[["day_segment"]]
requested_features <-  snakemake@params[["features"]]
features = data.frame(local_date = character(), stringsAsFactors = FALSE)

for(requested_feature in requested_features){
  feature <- compute_bluetooth_feature(data, requested_feature, day_segment)
  features <- merge(features, feature, by="local_date", all = TRUE)
}

features <- features %>% mutate_at(vars(contains("countscansmostuniquedevice")), list( ~ replace_na(., 0)))

write.csv(features, snakemake@output[[1]], row.names = FALSE)