source("packrat/init.R")

library(dplyr)

filter_by_day_segment <- function(data, day_segment) {
  if(day_segment %in% c("morning", "afternoon", "evening", "night"))
    data <- data %>% filter(local_day_segment == day_segment)

  return(data %>% group_by(local_date))
}

compute_bluetooth_metric <- function(data, metric, day_segment){
  if(metric %in% c("countscans", "uniquedevices")){
    data <- data %>% filter_by_day_segment(day_segment)
    data <- switch(metric,
              "countscans" = data %>% summarise(!!paste("bluetooth", day_segment, metric, sep = "_") := n()),
              "uniquedevices" = data %>% summarise(!!paste("bluetooth", day_segment, metric, sep = "_") := n_distinct(bt_address)))
    return(data)
   } else if(metric == "countscansmostuniquedevice"){
     # Get the most scanned device
    data <- data %>% group_by(bt_address) %>% 
      mutate(N=n()) %>% 
      ungroup() %>%
      filter(N == max(N))
    return(data %>% 
             filter_by_day_segment(day_segment) %>%
             summarise(!!paste("bluetooth", day_segment, metric, sep = "_") := n()))
  }
}

data <- read.csv(snakemake@input[[1]], stringsAsFactors = FALSE)
day_segment <- snakemake@params[["day_segment"]]
metrics <-  snakemake@params[["metrics"]]
features = data.frame(local_date = character(), stringsAsFactors = FALSE)

for(metric in metrics){
  feature <- compute_bluetooth_metric(data, metric, day_segment)
  features <- merge(features, feature, by="local_date", all = TRUE)
}

write.csv(features, snakemake@output[[1]], row.names = FALSE)