source("packrat/init.R")

library(dplyr)

filter_by_day_segment <- function(data, day_segment) {
  if(day_segment %in% c("morning", "afternoon", "evening", "night"))
    data <- data %>% filter(local_day_segment == day_segment)

  return(data %>% group_by(local_date))
}

compute_sms_feature <- function(sms, metric, day_segment){
  sms <- sms %>% filter_by_day_segment(day_segment)
  feature <- switch(metric,
       "count" = sms %>% summarise(!!paste("com", "sms", sms_type, day_segment, metric, sep = "_") := n()),
       "distinctcontacts" = sms %>% summarise(!!paste("com", "sms", sms_type, day_segment, metric, sep = "_") := n_distinct(trace)))
  return(feature)
}

sms <-  read.csv(snakemake@input[[1]])
day_segment <- snakemake@params[["day_segment"]]
metrics <-  snakemake@params[["metrics"]]
sms_type <-  snakemake@params[["sms_type"]]
features = data.frame(local_date = character(), stringsAsFactors = FALSE)

sms <- sms %>% filter(message_type == ifelse(sms_type == "received", "1", "2"))

for(metric in metrics){
  feature <- compute_sms_feature(sms, metric, day_segment)
  features <- merge(features, feature, by="local_date", all = TRUE)
}

write.csv(features, snakemake@output[[1]], row.names = FALSE)
