source("packrat/init.R")

library(dplyr)

sms <-  read.csv(snakemake@input[[1]])
day_segment <- snakemake@params[["day_segment"]]
metric <-  snakemake@params[["metric"]]
sms_type <-  snakemake@params[["sms_type"]]
output_file <- snakemake@output[[1]]

metrics <- sms %>% filter(message_type == ifelse(sms_type == "received", "1", "2"))

if(day_segment == "daily"){
  metrics <- metrics %>% group_by(local_date)
} else {
  metrics <- metrics %>% filter(day_segment == local_day_segment) %>% group_by(local_date)
}

metrics <- switch(metric,
       "count" = metrics %>% summarise(!!paste("com", "sms", sms_type, day_segment, metric, sep = "_") := n()),
       "distinctcontacts" = metrics %>% summarise(!!paste("com", "sms", sms_type, day_segment, metric, sep = "_") := n_distinct(trace)))

write.csv(na.omit(metrics), output_file, row.names = F)