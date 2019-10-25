source("packrat/init.R")

library(dplyr)
library(entropy)
library(robustbase)

calls <-  read.csv(snakemake@input[[1]])
day_segment <- snakemake@params[["day_segment"]]
metric <-  snakemake@params[["metric"]]
type <- snakemake@params[["call_type"]]
output_file <- snakemake@output[[1]]

metrics <- calls %>% filter(call_type == ifelse(type == "incoming", "1", ifelse(type == "outgoing", "2", "3")))

if(day_segment == "daily"){
  metrics <- metrics %>% group_by(local_date)
} else {
  metrics <- metrics %>% filter(day_segment == local_day_segment) %>% group_by(local_date)
}

metrics <- switch(metric,
       "count" = metrics %>% summarise(!!paste("com", "call", type, day_segment, metric, sep = "_") := n()),
       "distinctcontacts" = metrics %>% summarise(!!paste("com", "call", type, day_segment, metric, sep = "_") := n_distinct(trace)),
       "meanduration" = metrics %>% summarise(!!paste("com", "call", type, day_segment, metric, sep = "_") := mean(call_duration)),
       "sumduration" = metrics %>% summarise(!!paste("com", "call", type, day_segment, metric, sep = "_") := sum(call_duration)),
       "hubermduration" = metrics %>% summarise(!!paste("com", "call", type, day_segment, metric, sep = "_") := huberM(call_duration)$mu),
       "varqnduration" = metrics %>% summarise(!!paste("com", "call", type, day_segment, metric, sep = "_") := Qn(call_duration)),
       "entropyduration" = metrics %>% summarise(!!paste("com", "call", type, day_segment, metric, sep = "_") := entropy.MillerMadow(call_duration)))

write.csv(na.omit(metrics), output_file, row.names = F)
