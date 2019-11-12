source("packrat/init.R")

library(dplyr)
library(entropy)
library(robustbase)

filter_by_day_segment <- function(data, day_segment) {
  if(day_segment %in% c("morning", "afternoon", "evening", "night"))
    data <- data %>% filter(local_day_segment == day_segment)

  return(data %>% group_by(local_date))
}

Mode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

compute_call_feature <- function(calls, metric, day_segment){
  if(metric == "countmostfrequentcontact"){
     # Get the most frequent contact
    calls <- calls %>% group_by(trace) %>% 
      mutate(N=n()) %>% 
      ungroup() %>%
      filter(N == max(N))

    return(calls %>% 
             filter_by_day_segment(day_segment) %>%
             summarise(!!paste("call", type, day_segment, metric, sep = "_") := n()))
  } else {
    calls <- calls %>% filter_by_day_segment(day_segment)
    feature <- switch(metric,
        "count" = calls %>% summarise(!!paste("call", type, day_segment, metric, sep = "_") := n()),
        "distinctcontacts" = calls %>% summarise(!!paste("call", type, day_segment, metric, sep = "_") := n_distinct(trace)),
        "meanduration" = calls %>% summarise(!!paste("call", type, day_segment, metric, sep = "_") := mean(call_duration)),
        "sumduration" = calls %>% summarise(!!paste("call", type, day_segment, metric, sep = "_") := sum(call_duration)),
        "minduration" = calls %>% summarise(!!paste("call", type, day_segment, metric, sep = "_") := min(call_duration)),
        "maxduration" = calls %>% summarise(!!paste("call", type, day_segment, metric, sep = "_") := max(call_duration)),
        "stdduration" = calls %>% summarise(!!paste("call", type, day_segment, metric, sep = "_") := sd(call_duration)),
        "modeduration" = calls %>% summarise(!!paste("call", type, day_segment, metric, sep = "_") := Mode(call_duration)),
        "hubermduration" = calls %>% summarise(!!paste("call", type, day_segment, metric, sep = "_") := huberM(call_duration)$mu),
        "varqnduration" = calls %>% summarise(!!paste("call", type, day_segment, metric, sep = "_") := Qn(call_duration)),
        "entropyduration" = calls %>% summarise(!!paste("call", type, day_segment, metric, sep = "_") := entropy.MillerMadow(call_duration)),
        "timefirstcall" = calls %>% summarise(!!paste("call", type, day_segment, metric, sep = "_") := first(local_hour) + (first(local_minute)/60)),
        "timelastcall" = calls %>% summarise(!!paste("call", type, day_segment, metric, sep = "_") := last(local_hour) + (last(local_minute)/60)))
    return(feature)
  }
}

calls <-  read.csv(snakemake@input[[1]], stringsAsFactors = FALSE)
day_segment <- snakemake@params[["day_segment"]]
metrics <-  snakemake@params[["metrics"]]
type <- snakemake@params[["call_type"]]
features = data.frame(local_date = character(), stringsAsFactors = FALSE)

calls <- calls %>% filter(call_type == ifelse(type == "incoming", "1", ifelse(type == "outgoing", "2", "3")))

for(metric in metrics){
  feature <- compute_call_feature(calls, metric, day_segment)
  features <- merge(features, feature, by="local_date", all = TRUE)
}

write.csv(features, snakemake@output[[1]], row.names = FALSE)
