source("packrat/init.R")

library(tidyr)
library(purrr)
library(dplyr)

metric_files  <- snakemake@input[["metric_files"]]
phone_valid_sensed_days  <- snakemake@input[["phone_valid_sensed_days"]]
days_to_include <- snakemake@input[["days_to_include"]]
source <- snakemake@params[["source"]]

metrics_for_individual_model <- metric_files %>%
  map(read.csv, stringsAsFactors = F, colClasses = c(local_date = "character")) %>%
  reduce(full_join, by="local_date")

if(!is.null(phone_valid_sensed_days) && source %in% c("phone_metrics", "phone_fitbit_metrics")){
    metrics_for_individual_model <- merge(metrics_for_individual_model, read.csv(phone_valid_sensed_days), by="local_date") %>% select(-valid_hours)
}

if(!is.null(days_to_include)){
  metrics_for_individual_model <- merge(metrics_for_individual_model, read.csv(days_to_include), by="local_date")
}

write.csv(metrics_for_individual_model, snakemake@output[[1]], row.names = FALSE)