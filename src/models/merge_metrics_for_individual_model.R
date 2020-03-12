source("packrat/init.R")

library(tidyr)
library(purrr)
library(dplyr)

metric_files  <- snakemake@input[["metric_files"]]
phone_valid_sensed_days  <- read.csv(snakemake@input[["phone_valid_sensed_days"]])
drop_valid_sensed_days <- snakemake@params[["drop_valid_sensed_days"]]
source <- snakemake@params[["source"]]

metrics_for_individual_model <- metric_files %>%
  map(read.csv, stringsAsFactors = F, colClasses = c(local_date = "character")) %>%
  reduce(full_join, by="local_date")

if(drop_valid_sensed_days && source == "phone_metrics"){
    metrics_for_individual_model <- merge(metrics_for_individual_model, phone_valid_sensed_days, by="local_date") %>% select(-valid_hours)
  }

write.csv(metrics_for_individual_model, snakemake@output[[1]], row.names = FALSE)