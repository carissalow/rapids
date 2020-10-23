source("renv/activate.R")

library(tidyr)
library(purrr)
library("dplyr", warn.conflicts = F)
library("methods")
library("mgm")
library("qgraph")
library("dplyr", warn.conflicts = F)
library("scales")
library("ggplot2")
library("purrr")
library("tidyr")
library("reshape2")

feature_files  <- snakemake@input[["feature_files"]]
phone_valid_sensed_days  <- snakemake@input[["phone_valid_sensed_days"]]
days_to_include <- snakemake@input[["days_to_include"]]
source <- snakemake@params[["source"]]

features_for_individual_model <- feature_files %>%
  map(read.csv, stringsAsFactors = F, colClasses = c(local_date = "character")) %>%
  reduce(full_join, by="local_date")

if(!is.null(phone_valid_sensed_days) && source %in% c("phone_features", "phone_fitbit_features")){
    valid_days <- read.csv(phone_valid_sensed_days)
    valid_days <- valid_days[valid_days$is_valid_sensed_day == TRUE, ]
    features_for_individual_model <- merge(features_for_individual_model, valid_days, by="local_date") %>% select(-valid_sensed_hours, -is_valid_sensed_day)
}

if(!is.null(days_to_include)){
  features_for_individual_model <- merge(features_for_individual_model, read.csv(days_to_include), by="local_date")
}

write.csv(features_for_individual_model, snakemake@output[[1]], row.names = FALSE)