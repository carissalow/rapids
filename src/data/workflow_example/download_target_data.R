source("renv/activate.R")
library("dplyr", warn.conflicts = F)
library(readr)
library(stringr)
library(yaml)
library(lubridate)


participant_file <- snakemake@input[["participant_file"]]
sensor_file <- snakemake@output[[1]]

participant <- read_yaml(participant_file)
record_id <- participant$PHONE$LABEL

target_data <- read.csv(snakemake@input[["data"]])
target_data <- target_data[target_data$record_id == record_id, ]

target_data$local_date_time <-  paste(target_data$local_date, "00:00:00")
#target_data <- target_data %>% rename(local_date_time = local_date)

target_data$timestamp <- 0

write_csv(target_data, sensor_file)
