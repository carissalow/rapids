source("renv/activate.R")
library("dplyr", warn.conflicts = F)
library(readr)
library(stringr)
library(yaml)


participant_file <- snakemake@input[["participant_file"]]

sensor_file <- snakemake@output[[1]]

participant <- read_yaml(participant_file)
record_id <- participant$PHONE$LABEL

demographic_data = read.csv(snakemake@input[["data"]])
demographic_data = demographic_data[demographic_data$record_id == record_id, ]

write_csv(demographic_data, sensor_file)
