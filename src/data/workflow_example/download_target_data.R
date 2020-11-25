source("renv/activate.R")
library(RMySQL)
library("dplyr", warn.conflicts = F)
library(readr)
library(stringr)
library(yaml)
library(lubridate)


participant_file <- snakemake@input[["participant_file"]]
source <- snakemake@params[["source"]]
table <- snakemake@params[["table"]]
sensor_file <- snakemake@output[[1]]

participant <- read_yaml(participant_file)
record_id <- participant$PHONE$LABEL

dbEngine <- dbConnect(MySQL(), default.file = "./.env", group = source$DATABASE_GROUP)
query <- paste0("SELECT * FROM ", table, " WHERE record_id = '", record_id, "'")
sensor_data <- dbGetQuery(dbEngine, query)
dbDisconnect(dbEngine)

# generate timestamp based on local_date
sensor_data$timestamp <- as.numeric(ymd_hms(paste(sensor_data$local_date, "00:00:00"), tz=source$TIMEZONE, quiet=TRUE)) * 1000

write_csv(sensor_data, sensor_file)
