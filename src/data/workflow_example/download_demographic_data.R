source("renv/activate.R")
library(RMySQL)
library("dplyr", warn.conflicts = F)
library(readr)
library(stringr)
library(yaml)


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

write_csv(sensor_data, sensor_file)
