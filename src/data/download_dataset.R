source("packrat/init.R")

library(RMySQL)
library(stringr)
library(dplyr)

participant <- snakemake@input[[1]]
group <- snakemake@params[["group"]]
table <- snakemake@params[["table"]]
sensor_file <- snakemake@output[[1]]

device_id <- readLines(participant, n=1)
rmysql.settingsfile <- "./.env"

stopDB <- dbConnect(MySQL(), default.file = rmysql.settingsfile, group = group)
query <- paste("SELECT * FROM ", table, " WHERE device_id LIKE '", device_id, "'", sep = "")
sensor_data <- dbGetQuery(stopDB, query)
sensor_data <- sensor_data[order(sensor_data$timestamp),]

# Droping duplicates on all columns except for _id
sensor_data <- sensor_data %>% distinct(!!!syms(setdiff(names(sensor_data), "_id")))
write.table(sensor_data, sensor_file, row.names = FALSE, quote = FALSE, sep = ",")
dbDisconnect(stopDB)