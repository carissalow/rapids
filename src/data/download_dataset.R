source("packrat/init.R")

library(RMySQL)
library(stringr)
library(dplyr)

participant <- snakemake@input[[1]]
group <- snakemake@params[["group"]]
table <- snakemake@params[["table"]]
sensor_file <- snakemake@output[[1]]

device_ids <- readLines(participant, n=1)
unified_device_id <- tail(strsplit(device_ids, ",")[[1]], 1)
rmysql.settingsfile <- "./.env"

stopDB <- dbConnect(MySQL(), default.file = rmysql.settingsfile, group = group)
query <- paste0("SELECT * FROM ", table, " WHERE device_id IN ('", gsub(",", "','", device_ids), "')")
sensor_data <- dbGetQuery(stopDB, query)
sensor_data <- sensor_data %>% 
    arrange(timestamp) %>% 
    mutate(device_id = unified_device_id)

# Droping duplicates on all columns except for _id
sensor_data <- sensor_data %>% distinct(!!!syms(setdiff(names(sensor_data), "_id")))
write.csv(sensor_data, sensor_file, row.names = FALSE)
dbDisconnect(stopDB)