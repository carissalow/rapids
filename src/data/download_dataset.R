if (exists("snakemake"))
    source("packrat/init.R")

library(RMySQL)
library(stringr)
library(dplyr)

participant <- snakemake@input[[1]]
group <- snakemake@params[[1]]
sensor_file <- snakemake@output[[1]]

device_id <- readLines(participant, n=1)
rmysql.settingsfile <- "./.env"
sensor <- tools::file_path_sans_ext(basename(sensor_file))

stopDB <- dbConnect(MySQL(), default.file = rmysql.settingsfile, group = group)
query <- paste("SELECT * FROM ", sensor, " WHERE device_id LIKE '", device_id, "'", sep = "")
sensor_data <- dbGetQuery(stopDB, query)
sensor_data <- sensor_data[order(sensor_data$timestamp),]

# Droping duplicates on all columns except for _id
sensor_data <- sensor_data %>% distinct(!!!syms(setdiff(names(sensor_data), "_id")))
write.table(sensor_data, sensor_file, row.names = FALSE, quote = FALSE, sep = ",")
dbDisconnect(stopDB)