source("packrat/init.R")

library(RMySQL)
library(stringr)
library(dplyr)

participant <- snakemake@input[[1]]
group <- snakemake@params[["group"]]
table <- snakemake@params[["table"]]
timezone <- snakemake@params[["timezone"]]
sensor_file <- snakemake@output[[1]]

device_ids <- readLines(participant, n=1)
unified_device_id <- tail(strsplit(device_ids, ",")[[1]], 1)

start_date <- strsplit(readLines(participant, n=4)[4], ",")[[1]][1]
end_date <- strsplit(readLines(participant, n=4)[4], ",")[[1]][2]
start_datetime_utc = format(as.POSIXct(paste0(start_date, " 00:00:00"),format="%Y/%m/%d %H:%M:%S",origin="1970-01-01",tz=timezone), tz="UTC")
end_datetime_utc = format(as.POSIXct(paste0(end_date, " 23:59:59"),format="%Y/%m/%d %H:%M:%S",origin="1970-01-01",tz=timezone), tz="UTC")

rmysql.settingsfile <- "./.env"

stopDB <- dbConnect(MySQL(), default.file = rmysql.settingsfile, group = group)
query <- paste0("SELECT * FROM ", table, " WHERE device_id IN ('", gsub(",", "','", device_ids), "')")
if(!(is.na(start_datetime_utc)) && !(is.na(end_datetime_utc)) && start_datetime_utc < end_datetime_utc){
    query <- paste0(query, "AND timestamp BETWEEN 1000*UNIX_TIMESTAMP('", start_datetime_utc, "') AND 1000*UNIX_TIMESTAMP('", end_datetime_utc, "')")
}
sensor_data <- dbGetQuery(stopDB, query)
sensor_data <- sensor_data %>% 
    arrange(timestamp) %>% 
    mutate(device_id = unified_device_id)

# Droping duplicates on all columns except for _id
sensor_data <- sensor_data %>% distinct(!!!syms(setdiff(names(sensor_data), "_id")))
write.csv(sensor_data, sensor_file, row.names = FALSE)
dbDisconnect(stopDB)
