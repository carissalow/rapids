source("renv/activate.R")
source("src/data/unify_utils.R")
library(RMySQL)
library(stringr)
library(dplyr)
library(readr)

validate_deviceid_platforms <- function(device_ids, platforms){
  if(length(device_ids) == 1){
    if(length(platforms) > 1 || (platforms != "android" && platforms != "ios"))
      stop(paste0("If you have 1 device_id, its platform should be 'android' or 'ios' but you typed: '", paste0(platforms, collapse = ","), "'. Participant file: ", participant))
  } else if(length(device_ids) > 1 && length(platforms) == 1){
    if(platforms != "android" && platforms != "ios" && platforms != "multiple")
      stop(paste0("If you have more than 1 device_id, platform should be 'android', 'ios' OR 'multiple' but you typed: '", paste0(platforms, collapse = "s,"), "'. Participant file: ", participant))
  } else if(length(device_ids) > 1 && length(platforms) > 1){
    if(length(device_ids) != length(platforms))
      stop(paste0("The number of device_ids should match the number of platforms. Participant file:", participant))
    if(all(intersect(c("android", "ios"), unique(platforms)) != c("android", "ios")))
      stop(paste0("If you have more than 1 device_id and more than 1 platform, the platforms should be a mix of 'android' AND 'ios' but you typed: '", paste0(platforms, collapse = ","), "'. Participant file: ", participant))
  }
}

is_multiplaform_participant <- function(dbEngine, device_ids, platforms){
  # Multiple android and ios platforms or the same platform (android, ios) for multiple devices
  if((length(device_ids) > 1 && length(platforms) > 1) || (length(device_ids) > 1 && length(platforms) == 1 && (platforms == "android" || platforms == "ios"))){
    return(TRUE)
  }
  # Multiple platforms for multiple devices, we search the platform for every device in the aware_device table
  if(length(device_ids) > 1 && length(platforms) == 1 && platforms == "multiple"){
    devices_platforms <- dbGetQuery(dbEngine, paste0("SELECT device_id,brand FROM aware_device WHERE device_id IN ('", paste0(device_ids, collapse = "','"), "')"))
    platforms <- devices_platforms %>% distinct(brand) %>% pull(brand)
    # Android phones have different brands so we check that we got at least two different platforms and one of them is iPhone
    if(length(platforms) > 1 && "iPhone" %in% platforms){
      return(TRUE)
    }
  }
  return(FALSE)
}

participant <- snakemake@input[[1]]
group <- snakemake@params[["group"]]
table <- snakemake@params[["table"]]
timezone <- snakemake@params[["timezone"]]
aware_multiplatform_tables <- str_split(snakemake@params[["aware_multiplatform_tables"]], ",")[[1]]
unifiable_tables = snakemake@params[["unifiable_sensors"]]
sensor_file <- snakemake@output[[1]]

device_ids <- strsplit(readLines(participant, n=1), ",")[[1]]
unified_device_id <- tail(device_ids, 1)
platforms <- strsplit(readLines(participant, n=2)[[2]], ",")[[1]]
validate_deviceid_platforms(device_ids, platforms)

# Read start and end date from the participant file to filter data within that range
start_date <- strsplit(readLines(participant, n=4)[4], ",")[[1]][1]
end_date <- strsplit(readLines(participant, n=4)[4], ",")[[1]][2]
start_datetime_utc = format(as.POSIXct(paste0(start_date, " 00:00:00"),format="%Y/%m/%d %H:%M:%S",origin="1970-01-01",tz=timezone), tz="UTC")
end_datetime_utc = format(as.POSIXct(paste0(end_date, " 23:59:59"),format="%Y/%m/%d %H:%M:%S",origin="1970-01-01",tz=timezone), tz="UTC")

dbEngine <- dbConnect(MySQL(), default.file = "./.env", group = group)

# Get existent columns in table
available_columns <- colnames(dbGetQuery(dbEngine, paste0("SELECT * FROM ", table, " LIMIT 1")))

if("device_id" %in% available_columns){
  if(is_multiplaform_participant(dbEngine, device_ids, platforms)){
    sensor_data <- unify_raw_data(dbEngine, table, start_datetime_utc, end_datetime_utc, aware_multiplatform_tables, unifiable_tables, device_ids, platforms)
  }else {
    query <- paste0("SELECT * FROM ", table, " WHERE device_id IN ('", paste0(device_ids, collapse = "','"), "')")
    if("timestamp" %in% available_columns && !(is.na(start_datetime_utc)) && !(is.na(end_datetime_utc)) && start_datetime_utc < end_datetime_utc)
      query <- paste0(query, "AND timestamp BETWEEN 1000*UNIX_TIMESTAMP('", start_datetime_utc, "') AND 1000*UNIX_TIMESTAMP('", end_datetime_utc, "')")
    sensor_data <- dbGetQuery(dbEngine, query)
  }
  
  if("timestamp" %in% available_columns)
    sensor_data <- sensor_data %>% arrange(timestamp)
  
  # Unify device_id
  sensor_data <- sensor_data %>% mutate(device_id = unified_device_id)
  
  # Droping duplicates on all columns except for _id or id
  sensor_data <- sensor_data %>% distinct(!!!syms(setdiff(names(sensor_data), c("_id", "id"))))
  
} else 
    stop(paste0("Table ", table, "does not have a device_id column (Aware ID) to link its data to a participant"))

write_csv(sensor_data, sensor_file)
dbDisconnect(dbEngine)