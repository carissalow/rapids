source("renv/activate.R")
source("src/data/unify_utils.R")
library(RMySQL)
library(stringr)
library(dplyr)
library(readr)
library(yaml)
library(lubridate)
options(scipen=999)

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

get_timestamp_filter <- function(device_ids, participant, timezone){
    # Read start and end date from the participant file to filter data within that range
    start_date <- ymd_hms(paste(participant$PHONE$START_DATE,"00:00:00"), tz=timezone, quiet=TRUE)
    end_date <- ymd_hms(paste(participant$PHONE$END_DATE, "23:59:59"), tz=timezone, quiet=TRUE)
    start_timestamp = as.numeric(start_date) * 1000
    end_timestamp = as.numeric(end_date) * 1000
    if(is.na(start_timestamp)){
      message(paste("PHONE[START_DATE] was not provided or failed to parse (", participant$PHONE$START_DATE,"), all data for", paste0(device_ids, collapse=","),"is returned"))
      return("")
    }else if(is.na(end_timestamp)){
      message(paste("PHONE[END_DATE] was not provided or failed to parse (", participant$PHONE$END_DATE,"), all data for", paste0(device_ids, collapse=","),"is returned"))
      return("")
    } else if(start_timestamp > end_timestamp){
      stop(paste("Start date has to be before end date in PHONE[TIME_SPAN] (",start_date,",", date(end_date),"), all data for", paste0(device_ids, collapse=","),"is returned"))
      return("")
    } else {
      message(paste("Filtering data between", start_date, "and", end_date, "in", timezone, "for",paste0(device_ids, collapse=",")))
      return(paste0("AND timestamp BETWEEN ", start_timestamp, " AND ", end_timestamp))
    }
}

participant_file <- snakemake@input[[1]]
source <- snakemake@params[["source"]]
group <- source$DATABASE_GROUP
table <- snakemake@params[["table"]]
sensor <- snakemake@params[["sensor"]]
timezone <- snakemake@params[["timezone"]]
aware_multiplatform_tables <- str_split(snakemake@params[["aware_multiplatform_tables"]], ",")[[1]]
sensor_file <- snakemake@output[[1]]

participant <- read_yaml(participant_file)
if(! "PHONE" %in% names(participant)){
  stop(paste("The following participant file does not have a PHONE section, create one manually or automatically (see the docs):", participant_file))
}
device_ids <- participant$PHONE$DEVICE_IDS
unified_device_id <- tail(device_ids, 1)
platforms <- participant$PHONE$PLATFORMS
validate_deviceid_platforms(device_ids, platforms)
timestamp_filter <- get_timestamp_filter(device_ids, participant, timezone)

dbEngine <- dbConnect(MySQL(), default.file = "./.env", group = group)

if(is_multiplaform_participant(dbEngine, device_ids, platforms)){
  sensor_data <- unify_raw_data(dbEngine, table, sensor, timestamp_filter, aware_multiplatform_tables, device_ids, platforms)
}else {
  # table has two elements for conversation and activity recognition (they store data on a different table for ios and android)
  if(length(table) > 1)
    table <- table[[toupper(platforms[1])]]
  query <- paste0("SELECT * FROM ", table, " WHERE ",source$DEVICE_ID_COLUMN," IN ('", paste0(device_ids, collapse = "','"), "')", timestamp_filter)
  sensor_data <- dbGetQuery(dbEngine, query) %>%
    rename(device_id = source$DEVICE_ID_COLUMN)
}

sensor_data <- sensor_data %>% arrange(timestamp)

# Unify device_id
sensor_data <- sensor_data %>% mutate(device_id = unified_device_id)

# Droping duplicates on all columns except for _id or id
sensor_data <- sensor_data %>% distinct(!!!syms(setdiff(names(sensor_data), c("_id", "id"))))

write_csv(sensor_data, sensor_file)
dbDisconnect(dbEngine)