
library(tibble)
library(dplyr)
library(tidyr)
library(purrr)
library(yaml)
options(scipen = 999)

buils_tz_intervals <- function(tz_codes){
  tz_codes <- tz_codes %>%
    group_by(device_id) %>% 
    mutate(end_timestamp = lead(timestamp)) %>% 
    ungroup() %>% 
    replace_na(list(end_timestamp = as.numeric(Sys.time())*1000))
  return(tz_codes)
}

filter_tz_per_device <- function(device_id, tz_codes, default, IF_MISSING_TZCODE){
  device_tz_codes <- tz_codes %>% filter(device_id == !!device_id) %>% select(-device_id)

  if(nrow(device_tz_codes) > 0)
    return(device_tz_codes)
  else if(IF_MISSING_TZCODE == "STOP")
    stop(paste("The device id '", device_id, "' does not have any time zone codes in your [MULTIPLE][TZCODES_FILE], add one or set IF_MISSING_TZCODE to 'USE_DEFAULT'"))
  else if(IF_MISSING_TZCODE == "USE_DEFAULT")
    return(data.frame(timestamp = c(0), tzcode = default, end_timestamp = as.numeric(Sys.time())*1000))
  
  stop("We should have obtained the time zones for a device, stop the execution or use the default tz but this didn't happen. Create an issue on Github")
}

assign_tz_code <- function(data, tz_codes){

  for(i in 1:nrow(tz_codes)) {
    start_timestamp <- tz_codes[[i, "timestamp"]]
    end_timestamp <- tz_codes[[i, "end_timestamp"]]
    time_zone <- trimws(tz_codes[[i, "tzcode"]], which="both")
    
    data$local_timezone <- if_else(start_timestamp <= data$timestamp & data$timestamp < end_timestamp, time_zone, data$local_timezone)
  }
  return(data %>% filter(!is.na(local_timezone)))
  
}

validate_single_tz_per_fitbit_device <- function(tz_codes, INFER_FROM_SMARTPHONE_TZ){
  
  if(INFER_FROM_SMARTPHONE_TZ)
    stop("If [TIMEZONE][MULTIPLE][FITBIT][INFER_FROM_SMARTPHONE_TZ] is True (you want to infer Fitbit time zones with smartphone data), you need to set ALLOW_MULTIPLE_TZ_PER_DEVICE to True. However, read the docs to understand why this can be innacurate")
  
  tz_per_device <- tz_codes %>% group_by(device_id) %>% summarise(n = n(), .groups = "drop_last") %>% filter(n > 1)
  
  if(nrow(tz_per_device) > 0)
    stop(paste("The following Fitbit device ids have more than one time zone change which is not allowed if [TIMEZONE][MULTIPLE][FITBIT][ALLOW_MULTIPLE_TZ_PER_DEVICE] is False:", paste(tz_per_device %>% pull(device_id), collapse = ",")))
  
  zero_ts <- tz_codes %>% filter(timestamp > 0)

  if(nrow(zero_ts) > 0)
    stop(paste("The following Fitbit device ids have a time zone change with a timestamp bigger than 0 which is not allowed if [TIMEZONE][MULTIPLE][FITBIT][ALLOW_MULTIPLE_TZ_PER_DEVICE] is False: ", paste(zero_ts %>% pull(device_id), collapse = ",")))
  
}

validate_devies_exist_in_participant_file <- function(devices, device_type, pid, participant_file){
  if(length(devices) == 0)
    stop("[TIMEZONE][MULTIPLE][FITBIT][ALLOW_MULTIPLE_TZ_PER_DEVICE] is True (you want to infer Fitbit time zones with smartphone data), however participant ", pid," does not have any [",device_type,"][DEVICE_IDS] in ", participant_file)
  
}

# TODO include CSV timezone file in rule
multiple_time_zone_assignment <- function(sensor_data, timezone_parameters, device_type, pid, participant_file){
    tz_codes <- read.csv(timezone_parameters$MULTIPLE$TZCODES_FILE)
    default <- timezone_parameters$MULTIPLE$DEFAULT_TZCODE
    IF_MISSING_TZCODE <- timezone_parameters$MULTIPLE$IF_MISSING_TZCODE
    ALLOW_MULTIPLE_TZ_PER_DEVICE <- timezone_parameters$MULTIPLE$FITBIT$ALLOW_MULTIPLE_TZ_PER_DEVICE
    INFER_FROM_SMARTPHONE_TZ <- timezone_parameters$MULTIPLE$FITBIT$INFER_FROM_SMARTPHONE_TZ

    participant_data <- read_yaml(participant_file)
    phone_ids <- participant_data$PHONE$DEVICE_IDS
    fitbit_ids <- participant_data$FITBIT$DEVICE_IDS

    if(device_type == "fitbit"){
        if(!ALLOW_MULTIPLE_TZ_PER_DEVICE){
            validate_single_tz_per_fitbit_device(tz_codes, INFER_FROM_SMARTPHONE_TZ)
        } else if(INFER_FROM_SMARTPHONE_TZ){
            validate_devies_exist_in_participant_file(phone_ids, "PHONE", pid, participant_file)
            validate_devies_exist_in_participant_file(fitbit_ids, "FITBIT", pid, participant_file)
            unified_device_id <- paste0("unified_device_id", pid)
            
            sensor_data <- sensor_data %>% mutate(device_id = if_else(device_id %in% phone_ids, unified_device_id, device_id))
            tz_codes <- tz_codes %>% mutate(device_id = if_else(device_id %in% fitbit_ids, unified_device_id, device_id))
        }
    }
    
    tz_intervals <- buils_tz_intervals(tz_codes)
    sensor_data <- sensor_data %>% mutate(local_timezone = NA_character_)

    if(nrow(sensor_data) > 0){
      sensor_data <- sensor_data %>%
          group_by(device_id) %>% 
          nest() %>% 
          mutate(tz_codes_per_device = map(device_id, filter_tz_per_device, tz_intervals, default, IF_MISSING_TZCODE)) %>% 
          mutate(data = map2(data, tz_codes_per_device, assign_tz_code )) %>% 
          select(-tz_codes_per_device) %>% 
          unnest(cols = data)
    }
    return(sensor_data)
}
