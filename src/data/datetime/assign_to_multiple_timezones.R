library(tibble)
library(dplyr)
library(tidyr)
library(purrr)
library(yaml)
library(glue)
library(lubridate)
options(scipen = 999)

buils_tz_intervals <- function(tz_codes, device_type){
  tz_codes <- tz_codes %>%
    group_by(device_id) %>%
    arrange(timestamp)
  
  if(device_type == "fitbit" )
    tz_codes <- tz_codes %>%
      mutate(end_timestamp = lead(timestamp), end_local_date_time = lead(local_date_time)) %>% 
      ungroup() %>% 
      replace_na(list(end_timestamp = as.numeric(Sys.time())*1000, end_local_date_time = format(Sys.time(), format="%Y-%m-%d %H:%M:%S") )) %>% 
      mutate(local_date_time = lubridate::ymd_hms(local_date_time), end_local_date_time = lubridate::ymd_hms(end_local_date_time))
  else
    tz_codes <- tz_codes %>%
      mutate(end_timestamp = lead(timestamp)) %>% 
      ungroup() %>% 
      replace_na(list(end_timestamp = as.numeric(Sys.time())*1000 ))
  return(tz_codes)
}


assign_tz_code <- function(data, device_id, tz_codes, device_type){
  tz_codes <- tz_codes %>% filter(device_id == !!device_id) %>% select(-device_id)
  
  if(device_type == "fitbit" && all(data$timestamp == 0) && "local_date_time" %in% colnames(data)){
    # Only do this for Fitbit raw data, other devices and Fitbit sleep episodes can use timestamps
    column <- "local_date_time"
    data <- data %>% mutate(local_date_time=lubridate::ymd_hms(local_date_time))
  } else
    column <- "timestamp"
  
  for(i in 1:nrow(tz_codes)) {
    start <- tz_codes[[i, column]]
    end <- tz_codes[[i, paste0("end_", column)]]
    time_zone <- trimws(tz_codes[[i, "tzcode"]], which="both")

    data$local_timezone <- if_else(start <= data[[column]] & data[[column]] < end, time_zone, data$local_timezone)
  }
  
  if(column == "local_date_time")
    data$local_date_time <- format(data$local_date_time, format="%Y-%m-%d %H:%M:%S")

  if(any(is.na(data$local_timezone))){
    data_without_timezones <- data %>% filter(is.na(local_timezone))
    warning(glue("We are discarding data for device {device_id} because some rows could not be assigned to a timezone. 
                    This happens because the timezones in TZCODES_FILE are not covering the total period this device was sensing data. 
                    You have two possible fixes
                    - A) Add another timezone entry to TZCODES_FILE that covers the missing period from {min(data_without_timezones[[column]])} to {max(data_without_timezones[[column]])}.
                          Note there will be multiple chunks because we read the sensor data in batches but you could combine the suggested individual missing periods into a single one. 
                    - B) Assign this device id to a single timezone for the entire duration of the study. 
                          To do this, in TZCODES_FILE, replace any timezones for {device_id} with following line:
                              {device_id},{time_zone},0 \n\n
                  "))
  }

  return(data %>% filter(!is.na(local_timezone)))
  
}

validate_single_tz_per_fitbit_device <- function(tz_codes, INFER_FROM_SMARTPHONE_TZ){
  
  if(INFER_FROM_SMARTPHONE_TZ)
    stop("If [TIMEZONE][MULTIPLE][FITBIT][INFER_FROM_SMARTPHONE_TZ] is True (you want to infer Fitbit time zones with smartphone data), you need to set ALLOW_MULTIPLE_TZ_PER_DEVICE to True. However, read the docs to understand why this can be innacurate")
  
  tz_per_device <- tz_codes %>% group_by(device_id) %>% summarise(n = n(), .groups = "drop_last") %>% filter(n > 1)
  if(nrow(tz_per_device) > 0)
    stop("If [TIMEZONE][MULTIPLE][FITBIT][ALLOW_MULTIPLE_TZ_PER_DEVICE] is False, every fitbit device id in [MULTIPLE][TZCODES_FILE] must have only one timezone with a timestamp equal to 0. The following device ids do not comply:", paste(tz_per_device %>% pull(device_id), collapse = ","))
  
  zero_ts <- tz_codes %>% filter(timestamp > 0)
  if(nrow(zero_ts) > 0)
    stop("If [TIMEZONE][MULTIPLE][FITBIT][ALLOW_MULTIPLE_TZ_PER_DEVICE] is False, every fitbit device id in [MULTIPLE][TZCODES_FILE] must have only one timezone with a timestamp equal to 0. The following device ids do not comply:", paste(zero_ts %>% pull(device_id), collapse = ","))
  
}

infer_tz_codes_from_phones <- function(data_device_ids, tz_codes, participant_file){
  participant_data <- read_yaml(participant_file)
  phone_device_ids <- participant_data$PHONE$DEVICE_IDS
  phone_tz_codes <- tz_codes %>% filter(device_id %in% phone_device_ids)
  
  if(nrow(phone_tz_codes) == 0)
    warning("The PHONE device ids that we were supposed to use to infer fitbit devices timezones do not have timezone data in [MULTIPLE][TZCODES_FILE]. ",
            " Problematic phone devices ids: ", paste0(phone_device_ids, collapse = ", "))
  
  data_tz_codes <- NULL
  for(data_device_id in data_device_ids)
    data_tz_codes <- bind_rows(data_tz_codes, phone_tz_codes %>% mutate(device_id = data_device_id) %>% arrange(timestamp))
  data_tz_codes
}

get_devices_ids <- function(participant_data){
  devices_ids = c()
  for(device in participant_data)
    for(attribute in names(device))
      if(attribute == "DEVICE_IDS")
        devices_ids <- c(devices_ids, device[[attribute]])
      return(devices_ids)
}

get_participant_most_common_tz <- function(tz_codes_file, participant_file){
  tz_codes <- read.csv(tz_codes_file)
  participant_device_ids <- get_devices_ids(read_yaml(participant_file))
  
  participant_tz_codes <- tz_codes %>% filter(device_id %in% participant_device_ids)
  most_common_tz <- buils_tz_intervals(participant_tz_codes, "all") %>% 
    mutate(duration = end_timestamp - timestamp) %>% 
    filter(duration == max(duration)) %>% 
    head(1) %>% 
    pull(tzcode)

  if(length(most_common_tz)==0)
    most_common_tz <- "UTC"
  return(most_common_tz)
}

multiple_time_zone_assignment <- function(sensor_data, timezone_parameters, device_type, pid, participant_file){
  if(nrow(sensor_data) == 0)
    return(sensor_data %>% mutate(local_timezone = NA_character_))
  
  tz_codes <- read.csv(timezone_parameters$MULTIPLE$TZCODES_FILE)
  default <- timezone_parameters$MULTIPLE$DEFAULT_TZCODE
  IF_MISSING_TZCODE <- timezone_parameters$MULTIPLE$IF_MISSING_TZCODE
  ALLOW_MULTIPLE_TZ_PER_DEVICE <- timezone_parameters$MULTIPLE$FITBIT$ALLOW_MULTIPLE_TZ_PER_DEVICE
  INFER_FROM_SMARTPHONE_TZ <- timezone_parameters$MULTIPLE$FITBIT$INFER_FROM_SMARTPHONE_TZ
  data_device_ids <- sensor_data %>% distinct(device_id) %>% pull(device_id)
  
  if(INFER_FROM_SMARTPHONE_TZ && device_type == "fitbit")
    data_tz_codes <- infer_tz_codes_from_phones(data_device_ids, tz_codes, participant_file)
  else
    data_tz_codes <- tz_codes %>% filter(device_id %in% data_device_ids)

  # Check if we have timezones for all device ids
  if(length(unique(data_tz_codes$device_id)) < length(data_device_ids)){
    if(IF_MISSING_TZCODE == "STOP")
      stop(glue("One or more device ids do not have any time zone codes in your [MULTIPLE][TZCODES_FILE].",
                "You can add one or set [MULTIPLE][IF_MISSING_TZCODE] to 'USE_DEFAULT'. The missing device ids are [{ids}]",
                ids=paste0(setdiff(data_device_ids, data_tz_codes %>% pull(device_id)), collapse = ",")))
    else if(IF_MISSING_TZCODE == "USE_DEFAULT"){
      warning("Using DEFAULT time zone for ", paste0(setdiff(data_device_ids, data_tz_codes %>% pull(device_id)), collapse = ","))
      default_tz_codes <- data.frame(timestamp = rep_along(c(data_device_ids),0), tzcode = rep_along(c(data_device_ids),default), device_id=data_device_ids) %>% 
        filter(!device_id %in% data_tz_codes$device_id)
      data_tz_codes <- bind_rows(data_tz_codes, default_tz_codes)
    }
  }
  
  if(device_type == "fitbit"){
    if(!ALLOW_MULTIPLE_TZ_PER_DEVICE)
      validate_single_tz_per_fitbit_device(data_tz_codes, INFER_FROM_SMARTPHONE_TZ)
    # We only use datetimes for raw Fitbit data  
    data_tz_codes <- data_tz_codes %>% 
      group_by(tzcode) %>% 
      nest() %>% 
      mutate(data = map2(data, tzcode, function(nested_data, tz){
        nested_data %>% mutate(local_date_time = format(as_datetime(timestamp / 1000, tz=tz), format="%Y-%m-%d %H:%M:%S"))  
      })) %>% 
      unnest(cols=everything()) %>% 
      ungroup()
  }
  
  tz_intervals <- buils_tz_intervals(data_tz_codes, device_type)
  sensor_data <- sensor_data %>% mutate(local_timezone = NA_character_)
  
  if(nrow(sensor_data) > 0){
    sensor_data <- sensor_data %>%
      group_by(device_id) %>% 
      nest() %>% 
      mutate(data = map2(data, device_id, assign_tz_code, tz_intervals, device_type)) %>% 
      unnest(cols = data) %>% 
      ungroup()
  }
  
  return(sensor_data)
}