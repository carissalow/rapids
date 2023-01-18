source("renv/activate.R")
library("tidyverse")
library("readr")
library("tidyr")
library("lubridate")
library("yaml")

source("src/data/datetime/assign_to_time_segment.R")
source("src/data/datetime/assign_to_multiple_timezones.R")

split_local_date_time <- function(data){
  data <- data %>% 
    separate(local_date_time, c("local_date","local_time"), "\\s", remove = FALSE) %>%
    separate(local_time, c("local_hour", "local_minute"), ":", remove = FALSE, extra = "drop") %>%
    mutate(local_hour = as.numeric(local_hour),
           local_minute = as.numeric(local_minute))

  return(data)
}

is_valid_timezone <- function(timezone) {
  return(timezone %in% (OlsonNames()))
}

validate_user_timezones <- function(timezone_parameters){
  if(!timezone_parameters$TYPE %in% c("SINGLE", "MULTIPLE"))
    stop("Invalid [TIMEZONE][TYPE], only valid options are SINGLE or MULTIPLE")

  if(timezone_parameters$TYPE == "SINGLE"){
    if(!is_valid_timezone(timezone_parameters$SINGLE$TZCODE))
      stop(paste("[TIMEZONE][SINGLE][TZCODE] is not a valid timezone: ", timezone_parameters$SINGLE$TZCODE))
  } else if(timezone_parameters$TYPE == "MULTIPLE"){
    tz_codes <- read.csv(timezone_parameters$MULTIPLE$TZCODES_FILE)
    valid_file_columns <- c("device_id", "timestamp", "tzcode")

    if(length(colnames(tz_codes)) != length(valid_file_columns) || !setequal(colnames(tz_codes), valid_file_columns))
      stop(paste("[TIMEZONE][MULTIPLE][TZCODES_FILE] has does not have the required columns. You provided",paste(colnames(tz_codes), collapse=","),"but we need",paste(valid_file_columns, collapse=",")))
    
    invalid_tz_codes <- tz_codes %>% 
      mutate(row = (1:n()) + 1, 
              tzcode = trimws(tzcode, which="both"),
              is_valid = is_valid_timezone(tzcode)) %>% 
      filter(is_valid == FALSE)
    if(nrow(invalid_tz_codes) > 0)
      stop(paste("[TIMEZONE][MULTIPLE][TZCODES_FILE] has invalid time zone codes. In file ", timezone_parameters$MULTIPLE$TZCODES_FILE, ".\nAffected rows=[", paste(invalid_tz_codes %>% pull(row),collapse=","), "], with invalid codes=[", paste(invalid_tz_codes %>% pull(tzcode),collapse=",") ,"]"))
    

  }
}

create_mising_temporal_column <- function(data, device_type){
  if(device_type == "fitbit" && all(data$timestamp == 0) && "local_date_time" %in% colnames(data)){
      # For fibit we infere timestamp from Fitbit's local date time only right after pulling it
      if(nrow(data) == 0)
        return(data %>% mutate(timestamp = NA_real_) %>% relocate(local_timezone)) # move local_timezone to the beginning of df to match column positions when data is not empty (group > nest)
      if(any(is.na(parse_date_time(data$local_date_time, orders= c("%Y/%m/%d %H:%M:%S","%Y-%m-%d %H:%M:%S"), exact=T))))
        stop("One or more values in the local_date_time column do not have the expected format: yyyy-mm-dd hh:mm:ss or yyyy/mm/dd hh:mm:ss")
      return(data %>% 
          group_by(local_timezone) %>% 
          nest() %>% 
          mutate(data = map2(data, local_timezone, function(nested_data, tz){
            return(nested_data %>%  mutate(timestamp = as.numeric(ymd_hms(local_date_time, tz=tz)) * 1000) %>% drop_na(timestamp))
            })) %>% 
          unnest(cols = everything()) %>% 
          ungroup())
    } else {
      # For the rest of devices we infere local date time from timestamp
      if(nrow(data) == 0)
        return(data %>% mutate(local_date_time = NA_character_) %>% relocate(local_timezone)) # move local_timezone to the beginning of df to match column positions when data is not empty (group > nest)
      return(data %>% 
          group_by(local_timezone) %>% 
          nest() %>% 
          mutate(data = map2(data, local_timezone, function(nested_data, tz){
            return(nested_data %>%  mutate(local_date_time = format(as_datetime(timestamp / 1000, tz=tz), format="%Y-%m-%d %H:%M:%S")) %>% drop_na(local_date_time) )
            })) %>% 
          unnest(cols = everything()) %>% 
          ungroup())
    }
}

filter_wanted_dates <- function(output, participant_file, device_type){
  participant_data <- read_yaml(participant_file)
  device_type <- toupper(device_type)
  start_date <- participant_data[[device_type]]$START_DATE
  end_date <- participant_data[[device_type]]$END_DATE

  if(!is.null(start_date)){
    start_date <- parse_date_time(start_date, orders = c("ymd", "ymdhMs", "ymdhM", "ymdh"))
    if(is.na(start_date))
      stop(paste0("[",device_type, "][START_DATE] does not have one of these valid formats: [ymd, ymd hms, ymd hm, ymd h], you typed: '", participant_data[[device_type]]$START_DATE, "' in ", participant_file))
    output <- output %>% filter(ymd_hms(local_date_time) >= start_date)
  }
  if(!is.null(end_date)){
    end_date <- parse_date_time(end_date, orders = c("ymd", "ymdhMs", "ymdhM", "ymdh"))
    if(is.na(end_date))
      stop(paste0("[",device_type, "][END_DATE] does not have one of these valid formats: [ymd, ymd hms, ymd hm, ymd h], you typed: '", participant_data[[device_type]]$END_DATE, "' in ", participant_file))
    output <- output %>% filter(ymd_hms(local_date_time) <= end_date)
  }

  return(output)
}

process_by_chunks <- function(timezone_parameters, device_type, pid, participant_file, time_segments, time_segments_type, include_past_periodic_segments, output_file_path){
  
  function(x, pos, acc){

    if(timezone_parameters$TYPE == "SINGLE"){
      output <- x %>% mutate(local_timezone = timezone_parameters$SINGLE$TZCODE)
      most_common_tz <- timezone_parameters$SINGLE$TZCODE
    }
    else if(timezone_parameters$TYPE == "MULTIPLE"){
      output <- multiple_time_zone_assignment(x, timezone_parameters, device_type, pid, participant_file)
      most_common_tz <- get_participant_most_common_tz(timezone_parameters$MULTIPLE$TZCODES_FILE, participant_file) # in assign_to_multiple_timezones.R
    }

    output %<>%  
      create_mising_temporal_column(device_type)  %>% 
      split_local_date_time() %>% 
      assign_to_time_segment(time_segments, time_segments_type, include_past_periodic_segments, most_common_tz) %>% 
      filter_wanted_dates(participant_file, device_type) %>% 
      arrange(timestamp)
    
    output %>% write_csv(output_file_path, append=ifelse(pos > 1, T, F))
    
    return(acc + 1)
  }
}

readable_datetime <- function(){

  time_segments <- read.csv(snakemake@input[["time_segments"]])
  sensor_input <- snakemake@input[["sensor_input"]]
  participant_file <- snakemake@input[["pid_file"]]
  device_type <- snakemake@params[["device_type"]]
  timezone_parameters <- snakemake@params[["timezone_parameters"]]
  pid <- snakemake@params[["pid"]]
  time_segments_type <- snakemake@params[["time_segments_type"]]
  include_past_periodic_segments <- snakemake@params[["include_past_periodic_segments"]]
  output_file_path <- snakemake@output[[1]]

  if (device_type != "fitbit") {
    column_specifications <- cols()
  } else if (!grepl("sleep_summary", sensor_input)) {
    column_specifications <- cols(local_date_time = col_character())
  } else {
    column_specifications <- cols(
      local_date_time = col_character(),
      local_start_date_time = col_character(),
      local_end_date_time = col_character()
    )
  } 

  validate_user_timezones(timezone_parameters)

  acc <- read_csv_chunked(sensor_input, 
                  callback=AccumulateCallback$new(process_by_chunks(timezone_parameters, device_type, pid, participant_file, time_segments, time_segments_type, include_past_periodic_segments, output_file_path), acc=0), 
                  chunk_size=10000,
                  col_names=TRUE,
                  col_types=column_specifications,
                  trim_ws=FALSE)
  
  if(acc == 0){
    column_names <- c("local_timezone", colnames(read.csv(sensor_input)), "local_date", "local_time", "local_hour", "local_minute", "assigned_segments")
    output <- setNames(data.frame(matrix(nrow=0, ncol=length(column_names))), column_names)
    output %>% write_csv(output_file_path, append=FALSE)
  }

}

readable_datetime()