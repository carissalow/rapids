source("renv/activate.R")
library("lubridate")
library("readr")
library("dplyr")
library("tidyr")
library("stringr")
library("yaml")

validate_periodic_segments <- function(segments){
  invalid_lengths <- segments %>% mutate(is_valid = str_detect(length, "^[[:space:]]*(\\d+?[d|D])??[[:space:]]*(\\d+?[h|H])??[[:space:]]*(\\d+?[m|M])??[[:space:]]*(\\d+?[s|S])??$"))
  if(any(!(invalid_lengths$is_valid)))
     stop("One or more rows in your periodic time segments file have an invalid length format (XXD XXH XXM XXS): ", 
          paste(invalid_lengths %>% filter(!is_valid) %>% pull(label), collapse = ", "))
     
  if(any(is.na(segments$length_period)))
    stop("One or more rows in your periodic time segments file have an invalid length value: ", 
         paste(segments %>% filter(is.na(length_period)) %>% pull(label), collapse = ","))
  
  if(any(is.na(segments$start_time_format)))
    stop("One or more rows in your periodic time segments file have an invalid start_time (HH:MM:SS): ", 
         paste(segments %>% filter(is.na(start_time_format)) %>% pull(label), collapse = ", "))
  
  longer_start_time <- segments %>% mutate(is_longer = start_time_format > period("23H 59M 59S"))
  if(any(longer_start_time$is_longer))
    stop("One or more rows in your periodic time segments file have a start_time longer than 23:59:59: ", 
         paste(longer_start_time %>% filter(is_longer) %>% pull(label), collapse = ", "))
  
  invalid_repeats_on <- segments %>% filter(!repeats_on %in% c("every_day", "wday", "mday", "qday","yday")) %>% pull(label)
  if(length(invalid_repeats_on) > 0)
    stop("One or more rows in your periodic time segments file have an invalid repeats_on: ", 
         paste(invalid_repeats_on, collapse = ","), 
         ". Valid values include: ", 
         paste(c("every_day", "wday", "mday", "qday","yday"), collapse = ", "))
  
  if(nrow(count(segments, label) %>% filter(n > 1)) > 0)
    stop("The values in the column 'label' should be unique but they are not: ", 
         paste(count(segments, label) %>% filter(n > 1) %>% pull(label), collapse = ", "), 
         ". Valid values include: ", 
         paste(c("every_day", "wday", "mday", "qday","yday"), collapse = ", "))
  
  if(nrow(filter(segments, length_period > repeats_on_period & repeats_on %in% c("mday", "qday", "yday"))))
  stop("We do not support mday, qday, or yday segments that overlap yet. Get in touch with the RAPIDS team if you'd like to have this functionality. Overlapping segments: ",
       paste((filter(segments, length_period > repeats_on_period)) %>% filter(repeats_on %in% c("mday", "qday", "yday")) %>% pull(label), collapse = ","))
  
  distinct_segments <- segments %>% distinct(across(-label), .keep_all=TRUE)
  if(nrow(segments) != nrow(distinct_segments))
    stop("Your periodic time segments file has ", nrow(segments) - nrow(distinct_segments), " duplicated row(s) (excluding label): ",
         paste(setdiff(segments %>% pull(label), distinct_segments %>% pull(label)), collapse = ","))
  
  invalid_repeats_value <-  segments %>% 
    mutate(is_invalid = case_when(repeats_on == "every_day" ~ repeats_value != 0,
                                  repeats_on == "wday" ~ repeats_value < 1 | repeats_value > 7,
                                  repeats_on == "mday" ~ repeats_value < 1 | repeats_value > 31,
                                  repeats_on == "qday" ~ repeats_value < 1 | repeats_value > 91,
                                  repeats_on == "yday" ~ repeats_value < 1 | repeats_value > 365))
  if(any(invalid_repeats_value$is_invalid))
    stop("One or more rows in your periodic time segments file have an invalid repeats_value (0 for every_day, [1,7] for wday, [1,31] for mday, [1,91] for qday, [1,366] for yday): ", 
         paste(invalid_repeats_value %>% filter(is_invalid) %>% pull(label), collapse = ", "))
  return(segments)

}

validate_periodic_columns <- function(segments){
  if(nrow(segments) == 0)
    stop("Your periodic time segments file is empty: ", segments_file)
  
  if(!identical(colnames(segments), c("label","start_time","length","repeats_on","repeats_value")))
    stop("Your periodic time segments file does not have the expected columns (label,start_time,length,repeats_on,repeats_value). Maybe you have a typo in the names?")
  return(segments)
}

prepare_periodic_segments <- function(segments){
  segments <- segments %>% 
    validate_periodic_columns() %>% 
    mutate(length_period = period(length),
           start_time_format = hms(start_time, quiet = TRUE),
           repeats_on_period = case_when(repeats_on == "every_day" ~ period("1D"),
                                         repeats_on == "wday" ~ period("7D"),
                                         repeats_on == "mday" ~ period("28D"),
                                         repeats_on == "qday" ~ period("95D"),
                                         repeats_on == "yday" ~ period("365D"))) %>% 
    validate_periodic_segments() %>% 
    mutate(new_segments = (length_period %/% repeats_on_period) + 1) %>% 
    uncount(weights = new_segments, .remove = FALSE, .id = "overlap_id") %>% 
    mutate(overlap_id = overlap_id -1,
           original_label = label,
           overlap_duration = paste0(overlap_id * repeats_on_period / days(1),"D"),
           label = paste0(label, "_RR", overlap_id, "SS")) %>% 
    select(label,start_time,length,repeats_on,repeats_value,overlap_duration,overlap_id,original_label)
  return(segments)
}

validate_frequency_segments <- function(segments){
  if(nrow(segments) == 0)
    stop("Your frequency time segments file is empty: ", segments_file)
  if(!identical(colnames(segments), c("label","length")))
    stop("Your frequency time segments file does not have the expected columns (label, length). Maybe you have a typo in the names?")
  if(nrow(segments) > 1)
    stop("Your frequency time segments file cannot have more than one row")
  if(any(is.na(segments$label)))
    stop("Your frequency time segments file has an empty or invalid label")
  if(nrow(segments %>% filter(!is.na(length) & length >= 1 & length <= 1440)) == 0)
    stop("Your frequency time segments file has an empty or invalid length (only numbers between [1,1440] are accepted), you typed: ", segments$length)
  return(segments)
}

prepare_frequency_segments <- function(segments){
  #FREQUENCY segments are just syntactic sugar for PERIODIC
  validate_frequency_segments(segments)
  stamp_fn <- stamp("23:10:00", orders = c("HMS"), quiet = TRUE)
  new_segments <- data.frame(start_time = seq.POSIXt(from = ymd_hms("2020-01-01 00:00:00"), 
                                                        to=ymd_hms("2020-01-02 00:00:00"), 
                                                        by=paste(segments$length, "min")))
  new_segments <- new_segments %>%
    head(-1) %>% 
    mutate(label = paste0(segments$label, str_pad(row_number()-1, width = 4, pad = "0")),
           start_time = stamp_fn(start_time),
           length = paste0((segments$length * 60)-1, "S"),
           repeats_on = "every_day",
           repeats_value=0,
           overlap_id = 0,
           original_label = label,
           overlap_duration = "0D")
    
}

get_devices_ids <- function(participant_data){
  devices_ids = c()
  for(device in participant_data)
    for(attribute in names(device))
      if(attribute == "DEVICE_IDS")
        devices_ids <- c(devices_ids, device[[attribute]])
  return(devices_ids)
}

validate_event_segments <- function(segments){
  if(nrow(segments) == 0)
    stop("The following time segments file is empty: ", segments_file)
  
  if(!identical(colnames(segments), c("label","event_timestamp","length","shift","shift_direction","device_id")))
    stop("Your periodic time segments file does not have the expected columns (label,event_timestamp,length,shift,shift_direction,device_id). Maybe you have a typo in the names?")
  
  invalid_lengths <- segments %>% mutate(is_valid = str_detect(length, "^[[:space:]]*(\\d+?[d|D])??[[:space:]]*(\\d+?[h|H])??[[:space:]]*(\\d+?[m|M])??[[:space:]]*(\\d+?[s|S])??$"))
  if(any(!(invalid_lengths$is_valid)))
    stop("One or more rows in your event time segments file have an invalid length format (XXD XXH XXM XXS): ", 
         paste(invalid_lengths %>% filter(!is_valid) %>% pull(label), collapse = ", "))
  
  invalid_shifts <- segments %>% mutate(is_valid = str_detect(shift, "^[[:space:]]*(\\d+?[d|D])??[[:space:]]*(\\d+?[h|H])??[[:space:]]*(\\d+?[m|M])??[[:space:]]*(\\d+?[s|S])??$"))
  if(any(!(invalid_shifts$is_valid)))
    stop("One or more rows in your event time segments file have an invalid shift format (XXD XXH XXM XXS): ", 
         paste(invalid_shifts %>% filter(!is_valid) %>% pull(label), collapse = ", "))
  
  invalid_shift_direction <- segments %>% filter(shift_direction < -1 | shift_direction > 1)
  if(nrow(invalid_shift_direction) > 0)
    stop("One or more rows in your event time segments file have an invalid shift direction (-1,0,1): ", 
         paste(invalid_shift_direction %>% pull(label), collapse = ", "))
  
  invalid_timestamps <- segments %>% filter(is.na(event_timestamp))
  if(nrow(invalid_timestamps) > 0)
    stop("One or more rows in your event time segments file have an empty timestamp: ", 
         paste(invalid_timestamps %>% pull(label), collapse = ", "))
  
  invalid_timestamps <- segments %>% filter(event_timestamp <= 999999999999)
  if(nrow(invalid_timestamps) > 0)
    stop("One or more rows in your event time segments file is not in milliseconds: ", 
         paste(invalid_timestamps %>% pull(label), collapse = ", "))
  
  distinct_segments <- segments %>% mutate(row_id = row_number()) %>% distinct(across(c(-label, -row_id)), .keep_all=TRUE)
  if(nrow(segments) != nrow(distinct_segments))
    stop("Your event time segments file has ", nrow(segments) - nrow(distinct_segments), " duplicated row(s) (excluding label). Duplicated row number(s): ",
         paste(setdiff(segments %>% mutate(row_id = row_number()) %>% pull(row_id), distinct_segments %>% pull(row_id)), collapse = ","))
    
  return(segments)
}

prepare_event_segments <- function(segments, participant_devices){
  new_segments <- segments%>%
    validate_event_segments() %>% 
    filter(device_id %in% participant_devices)
}

compute_time_segments <- function(){
  type = snakemake@params[["time_segments_type"]]
  pid = snakemake@params[["pid"]]
  segments_file <- snakemake@input[["segments_file"]]
  participant_file <- snakemake@input[["participant_file"]]
  message("Processing ",type, " time segments for ", pid,"'s ", participant_file)
  
  participant_data <- yaml::read_yaml(participant_file)
  participant_devices <- get_devices_ids(participant_data)
  if(length(participant_devices) == 0)
    stop("There are no device ids in this participant file for smartphones or wearables: ", participant_file)
  
  if(type == "FREQUENCY"){
    segments <- read_csv(segments_file, col_types = cols_only(label = "c", length = "i"), trim_ws = TRUE)
    new_segments <- prepare_frequency_segments(segments)
  } else if(type == "PERIODIC"){
    segments <- read_csv(segments_file, col_types = cols_only(label = "c", start_time = "c",length = "c",repeats_on = "c",repeats_value = "i"), trim_ws = TRUE)
    new_segments <- prepare_periodic_segments(segments)
  } else if(type == "EVENT"){
    segments <- read_csv(segments_file, col_types = cols_only(label = "c", event_timestamp = "d",length = "c",shift = "c",shift_direction = "i", device_id = "c"), trim_ws = TRUE)
    new_segments <- prepare_event_segments(segments, participant_devices)
  }
  
  write.csv(new_segments %>% select(label) %>% distinct(label), snakemake@output[["segments_labels_file"]], row.names = FALSE, quote = FALSE)
  write.csv(new_segments,snakemake@output[["segments_file"]], row.names = FALSE, quote = FALSE)
}

compute_time_segments()