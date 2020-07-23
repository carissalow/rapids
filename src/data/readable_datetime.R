source("renv/activate.R")

library("tidyverse")
library("readr")
library("lubridate")

input <- read.csv(snakemake@input[["sensor_input"]]) %>% arrange(timestamp)
day_segments <- read.csv(snakemake@input[["day_segments"]]) %>% filter(label != "daily") #daily is done by default by all scripts
sensor_output <- snakemake@output[[1]]
timezone_periods <- snakemake@params[["timezone_periods"]]
fixed_timezone <- snakemake@params[["fixed_timezone"]]

assign_to_day_segment <- function(data, day_segments){
  data <- data %>% mutate(local_day_segment = NA)
  
  # All segments belong to the same date, so we assume all days have the same segments
  if(length(unique(day_segments$local_date)) == 1){ 
    data <- data %>% mutate(local_time_obj = lubridate::hms(local_time))
    day_segments <- day_segments %>% mutate(start_time = lubridate::hm(start_time),
                                           end_time = lubridate::hm(end_time))
    for(row_id in 1:nrow(day_segments)){
      row = day_segments[row_id,]
      data <- data %>% mutate(local_day_segment = ifelse(local_time_obj >= row$start_time & local_time_obj <= row$end_time, row$label, local_day_segment))
    }
    data <- data %>% select(-local_time_obj)
  # Segments belong to different dates, so each day can have different segments
  }else{ 
    data <- data %>% mutate(local_date_time_obj = lubridate::ymd_hms(local_date_time))
    day_segments <- day_segments %>% mutate(start_local_date_time_obj = lubridate::ymd_hm(paste(local_date, start_time)),
                                            end_local_date_time_obj = lubridate::ymd_hm(paste(local_date, end_time)),
                                            date_time_interval = lubridate::interval(start_local_date_time_obj, end_local_date_time_obj))
    for(row_id in 1:nrow(day_segments)){
      row = day_segments[row_id,]
      data <- data %>% mutate(local_day_segment = ifelse(local_date_time_obj %within% row$date_time_interval, row$label, local_day_segment))
    }
    data <- data %>% select(-local_date_time_obj)
  }
  
  return(data)
}

split_local_date_time <- function(data, day_segments){
  split_data <- data %>% 
    separate(local_date_time, c("local_date","local_time"), "\\s", remove = FALSE) %>%
    separate(local_time, c("local_hour", "local_minute"), ":", remove = FALSE, extra = "drop") %>%
    mutate(local_hour = as.numeric(local_hour),
           local_minute = as.numeric(local_minute))
  
  split_data <- assign_to_day_segment(split_data, day_segments)
  return(split_data)
}

if(!is.null(timezone_periods)){
    timezones <- read_csv(timezone_periods)
    tz_starts <- timezones$start
    output <- input %>% 
                mutate(timezone = findInterval(timestamp / 1000, tz_starts), # Set an interval ID based on timezones' start column
                        timezone = ifelse(timezone == 0, 1, timezone), # Correct the first timezone ID
                        timezone = recode(timezone, !!! timezones$timezone), # Swap IDs for text labels
                        timezone = as.character(timezone)) %>%
                rowwise() %>%
                mutate(utc_date_time = as.POSIXct(timestamp/1000, origin="1970-01-01", tz="UTC"),
                        local_date_time = format(utc_date_time, tz = timezone, usetz = T))
    output <- split_local_date_time(output, day_segments)
    write.csv(output, sensor_output)
} else if(!is.null(fixed_timezone)){
    output <- input %>% 
                mutate(utc_date_time = as.POSIXct(timestamp/1000, origin="1970-01-01", tz="UTC"),
                        local_date_time = format(utc_date_time, tz = fixed_timezone, usetz = F))
    output <- split_local_date_time(output, day_segments)
    write_csv(output, sensor_output)
}
