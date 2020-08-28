source("renv/activate.R")
library("tidyverse")
library("readr")

source("src/data/assign_to_day_segment.R")

input <- read.csv(snakemake@input[["sensor_input"]]) %>% arrange(timestamp)
day_segments <- read.csv(snakemake@input[["day_segments"]])
day_segments_type <- snakemake@params[["day_segments_type"]]
sensor_output <- snakemake@output[[1]]
timezone_periods <- snakemake@params[["timezone_periods"]]
fixed_timezone <- snakemake@params[["fixed_timezone"]]

split_local_date_time <- function(data, day_segments){
  split_data <- data %>% 
    separate(local_date_time, c("local_date","local_time"), "\\s", remove = FALSE) %>%
    separate(local_time, c("local_hour", "local_minute"), ":", remove = FALSE, extra = "drop") %>%
    mutate(local_hour = as.numeric(local_hour),
           local_minute = as.numeric(local_minute))

  return(split_data)
}

if(!is.null(timezone_periods)){
  # TODO: Not active yet
  # timezones <- read_csv(timezone_periods)
  # tz_starts <- timezones$start
  # output <- input %>% 
  #   mutate(timezone = findInterval(timestamp / 1000, tz_starts), # Set an interval ID based on timezones' start column
  #          timezone = ifelse(timezone == 0, 1, timezone), # Correct the first timezone ID
  #          timezone = recode(timezone, !!! timezones$timezone), # Swap IDs for text labels
  #          timezone = as.character(timezone)) %>%
  #   rowwise() %>%
  #   mutate(utc_date_time = as.POSIXct(timestamp/1000, origin="1970-01-01", tz="UTC"),
  #          local_date_time = format(utc_date_time, tz = timezone, usetz = T, "%Y-%m-%d %H:%M:%S"))
  # output <- split_local_date_time(output, day_segments)
  # TODO: Implement day segment assigment with support for multiple timezones
  # output <- assign_to_day_segment(output, day_segments, day_segments_type, fixed_timezone)
  # write.csv(output, sensor_output)
} else if(!is.null(fixed_timezone)){
  output <- input %>% 
    mutate(utc_date_time = as.POSIXct(timestamp/1000, origin="1970-01-01", tz="UTC"),
           local_date_time = format(utc_date_time, tz = fixed_timezone, usetz = F, "%Y-%m-%d %H:%M:%S"))
  output <- split_local_date_time(output, day_segments)
  output <- assign_to_day_segment(output, day_segments, day_segments_type, fixed_timezone)
  write_csv(output, sensor_output)
}
