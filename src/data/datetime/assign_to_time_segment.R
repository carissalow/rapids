library("tidyverse")
library("glue")
library("lubridate", warn.conflicts = F)
options(scipen=999)

assign_rows_to_segments <- function(data, segments){
  # This function is used by all segment types, we use data.tables because they are fast
  data <- data.table::as.data.table(data)
  data[, assigned_segments := ""]
  for(i in seq_len(nrow(segments))) {
    segment <- segments[i,]
    data[segment$segment_start_ts<= timestamp & segment$segment_end_ts >= timestamp,
         assigned_segments := stringi::stri_c(assigned_segments, segment$segment_id, sep = "|")]
  }
  data[,assigned_segments:=substring(assigned_segments, 2)]
  data
}

assign_to_time_segment <- function(sensor_data, time_segments, time_segments_type, include_past_periodic_segments){
  
  if(nrow(sensor_data) == 0 || nrow(time_segments) == 0)
    return(sensor_data %>% mutate(assigned_segments = NA))
  
  if (time_segments_type == "FREQUENCY" || time_segments_type == "PERIODIC"){ #FREQUENCY segments are just syntactic sugar for PERIODIC
    source("src/data/datetime/assign_to_periodic_segments.R")
    sensor_data <- assign_to_periodic_segments(sensor_data, time_segments, include_past_periodic_segments)
    return(sensor_data)
    
  } else if ( time_segments_type == "EVENT"){
    source("src/data/datetime/assign_to_event_segments.R")
    sensor_data <- assign_to_event_segments(sensor_data, time_segments)
    return(sensor_data)
  }
}