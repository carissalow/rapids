library("dplyr", warn.conflicts = F)
library(tidyr)
library(readr)

compute_data_yield_features <- function(data, feature_name, time_segment, provider){
  data <- data %>% filter_data_by_segment(time_segment)
  features <- data %>%
    separate(timestamps_segment, into = c("start_timestamp", "end_timestamp"), convert = T, sep = ",") %>% 
    mutate(duration_minutes = (end_timestamp - start_timestamp) / 60000,
           timestamp_since_segment_start = timestamp - start_timestamp,
           minute_bin = timestamp_since_segment_start %/% 60000, # 60 * 1000
           hour_bin = timestamp_since_segment_start %/% 3600000) %>% # (60 * 60 * 1000)
    group_by(local_segment, hour_bin) %>% 
    summarise(minute_count = n_distinct(minute_bin),
              duration_minutes = first(duration_minutes),
              valid_hour = (minute_count/min(duration_minutes, 60)) > provider$MINUTE_RATIO_THRESHOLD_FOR_VALID_YIELDED_HOURS) %>% 
    group_by(local_segment) %>% 
    summarise(valid_yielded_minutes = sum(minute_count),
              valid_yielded_hours = sum(valid_hour == TRUE) / 1.0,
              duration_minutes = first(duration_minutes),
              duration_hours = duration_minutes / 60.0,
              ratiovalidyieldedminutes = min( valid_yielded_minutes / duration_minutes, 1),
              ratiovalidyieldedhours = if_else(duration_hours > 1, min( valid_yielded_hours / duration_hours, 1), valid_yielded_hours))
  return(features)
}



rapids_features <- function(sensor_data_files, time_segment, provider){
  
  yield_data <-  read_csv(sensor_data_files[["sensor_data"]], col_types = cols_only(timestamp ="d", assigned_segments = "c"))
  requested_features <- provider[["FEATURES"]]
  
  # Output dataframe
  features = data.frame(local_segment = character(), stringsAsFactors = FALSE)
  
  # The name of the features this function can compute
  base_features_names  <- c("ratiovalidyieldedminutes", "ratiovalidyieldedhours")
  
  # The subset of requested features this function can compute
  features_to_compute  <- intersect(base_features_names, requested_features)
  
  features <- compute_data_yield_features(yield_data, feature_name, time_segment, provider) %>% 
    select(c("local_segment", features_to_compute))
  
  return(features)
}