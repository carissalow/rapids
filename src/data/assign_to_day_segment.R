library("tidyverse")
library("lubridate")

assign_to_day_segment <- function(data, day_segments, day_segments_type, fixed_timezone){
  
  if(day_segments_type == "FREQUENCY_EVERY_DAY"){
    data <- data %>% mutate(local_date_time_obj = lubridate::parse_date_time(local_time, orders = c("HMS", "HM")))
    day_segments <- day_segments %>% mutate(start_time = lubridate::parse_date_time(start_time, orders = c("HMS", "HM")),
                                            end_time = start_time + minutes(length))
    
    # Create a new column for each day_segment
    for(row_id in 1:nrow(day_segments)){
      row = day_segments[row_id,]
      data <- data %>% mutate(!!paste("local_day_segment", row_id, sep = "_") := ifelse(local_date_time_obj >= row$start_time & local_date_time_obj < row$end_time, 
                                                                                        paste0("[", 
                                                                                               row$label, "_", 
                                                                                               local_date, "_",
                                                                                               paste(str_pad(hour(row$start_time),2, pad="0"), str_pad(minute(row$start_time),2, pad="0"), str_pad(second(row$start_time),2, pad="0"),sep =":"),
                                                                                               "]"), NA))
    }
    
  } else if (day_segments_type == "INTERVAL_EVERY_DAY"){
    
    data_dates <- data %>% select(local_date) %>% distinct(local_date)
    inferred_day_segments <- crossing(day_segments, data_dates) %>% 
      mutate(start_local_date_time_obj = lubridate::parse_date_time(paste(local_date, start_time), orders = c("Ymd HMS", "Ymd HM"), tz = fixed_timezone),
             end_local_date_time_obj = start_local_date_time_obj + lubridate::period(length),
             date_time_interval = lubridate::interval(start_local_date_time_obj, end_local_date_time_obj)) %>% 
      group_by(label, local_date) %>% 
      mutate(group_start_datetime = lubridate::parse_date_time(paste(local_date, start_time), orders = c("Ymd HMS", "Ymd HM"), tz = fixed_timezone),
             group_end_datetime = group_start_datetime + lubridate::period(length),
             group_start_datetime = min(group_start_datetime),
             group_end_datetime = max(group_end_datetime)) %>% 
      ungroup()

    
    data <- data %>% mutate(local_date_time_obj = lubridate::ymd_hms(local_date_time, tz = fixed_timezone))
    
    # Create a new column for each day_segment
    for(row_id in 1:nrow(inferred_day_segments)){
      row = inferred_day_segments[row_id,]
      data <- data %>% mutate(!!paste("local_day_segment", row_id, sep = "_") := ifelse(local_date_time_obj %within% row$date_time_interval, 
                                                                                        paste0("[", 
                                                                                               paste(sep= "#",
                                                                                                     row$label,
                                                                                                     lubridate::date(row$group_start_datetime),
                                                                                                     paste(str_pad(hour(row$group_start_datetime),2, pad="0"), str_pad(minute(row$group_start_datetime),2, pad="0"), str_pad(second(row$group_start_datetime),2, pad="0"),sep =":"),
                                                                                                     lubridate::date(row$group_end_datetime),
                                                                                                     paste(str_pad(hour(row$group_end_datetime),2, pad="0"), str_pad(minute(row$group_end_datetime),2, pad="0"), str_pad(second(row$group_end_datetime),2, pad="0"),sep =":")
                                                                                               ),
                                                                                               "]"), NA))
    }
    
  
  } else if ( day_segments_type == "INTERVAL_FLEXIBLE_DAY"){
    data <- data %>% mutate(local_date_time_obj = lubridate::ymd_hms(local_date_time, tz = fixed_timezone))
    day_segments <- day_segments %>% mutate(shift = ifelse(shift == "0", "0seconds", shift),
                                            start_local_date_time_obj = lubridate::ymd_hms(start_date_time, tz = fixed_timezone) + (lubridate::period(shift) * ifelse(shift_direction >= 0, 1, -1)),
                                            end_local_date_time_obj = start_local_date_time_obj + lubridate::period(length),
                                            date_time_interval = lubridate::interval(start_local_date_time_obj, end_local_date_time_obj))
    
    # Create a new column for each day_segment
    for(row_id in 1:nrow(day_segments)){
      row = day_segments[row_id,]
      print(row$length)
      data <- data %>% mutate(!!paste("local_day_segment", row_id, sep = "_") := ifelse(local_date_time_obj %within% row$date_time_interval, 
                                                                                        paste0("[", 
                                                                                              paste(sep= "#",
                                                                                                row$label,
                                                                                                lubridate::date(row$start_local_date_time_obj),
                                                                                                paste(str_pad(hour(row$start_local_date_time_obj),2, pad="0"), str_pad(minute(row$start_local_date_time_obj),2, pad="0"), str_pad(second(row$start_local_date_time_obj),2, pad="0"),sep =":"),
                                                                                                lubridate::date(row$end_local_date_time_obj),
                                                                                                paste(str_pad(hour(row$end_local_date_time_obj),2, pad="0"), str_pad(minute(row$end_local_date_time_obj),2, pad="0"), str_pad(second(row$end_local_date_time_obj),2, pad="0"),sep =":")
                                                                                                ),
                                                                                              "]"), NA))
    }
  }
  
  # Join all day_segments in a single column
  data <- data %>% 
    unite("assigned_segments", starts_with("local_day_segment"), sep = "|", na.rm = TRUE) %>% 
    select(-local_date_time_obj)

  return(data)
}