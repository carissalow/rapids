library("tidyverse")
library("lubridate")

find_segments_frequency <- function(local_date, local_time_obj, segments){
  return(paste(segments %>% 
                 mutate(in_segment = local_time_obj >= segment_start & local_time_obj <= segment_end) %>% 
                 filter(in_segment == TRUE) %>% 
                 mutate(segment_id =  paste0("[", 
                                                label, "#",
                                                local_date, "#",
                                                paste(str_pad(hour(segment_start),2, pad="0"), str_pad(minute(segment_start),2, pad="0"), str_pad(second(segment_start),2, pad="0"),sep =":"), "#",
                                                local_date, "#",
                                                paste(str_pad(hour(segment_end),2, pad="0"), str_pad(minute(segment_end),2, pad="0"), str_pad(second(segment_end),2, pad="0"),sep =":"),
                                                "]")) %>% 
                 pull(segment_id), collapse = "|"))
}

find_segments_periodic <- function(date_time, segments){
  return(paste(segments[[1]] %>% 
                 select(segment_interval, segment_id) %>% 
                 mutate(in_segment = date_time %within% segment_interval) %>% 
                 filter(in_segment == TRUE) %>% 
                 pull(segment_id), collapse = "|"))
}

find_segments_event <- function(timestamp, segments){
  return(paste(segments %>% 
                 mutate(in_segment = timestamp >= segment_start & timestamp <= segment_end) %>% 
                 filter(in_segment == TRUE) %>% 
                 pull(segment_id), collapse = "|"))
}

assign_to_day_segment <- function(sensor_data, day_segments, day_segments_type){

  if(day_segments_type == "FREQUENCY"){ #FREQUENCY
    
    day_segments <- day_segments %>% mutate(segment_start = lubridate::parse_date_time(start_time, orders = c("HMS", "HM")),
                                            segment_end = segment_start + minutes(length))
    sensor_data <- sensor_data %>% mutate(local_time_obj = lubridate::parse_date_time(local_time, orders = c("HMS", "HM")),
                                          assigned_segments = map2_chr(local_date, local_time_obj, ~find_segments_frequency(.x, .y, day_segments))) %>% select(-local_time_obj)
    
  } else if (day_segments_type == "PERIODIC"){ #PERIODIC
    
    sensor_data <- sensor_data %>%
      mutate(row_n = row_number()) %>% 
      group_by(local_timezone) %>% 
      nest() %>% 
      # get existent days that we need to start segments from
      mutate(existent_dates = map(data, ~.x %>% 
                                    distinct(local_date) %>% 
                                    mutate(local_date_obj = lubridate::ymd(local_date, tz = local_timezone), 
                                           every_day = 0,
                                           wday = wday(local_date_obj, week_start = 1),
                                           mday = mday(local_date_obj),
                                           qday = qday(local_date_obj),
                                           yday = yday(local_date_obj)
                                           ) %>%  select(local_date, every_day, wday, mday, qday, yday)),
             # build the actual day segments taking into account the users requested leangth and repeat schedule
             inferred_day_segments = map(existent_dates,
                                         ~ crossing(day_segments, .x) %>%
                                           pivot_longer(cols = c(every_day,wday, mday, qday, yday), names_to = "day_type", values_to = "day_value") %>%
                                           filter(repeats_on == day_type & repeats_value == day_value) %>%
                                           mutate(segment_start = (lubridate::parse_date_time(paste(local_date, start_time), orders = c("Ymd HMS", "Ymd HM"), tz = local_timezone)),
                                                  segment_end = segment_start + lubridate::period(length),
                                                  segment_interval = lubridate::interval(segment_start, segment_end),
                                                  segment_id = paste0("[",
                                                                         paste(sep= "#",
                                                                               label,
                                                                               lubridate::date(int_start(segment_interval)),
                                                                               paste(str_pad(hour(int_start(segment_interval)),2, pad="0"),
                                                                                     str_pad(minute(int_start(segment_interval)),2, pad="0"),
                                                                                     str_pad(second(int_start(segment_interval)),2, pad="0"),sep =":"),
                                                                               lubridate::date(int_end(segment_interval)),
                                                                               paste(str_pad(hour(int_end(segment_interval)),2, pad="0"),
                                                                                     str_pad(minute(int_end(segment_interval)),2, pad="0"),
                                                                                     str_pad(second(int_end(segment_interval)),2, pad="0"),sep =":")
                                                                         ),
                                                                         "]")) %>%
                                           select(segment_interval, label, segment_id)),
             # loop thorugh every day segment and assigned it to the rows that fall within its start and end
             data = map2(data, inferred_day_segments, ~ .x %>% mutate(row_date_time = lubridate::ymd_hms(local_date_time, tz = local_timezone),
                                                                      assigned_segments = map_chr(row_date_time, ~find_segments_periodic(.x, inferred_day_segments))) %>% 
                                                               select(-row_date_time))
             ) %>% 
      select(-existent_dates, -inferred_day_segments) %>% 
      unnest(cols = data) %>% 
      arrange(row_n) %>%
      select(-row_n)
      
    
  } else if ( day_segments_type == "EVENT"){
    
    most_common_tz <- sensor_data %>% count(local_timezone) %>% slice(which.max(n)) %>% pull(local_timezone)
    day_segments <- day_segments %>% mutate(shift = ifelse(shift == "0", "0seconds", shift),
                                            segment_start = event_timestamp + (as.integer(seconds(lubridate::duration(shift))) * ifelse(shift_direction >= 0, 1, -1) * 1000),
                                            segment_end = segment_start + (as.integer(seconds(lubridate::duration(length))) * 1000),
                                            segment_start_datetime = lubridate::as_datetime(segment_start/1000, tz = most_common_tz), # these start and end datetime objects are for labeling only
                                            segment_end_datetime = lubridate::as_datetime(segment_end/1000, tz = most_common_tz),
                                            segment_id = paste0("[",
                                                                paste(sep= "#",
                                                                      label,
                                                                      lubridate::date(segment_start_datetime),
                                                                      paste(str_pad(hour(segment_start_datetime),2, pad="0"),
                                                                            str_pad(minute(segment_start_datetime),2, pad="0"),
                                                                            str_pad(second(segment_start_datetime),2, pad="0"),sep =":"),
                                                                      lubridate::date(segment_end_datetime),
                                                                      paste(str_pad(hour(segment_end_datetime),2, pad="0"),
                                                                            str_pad(minute(segment_end_datetime),2, pad="0"),
                                                                            str_pad(second(segment_end_datetime),2, pad="0"),sep =":")
                                                                ),
                                                                "]")) %>% 
      select(-segment_start_datetime, -segment_end_datetime)
    
    sensor_data <- sensor_data %>% mutate(assigned_segments = map_chr(timestamp, ~find_segments_event(.x, day_segments)))
  }
  
  return(sensor_data)
}