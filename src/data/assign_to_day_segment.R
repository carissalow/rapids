library("tidyverse")
library("lubridate")

assign_to_day_segment <- function(sensor_data, day_segments, day_segments_type){

  if(day_segments_type == "FREQUENCY"){ #FREQUENCY
    sensor_data <- sensor_data %>% mutate(local_date_time_obj = lubridate::parse_date_time(local_time, orders = c("HMS", "HM")))
    day_segments <- day_segments %>% mutate(start_time = lubridate::parse_date_time(start_time, orders = c("HMS", "HM")),
                                            end_time = start_time + minutes(length))
    
    # Create a new column for each day_segment
    for(row_id in 1:nrow(day_segments)){
      row = day_segments[row_id,]
      sensor_data <- sensor_data %>% mutate(!!paste("local_day_segment", row_id, sep = "_") := ifelse(local_date_time_obj >= row$start_time & local_date_time_obj < row$end_time, 
                                                                                        paste0("[", 
                                                                                               row$label, "#", 
                                                                                               local_date, "#",
                                                                                               paste(str_pad(hour(row$start_time),2, pad="0"), str_pad(minute(row$start_time),2, pad="0"), str_pad(second(row$start_time),2, pad="0"),sep =":"), "#",
                                                                                               local_date, "#",
                                                                                               paste(str_pad(hour(row$end_time),2, pad="0"), str_pad(minute(row$end_time),2, pad="0"), str_pad(second(row$end_time),2, pad="0"),sep =":"),
                                                                                               "]"), NA))
    }
    # Join all day_segments in a single column
    sensor_data <- sensor_data %>% 
      unite("assigned_segments", starts_with("local_day_segment"), sep = "|", na.rm = TRUE) %>% 
      select(-local_date_time_obj)
 
    
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
             data = map2(data, inferred_day_segments, function(nested_data, segments){
               nested_data <- nested_data %>% mutate(assigned_segments = NA_character_, row_date_time = lubridate::ymd_hms(local_date_time, tz = local_timezone))
               for(row_id in 1:nrow(segments)){
                 row = segments[row_id,]
                 nested_data <- nested_data %>%
                   mutate(assigned_segments_temp = if_else(row_date_time %within% row$segment_interval, row$segment_id, NA_character_)) %>%
                   unite(col = "assigned_segments", c(assigned_segments, assigned_segments_temp), na.rm = TRUE, sep = "") %>%
                   mutate(assigned_segments = str_replace(assigned_segments, pattern = "\\]\\[", replacement = "\\]\\|\\[")) # this replaces ][ with ]|[
               }

               return(nested_data %>% select(-row_date_time))
             })
             ) %>% 
      unnest(cols = data) %>% 
      arrange(row_n) %>%
      select(-row_n, -existent_dates, -inferred_day_segments)
      
    
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
    
    
    sensor_data <- sensor_data %>%
      mutate(row_n = row_number()) %>% 
      group_by(local_timezone) %>% 
      nest() %>% 
      mutate(data = map(data, function(nested_data){
        nested_data <- nested_data %>% mutate(assigned_segments = NA_character_)
        for(row_id in 1:nrow(day_segments)){
          row = day_segments[row_id,]
          nested_data <- nested_data %>%
            mutate(assigned_segments_temp = if_else(timestamp >= row$segment_start & timestamp <= row$segment_end, row$segment_id, NA_character_)) %>%
            unite(col = "assigned_segments", c(assigned_segments, assigned_segments_temp), na.rm = TRUE, sep = "") %>% 
            mutate(assigned_segments = str_replace(assigned_segments, pattern = "\\]\\[", replacement = "\\]\\|\\[")) #replace ][ with ]|[
        }
        
        return(nested_data)
      })) %>% 
      unnest(cols = data) %>% 
      arrange(row_n) %>% 
      select(-row_n)
    

  }