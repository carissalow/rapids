library("tidyverse")
library("lubridate")

find_segments_frequency <- function(local_date, local_time, segments){
  assigned_segments <- segments[segments$segment_start<= local_time & segments$segment_end >= local_time, ]

  return(stringi::stri_c(stringi::stri_c("[", 
                                assigned_segments[["label"]], "#",
                                local_date, "#",
                                assigned_segments[["segment_id_start_time"]], "#",
                                local_date, "#",
                                assigned_segments[["segment_id_end_time"]], 
                                "]"), collapse = "|"))
}

find_segments_periodic <- function(timestamp, segments){
  return(stringi::stri_c(segments[[1]][segments[[1]]$segment_start_ts<= timestamp & segments[[1]]$segment_end_ts >= timestamp, "segment_id"][["segment_id"]], collapse = "|"))
}

  # We might need to optimise the event function as well, filter, and pull are slow
find_segments_event <- function(timestamp, segments){
  return(stringi::stri_c(segments %>% 
                           filter(segment_start <= timestamp  & segment_end >= timestamp) %>% 
                           pull(segment_id), collapse = "|"))
}

assign_to_day_segment <- function(sensor_data, day_segments, day_segments_type, include_past_periodic_segments){
  if(nrow(sensor_data) == 0)
    return(sensor_data %>% mutate(assigned_segments = NA))

  if(day_segments_type == "FREQUENCY"){ #FREQUENCY
    
    day_segments <- day_segments %>% mutate(start_time = lubridate::hm(start_time),
                                            end_time = start_time + minutes(length) - seconds(1),
                                            segment_id_start_time = paste(str_pad(hour(start_time),2, pad="0"), str_pad(minute(start_time),2, pad="0"), str_pad(second(start_time),2, pad="0"),sep =":"),
                                            segment_id_end_time = paste(str_pad(hour(ymd("1970-01-01") + end_time),2, pad="0"), str_pad(minute(ymd("1970-01-01") + end_time),2, pad="0"), str_pad(second(ymd("1970-01-01") + end_time),2, pad="0"),sep =":"), # add ymd("1970-01-01") to get a real time instead of duration
                                            segment_start = as.numeric(start_time),
                                            segment_end = as.numeric(end_time))
    sensor_data <- sensor_data %>% mutate(local_time_obj = as.numeric(lubridate::hms(local_time)),
                                          assigned_segments = map2_chr(local_date, local_time_obj, ~find_segments_frequency(.x, .y, day_segments))) %>% select(-local_time_obj)
    
  } else if (day_segments_type == "PERIODIC"){ #PERIODIC
    
    # We need to take into account segment start dates that could include the first day of data
    day_segments <- day_segments %>% mutate(length_duration = duration(length))
    wday_delay <- day_segments %>% mutate(length_duration = duration(length)) %>%  filter(repeats_on == "wday") %>% arrange(-length_duration) %>% pull(length_duration) %>% first()
    wday_delay <- if_else(is.na(wday_delay) | include_past_periodic_segments == FALSE, duration("0days"), wday_delay)
    
    mday_delay <- day_segments %>% mutate(length_duration = duration(length)) %>%  filter(repeats_on == "mday") %>% arrange(-length_duration) %>% pull(length_duration) %>% first()
    mday_delay <- if_else(is.na(mday_delay) | include_past_periodic_segments == FALSE, duration("0days"), mday_delay)
    
    qday_delay <- day_segments %>% mutate(length_duration = duration(length)) %>%  filter(repeats_on == "qday") %>% arrange(-length_duration) %>% pull(length_duration) %>% first()
    qday_delay <- if_else(is.na(qday_delay) | include_past_periodic_segments == FALSE, duration("0days"), qday_delay)
    
    yday_delay <- day_segments %>% mutate(length_duration = duration(length)) %>%  filter(repeats_on == "yday") %>% arrange(-length_duration) %>% pull(length_duration) %>% first()
    yday_delay <- if_else(is.na(yday_delay) | include_past_periodic_segments == FALSE, duration("0days"), yday_delay)
    
    sensor_data <- sensor_data %>%
      # mutate(row_n = row_number()) %>% 
      group_by(local_timezone) %>% 
      nest() %>% 
      # get existent days that we need to start segments from
      mutate(every_date = map(data, ~.x %>% 
                                    distinct(local_date) %>% 
                                    mutate(local_date_obj = date(lubridate::ymd(local_date, tz = local_timezone))) %>% 
                                    complete(local_date_obj = seq(min(local_date_obj), max(local_date_obj), by="days")) %>%
                                    mutate(local_date = replace_na(as.character(date(local_date_obj)))) %>% 
                                    mutate(every_day = 0)),
             week_dates = map(data, ~.x %>% 
                               distinct(local_date) %>% 
                                mutate(local_date_obj = date(lubridate::ymd(local_date, tz = local_timezone))) %>% 
                                complete(local_date_obj = seq(date(min(local_date_obj) - wday_delay), max(local_date_obj), by="days")) %>%
                               mutate(local_date = replace_na(as.character(date(local_date_obj)))) %>% 
                               mutate(wday = wday(local_date_obj, week_start = 1))  ), 
             month_dates = map(data, ~.x %>% 
                                   distinct(local_date) %>% 
                                   mutate(local_date_obj = date(lubridate::ymd(local_date, tz = local_timezone))) %>% 
                                   complete(local_date_obj = seq(date(min(local_date_obj) - mday_delay), max(local_date_obj), by="days")) %>%
                                   mutate(local_date = replace_na(as.character(date(local_date_obj)))) %>%
                                   mutate(mday = mday(local_date_obj))), 
             quarter_dates = map(data, ~.x %>% 
                                distinct(local_date) %>% 
                                mutate(local_date_obj = date(lubridate::ymd(local_date, tz = local_timezone))) %>% 
                                complete(local_date_obj = seq(date(min(local_date_obj) - qday_delay), max(local_date_obj), by="days")) %>%
                                mutate(local_date = replace_na(as.character(date(local_date_obj)))) %>%
                                mutate(qday = qday(local_date_obj)) ), 
             year_dates = map(data, ~.x %>% 
                                 distinct(local_date) %>% 
                                 mutate(local_date_obj = date(lubridate::ymd(local_date, tz = local_timezone))) %>% 
                                 complete(local_date_obj = seq(date(min(local_date_obj) - yday_delay), max(local_date_obj), by="days")) %>%
                                 mutate(local_date = replace_na(as.character(date(local_date_obj)))) %>%
                                 mutate(yday = yday(local_date_obj)) ),
             existent_dates = pmap(list(every_date, week_dates, month_dates, quarter_dates, year_dates),
                                    function(every_date, week_dates, month_dates, quarter_dates, year_dates) reduce(list(every_date, week_dates,month_dates, quarter_dates, year_dates), .f=full_join)),
             every_date = NULL,
             week_dates = NULL,
             month_dates = NULL,
             quarter_dates = NULL,
             year_dates = NULL,
             # build the actual day segments taking into account the users requested leangth and repeat schedule
             inferred_day_segments = map(existent_dates,
                                         ~ crossing(day_segments, .x) %>%
                                           pivot_longer(cols = c(every_day,wday, mday, qday, yday), names_to = "day_type", values_to = "day_value") %>%
                                           filter(repeats_on == day_type & repeats_value == day_value) %>%
                                           mutate(segment_id_start = lubridate::parse_date_time(paste(local_date, start_time), orders = c("Ymd HMS", "Ymd HM")), # The segment ids (label#start#end) are computed in UTC to avoid having different labels for instances of a segment that happen in different timezones
                                                  segment_id_end = segment_id_start + lubridate::duration(length),
                                                  segment_start_ts = as.numeric(lubridate::parse_date_time(paste(local_date, start_time), orders = c("Ymd HMS", "Ymd HM"), tz = local_timezone)), # The actual segments are computed using timestamps taking into account the timezone
                                                  segment_end_ts = segment_start_ts + as.numeric(lubridate::duration(length)),
                                                  segment_id = paste0("[",
                                                                      paste(sep= "#",
                                                                            label,
                                                                            lubridate::date(segment_id_start),
                                                                            paste(str_pad(hour(segment_id_start),2, pad="0"),
                                                                                  str_pad(minute(segment_id_start),2, pad="0"),
                                                                                  str_pad(second(segment_id_start),2, pad="0"),sep =":"),
                                                                            lubridate::date(segment_id_end),
                                                                            paste(str_pad(hour(segment_id_end),2, pad="0"),
                                                                                  str_pad(minute(segment_id_end),2, pad="0"),
                                                                                  str_pad(second(segment_id_end),2, pad="0"),sep =":")
                                                                      ),
                                                                      "]")) %>%
                                           select(segment_start_ts, segment_end_ts,  segment_id)),
             # loop thorugh every day segment and assigned it to the rows that fall within its start and end
             data = map2(data, inferred_day_segments, ~ .x %>% mutate(row_date_time = as.numeric(lubridate::ymd_hms(local_date_time, tz = local_timezone)),
                                                                      assigned_segments = map_chr(row_date_time, ~find_segments_periodic(.x, inferred_day_segments)),
                                                                      row_date_time = NULL))
      ) %>%
      select(-existent_dates, -inferred_day_segments) %>%
      unnest(cols = data) %>% 
      arrange(timestamp)
    
    
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