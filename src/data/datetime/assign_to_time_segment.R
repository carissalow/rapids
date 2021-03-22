library("tidyverse")
library("lubridate", warn.conflicts = F)
options(scipen=999)

day_type_delay <- function(time_segments, day_type, include_past_periodic_segments){
  delay <- time_segments %>% mutate(length_duration = duration(length)) %>%  filter(repeats_on == day_type) %>% arrange(-length_duration) %>% pull(length_duration) %>% first()
  return(if_else(is.na(delay) | include_past_periodic_segments == FALSE, duration("0days"), delay))
}

get_segment_dates <- function(data, local_timezone, day_type, delay){
  dates <-  data %>% 
            distinct(local_date) %>% 
            mutate(local_date_obj = date(lubridate::ymd(local_date, tz = local_timezone))) %>% 
            complete(local_date_obj = seq(date(min(local_date_obj) - delay), date(max(local_date_obj) + delay), by="days")) %>%
            mutate(local_date = replace_na(as.character(date(local_date_obj))))

  if(day_type == "every_day")
    dates <- dates %>% mutate(every_day = 0)
  else if (day_type == "wday")
    dates <- dates %>% mutate(wday = wday(local_date_obj, week_start = 1))
  else if (day_type == "mday")
    dates <- dates %>% mutate(mday = mday(local_date_obj))
  else if (day_type == "qday")
    dates <- dates %>% mutate(qday = qday(local_date_obj))
  else if (day_type == "yday")
    dates <- dates %>% mutate(yday = yday(local_date_obj))
  return(dates)
}

create_nonoverlapping_periodic_segments <- function(nested_inferred_time_segments){
  new_segments <- (data.frame(nested_inferred_time_segments %>% 
    group_by(original_label) %>%
    mutate(max_groups = max(overlap_id) + 1) %>% 
    # select(label, segment_id_start, segment_id_end, overlap_id, max_groups) %>% 
    nest() %>%  
    mutate(data = map(data, function(nested_data){
        nested_data <- nested_data %>% arrange( segment_id_start, segment_id_end) %>% 
        group_by(segment_id_start) %>% 
        mutate(n_id = ((cur_group_id()-1) %% max_groups)) %>% 
        filter(overlap_id == n_id) %>% 
        # select(label, segment_id_start, overlap_id, n_id) %>% 
        ungroup()
    })) %>% 
    unnest(cols = data) %>% 
    ungroup()
  ))
  return(new_segments)
}


assign_rows_to_segments <- function(nested_data, nested_inferred_time_segments){
  nested_data <- nested_data %>% mutate(assigned_segments = "")
  for(i in seq_len(nrow(nested_inferred_time_segments))) {
    segment <- nested_inferred_time_segments[i,]
    nested_data$assigned_segments <- ifelse(segment$segment_start_ts<= nested_data$timestamp & segment$segment_end_ts >= nested_data$timestamp,
                                            stringi::stri_c(nested_data$assigned_segments, segment$segment_id, sep = "|"), nested_data$assigned_segments)
  }
  nested_data$assigned_segments <- substring(nested_data$assigned_segments, 2)
  return(nested_data)
}

assign_rows_to_segments_frequency <- function(nested_data, nested_timezone, time_segments){
  for(i in 1:nrow(time_segments)) {
    segment <- time_segments[i,]
    nested_data$assigned_segments <- ifelse(segment$segment_start_ts<= nested_data$local_time_obj & segment$segment_end_ts >= nested_data$local_time_obj,
                                            # The segment_id is assambled on the fly because it depends on each row's local_date and timezone 
                                            stringi::stri_c("[",
                                                            segment[["label"]], "#",
                                                            nested_data$local_date, " ",
                                                            segment[["segment_id_start_time"]], ",",
                                                            nested_data$local_date, " ",
                                                            segment[["segment_id_end_time"]], ";",
                                                            as.numeric(lubridate::as_datetime(stringi::stri_c(nested_data$local_date, segment$segment_id_start_time), tz = nested_timezone)) * 1000, ",",
                                                            as.numeric(lubridate::as_datetime(stringi::stri_c(nested_data$local_date, segment$segment_id_end_time), tz = nested_timezone)) * 1000 + 999,
                                                            "]"),
                                            nested_data$assigned_segments)
  }
  return(nested_data)
}

assign_to_time_segment <- function(sensor_data, time_segments, time_segments_type, include_past_periodic_segments){
  
  if(nrow(sensor_data) == 0 || nrow(time_segments) == 0)
    return(sensor_data %>% mutate(assigned_segments = NA))
  
  if(time_segments_type == "FREQUENCY"){
    
    time_segments <- time_segments %>% mutate(start_time = lubridate::hm(start_time),
                                            end_time = start_time + minutes(length) - seconds(1),
                                            segment_id_start_time = paste(str_pad(hour(start_time),2, pad="0"), str_pad(minute(start_time),2, pad="0"), str_pad(second(start_time),2, pad="0"),sep =":"),
                                            segment_id_end_time = paste(str_pad(hour(ymd("1970-01-01") + end_time),2, pad="0"), str_pad(minute(ymd("1970-01-01") + end_time),2, pad="0"), str_pad(second(ymd("1970-01-01") + end_time),2, pad="0"),sep =":"), # add ymd("1970-01-01") to get a real time instead of duration
                                            segment_start_ts = as.numeric(start_time),
                                            segment_end_ts = as.numeric(end_time))
    
    sensor_data <- sensor_data %>% mutate(local_time_obj = as.numeric(lubridate::hms(local_time)),
                                          assigned_segments = "")
    
    sensor_data <- sensor_data %>%
      group_by(local_timezone) %>% 
      nest() %>% 
      mutate(data = map2(data, local_timezone, assign_rows_to_segments_frequency, time_segments)) %>% 
      unnest(cols = data) %>% 
      arrange(timestamp) %>% 
      select(-local_time_obj) %>% 
      ungroup()
    
    return(sensor_data)

    
  } else if (time_segments_type == "PERIODIC"){
    
    # We need to take into account segment start dates that could include the first day of data
    time_segments <- time_segments %>% mutate(length_duration = duration(length))
    every_day_delay <- duration("0days")
    wday_delay <- day_type_delay(time_segments, "wday", include_past_periodic_segments)
    mday_delay <- day_type_delay(time_segments, "mday", include_past_periodic_segments)
    qday_delay <- day_type_delay(time_segments, "qday", include_past_periodic_segments)
    yday_delay <- day_type_delay(time_segments, "yday", include_past_periodic_segments)
    
    sensor_data <- sensor_data %>%
      group_by(local_timezone) %>% 
      nest() %>% 
      # get existent days that we need to start segments from
      mutate(every_date = map2(data, local_timezone, get_segment_dates, "every_day", every_day_delay),
             week_dates = map2(data, local_timezone, get_segment_dates, "wday", wday_delay),
             month_dates = map2(data, local_timezone, get_segment_dates, "mday", mday_delay),
             quarter_dates = map2(data, local_timezone, get_segment_dates, "qday", qday_delay),
             year_dates = map2(data, local_timezone, get_segment_dates, "yday", yday_delay),
             existent_dates = pmap(list(every_date, week_dates, month_dates, quarter_dates, year_dates),
                                   function(every_date, week_dates, month_dates, quarter_dates, year_dates) reduce(list(every_date, week_dates,month_dates, quarter_dates, year_dates), .f=full_join)),
             # build the actual time segments taking into account the users requested length and repeat schedule
             inferred_time_segments = map(existent_dates,
                                         ~ crossing(time_segments, .x) %>%
                                           pivot_longer(cols = c(every_day,wday, mday, qday, yday), names_to = "day_type", values_to = "day_value") %>%
                                           filter(repeats_on == day_type & repeats_value == day_value) %>%
                                           # The segment ids (segment_id_start and segment_id_end) are computed in UTC to avoid having different labels for instances of a segment that happen in different timezones
                                           mutate(segment_id_start = lubridate::parse_date_time(paste(local_date, start_time), orders = c("Ymd HMS", "Ymd HM")) + period(overlap_duration),
                                                  segment_id_end = segment_id_start + lubridate::duration(length),
                                                  # The actual segments are computed using timestamps taking into account the timezone
                                                  segment_start_ts = as.numeric(lubridate::parse_date_time(paste(local_date, start_time), orders = c("Ymd HMS", "Ymd HM"), tz = local_timezone) + period(overlap_duration)) * 1000,
                                                  segment_end_ts = segment_start_ts + as.numeric(lubridate::duration(length)) * 1000 + 999,
                                                  segment_id = paste0("[",
                                                                      paste0(label,"#",
                                                                             paste0(lubridate::date(segment_id_start), " ",
                                                                                    paste(str_pad(hour(segment_id_start),2, pad="0"), str_pad(minute(segment_id_start),2, pad="0"), str_pad(second(segment_id_start),2, pad="0"),sep =":"), ",",
                                                                                    lubridate::date(segment_id_end), " ",
                                                                                    paste(str_pad(hour(segment_id_end),2, pad="0"), str_pad(minute(segment_id_end),2, pad="0"), str_pad(second(segment_id_end),2, pad="0"),sep =":")),";",
                                                                             paste0(segment_start_ts, ",", segment_end_ts)),
                                                                      "]")) %>% 
                                           # drop time segments with an invalid start or end time (mostly due to daylight saving changes, e.g. 2020-03-08 02:00:00 EST does not exist, clock jumps from 01:59am to 03:00am)
                                           drop_na(segment_start_ts, segment_end_ts)), 
              inferred_time_segments = map(inferred_time_segments, create_nonoverlapping_periodic_segments),
             data = map2(data, inferred_time_segments, assign_rows_to_segments)
      ) %>%
      select(-existent_dates, -inferred_time_segments, -every_date, -week_dates, -month_dates, -quarter_dates, -year_dates) %>%
      unnest(cols = data) %>% 
      arrange(timestamp)

  } else if ( time_segments_type == "EVENT"){

    sensor_data <- sensor_data %>% 
      group_by(local_timezone) %>% 
      nest() %>% 
      mutate(inferred_time_segments = map(local_timezone, function(tz){ 
                                         inferred <- time_segments %>% 
                                           mutate(shift = ifelse(shift == "0", "0seconds", shift),
                                                  segment_start_ts = event_timestamp + (as.integer(seconds(lubridate::duration(shift))) * ifelse(shift_direction >= 0, 1, -1) * 1000),
                                                  segment_end_ts = segment_start_ts + (as.integer(seconds(lubridate::duration(length))) * 1000),
                                                  # these start and end datetime objects are for labeling only
                                                  segment_id_start = lubridate::as_datetime(segment_start_ts/1000, tz = tz), 
                                                  segment_id_end = lubridate::as_datetime(segment_end_ts/1000, tz = tz),
                                                  segment_end_ts = segment_end_ts + 999,
                                                  segment_id = paste0("[",
                                                                      paste0(label,"#",
                                                                             paste0(lubridate::date(segment_id_start), " ",
                                                                                    paste(str_pad(hour(segment_id_start),2, pad="0"), str_pad(minute(segment_id_start),2, pad="0"), str_pad(second(segment_id_start),2, pad="0"),sep =":"), ",",
                                                                                    lubridate::date(segment_id_end), " ",
                                                                                    paste(str_pad(hour(segment_id_end),2, pad="0"), str_pad(minute(segment_id_end),2, pad="0"), str_pad(second(segment_id_end),2, pad="0"),sep =":")),";",
                                                                             paste0(segment_start_ts, ",", segment_end_ts)),
                                                                      "]"))
                                         # Check that for overlapping segments (not allowed because our resampling episode algorithm would have to have a second instead of minute granularity that increases storage and computation time)
                                         overlapping <-  inferred  %>% group_by(label) %>% arrange(segment_start_ts) %>% 
                                           mutate(overlaps = if_else(segment_start_ts <= lag(segment_end_ts), TRUE, FALSE),
                                                  overlapping_segments = paste(paste(lag(label), lag(event_timestamp), lag(length), lag(shift), lag(shift_direction), lag(device_id), sep = ","),"and",
                                                                                paste(label, event_timestamp, length, shift, shift_direction, device_id, sep = ",")))
                                         if(any(overlapping$overlaps, na.rm = TRUE)){
                                           stop(paste0("\n\nOne or more event time segments overlap for ",overlapping$device_id[[1]],", modify their lengths so they don't:\n", paste0(overlapping %>% filter(overlaps == TRUE) %>% pull(overlapping_segments), collapse = "\n"), "\n\n"))
                                         } else{
                                           return(inferred)
                                         }}),
             data = map2(data, inferred_time_segments, assign_rows_to_segments)) %>% 
      select(-inferred_time_segments) %>% 
      unnest(data) %>% 
      arrange(timestamp)
  }

  return(sensor_data %>% ungroup())
}