day_type_delay <- function(time_segments, day_type, include_past_periodic_segments){
  # Return a delay in days to consider or not the first row of data
  delay <- time_segments %>% 
    mutate(length_duration = duration(length)) %>%  
    filter(repeats_on == day_type) %>% arrange(-length_duration) %>% 
    pull(length_duration) %>% 
    first()
  return(if_else(is.na(delay) | include_past_periodic_segments == FALSE, duration("0days"), delay))
}

get_segment_dates <- function(data, local_timezone, day_type, delay){
  # Based on the data we are processing we extract unique dates to build segments
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

infer_existent_periodic_segments <- function(existent_dates, segments){
  # build the actual time segments taking into account the data and users' requested length and repeat schedule
  # segment datetime labels are computed on UTC
  crossing(segments, existent_dates) %>%
    pivot_longer(cols = c(every_day,wday, mday, qday, yday), names_to = "day_type", values_to = "day_value") %>%
    filter(repeats_on == day_type & repeats_value == day_value) %>%
    mutate(segment_id_start = lubridate::parse_date_time(paste(local_date, start_time), orders = c("Ymd HMS", "Ymd HM")) + period(overlap_duration),
            segment_id_end = segment_id_start + lubridate::duration(length))
}

dedup_nonoverlapping_periodic_segments <- function(nested_inferred_time_segments){
  # Overlapping segments exist when their length is longer than their repeating frequency, e.g. twoday segements starting on every day
  # In process_time_segments we decompose those segments into non-overlapping ones, e.g. twodayA +0days and twodayB +1days
  # This means that any date will have more than one non-overlapping instances, that we need to dedup
  # We choose alternating non-overlapping instances to guarantee any data row is only neeeded in one instance at a time
  # d1,r1,twoday0
  # d2,r2,twoday0 twoday1 
  # d3,r3,twoday1 twoday0 
  # d4,r4,twoday0 twoday1 
  new_segments <- data.frame(nested_inferred_time_segments %>% 
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
                                ungroup())
}



add_periodic_segment_timestamps_and_id <- function(segments, local_timezone){
  # segment timestamps are computed on the data's timezone(s)
  time_format_fn <- stamp("23:51:15", orders="HMS", quiet = TRUE)
  segments %>% mutate(segment_start_ts = as.numeric(lubridate::force_tz(segment_id_start, tzone = local_timezone)) * 1000,
                      segment_end_ts = segment_start_ts + as.numeric(lubridate::duration(length)) * 1000 + 999,
                      segment_id = glue("[{label}#{start_date} {start_time},{end_date} {end_time};{segment_start_ts},{segment_end_ts}]",
                                        start_date=lubridate::date(segment_id_start),
                                        start_time=time_format_fn(segment_id_start),
                                        end_date=lubridate::date(segment_id_end),
                                        end_time=time_format_fn(segment_id_end) )) %>% 
                drop_na(segment_start_ts, segment_end_ts)
}

assign_to_periodic_segments <- function(sensor_data, time_segments, include_past_periodic_segments){
  time_segments <- time_segments %>% mutate(length_duration = duration(length))
  every_day_delay <- duration("0days")
  wday_delay <- day_type_delay(time_segments, "wday", include_past_periodic_segments)
  mday_delay <- day_type_delay(time_segments, "mday", include_past_periodic_segments)
  qday_delay <- day_type_delay(time_segments, "qday", include_past_periodic_segments)
  yday_delay <- day_type_delay(time_segments, "yday", include_past_periodic_segments)
  
  sensor_data <- sensor_data %>%
    group_by(local_timezone) %>% 
    nest() %>% 
    mutate(every_date = map2(data, local_timezone, get_segment_dates, "every_day", every_day_delay),
           week_dates = map2(data, local_timezone, get_segment_dates, "wday", wday_delay),
           month_dates = map2(data, local_timezone, get_segment_dates, "mday", mday_delay),
           quarter_dates = map2(data, local_timezone, get_segment_dates, "qday", qday_delay),
           year_dates = map2(data, local_timezone, get_segment_dates, "yday", yday_delay),
           existent_dates = pmap(list(every_date, week_dates, month_dates, quarter_dates, year_dates), function(every_date, week_dates, month_dates, quarter_dates, year_dates) reduce(list(every_date, week_dates,month_dates, quarter_dates, year_dates), .f=full_join)),
           inferred_time_segments = map(existent_dates, infer_existent_periodic_segments, time_segments), 
           inferred_time_segments = map(inferred_time_segments, dedup_nonoverlapping_periodic_segments),
           inferred_time_segments = map(inferred_time_segments, add_periodic_segment_timestamps_and_id, local_timezone),
           data = map2(data, inferred_time_segments, assign_rows_to_segments)) %>%
    select(-existent_dates, -inferred_time_segments, -every_date, -week_dates, -month_dates, -quarter_dates, -year_dates) %>%
    unnest(cols = data) %>% 
    arrange(timestamp) %>% 
    ungroup()
  
  return(sensor_data)
}