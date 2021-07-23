get_existent_dates <- function(data, time_segments, include_past_periodic_segments){
  max_delay = max(time_segments$length_duration)
  max_delay <- (if_else(is.na(max_delay) | include_past_periodic_segments == FALSE, duration("0days"), max_delay))
  
  existent_dates <- data %>% 
    distinct(local_date, .keep_all = FALSE) %>% 
    mutate(local_date_obj = lubridate::ymd(local_date)) %>% 
    complete(local_date_obj = seq(date(min(local_date_obj) - max_delay), max(local_date_obj), by="days")) %>%
    mutate(local_date = replace_na(as.character(local_date_obj)),
           every_day = 0,
           wday = wday(local_date_obj, week_start = 1),
           mday = mday(local_date_obj),
           qday = qday(local_date_obj),
           yday = yday(local_date_obj)) %>% 
    select(-local_date_obj)
}

infer_existent_periodic_segments <- function(existent_dates, segments){
  # build the actual time segments taking into account the data and users' requested length and repeat schedule
  # segment datetime labels are computed on UTC
  crossing(segments, existent_dates) %>%
    pivot_longer(cols = c(every_day,wday, mday, qday, yday), names_to = "day_type", values_to = "day_value") %>%
    filter(repeats_on == day_type & repeats_value == day_value) %>%
    mutate(segment_id_start = lubridate::parse_date_time(paste(local_date, start_time), orders = c("Ymd HMS", "Ymd HM")) + period(overlap_duration),
            segment_id_end = segment_id_start + lubridate::period(length)) %>% 
    select(original_label, label, segment_id_start, segment_id_end, overlap_id, length)
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
  if(nrow(nested_inferred_time_segments) == 0)
    return(nested_inferred_time_segments)
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



add_periodic_segment_timestamps_and_id <- function(data, segments, local_timezone){
  # segment timestamps are computed on the data's timezone(s)
  time_format_fn <- stamp("23:51:15", orders="HMS", quiet = TRUE)
  segments %>% mutate(segment_id_start_tz = lubridate::force_tz(segment_id_start, tzone = local_timezone),
                      segment_start_ts = as.numeric(segment_id_start_tz) * 1000,
                      segment_end_ts = as.numeric(segment_id_start_tz + lubridate::period(length)) * 1000 + 999,
                      segment_id_start_tz = NULL,
                      segment_id = glue("[{label}#{start_date} {start_time},{end_date} {end_time};{segment_start_ts},{segment_end_ts}]",
                                        start_date=lubridate::date(segment_id_start),
                                        start_time=time_format_fn(segment_id_start),
                                        end_date=lubridate::date(segment_id_end),
                                        end_time=time_format_fn(segment_id_end) )) %>% 
                drop_na(segment_start_ts, segment_end_ts)
}

assign_to_periodic_segments <- function(sensor_data, time_segments, include_past_periodic_segments){
  time_segments <- time_segments %>% mutate(length_duration = duration(length))
  existent_dates <- get_existent_dates(sensor_data, time_segments, include_past_periodic_segments)
  inferred_segments <- infer_existent_periodic_segments(existent_dates, time_segments) %>%
    dedup_nonoverlapping_periodic_segments()

  sensor_data <- sensor_data %>%
    group_by(local_timezone) %>%
    nest() %>%
    mutate(localised_time_segments = map(data, add_periodic_segment_timestamps_and_id, inferred_segments, local_timezone),
          data = map2(data, localised_time_segments, assign_rows_to_segments)) %>%
    select(-localised_time_segments) %>%
    unnest(cols = data) %>%
    arrange(timestamp) %>%
    ungroup()
  
  return(sensor_data)
}