validate_overlapping_event_segments <- function(segments){
  # Check for overlapping segments (not allowed because our resampling episode algorithm would have to have a second instead of minute granularity that increases storage and computation time)
  overlapping <-  segments  %>% 
    group_by(label) %>% 
    arrange(segment_start_ts) %>% 
    mutate(overlaps = if_else(segment_start_ts <= lag(segment_end_ts), TRUE, FALSE),
           overlapping_segments = glue("a) [{lag(label)},\t{lag(event_timestamp)},\t{lag(length)},\t{lag(shift)},\t{lag(shift_direction)},\t{lag(device_id)}] \n",
                                       "b) [{label},\t{event_timestamp},\t{length},\t{shift},\t{shift_direction},\t{device_id}]"))
  if(any(overlapping$overlaps, na.rm = TRUE))
    stop("One or more event time segments overlap for ",overlapping$device_id[[1]],
         ", modify their lengths so they don't:\n", paste0(overlapping %>% filter(overlaps == TRUE) %>% pull(overlapping_segments), collapse = "\n"))
}

infer_event_segments <- function(tz, segments){ 
  time_format_fn <- stamp("23:51:15", orders="HMS", quiet = TRUE)
  inferred <- segments %>% 
    mutate(shift = ifelse(shift == "0", "0seconds", shift),
           segment_start_ts = event_timestamp + (as.integer(seconds(lubridate::duration(shift))) * ifelse(shift_direction >= 0, 1, -1) * 1000),
           segment_end_ts = segment_start_ts + (as.integer(seconds(lubridate::duration(length))) * 1000),
           segment_id_start = lubridate::as_datetime(segment_start_ts/1000, tz = tz), 
           segment_id_end = lubridate::as_datetime(segment_end_ts/1000, tz = tz),
           segment_end_ts = segment_end_ts + 999,
           segment_id = glue("[{label}#{start_date} {start_time},{end_date} {end_time};{segment_start_ts},{segment_end_ts}]",
                             start_date=lubridate::date(segment_id_start),
                             start_time=time_format_fn(segment_id_start),
                             end_date=lubridate::date(segment_id_end),
                             end_time=time_format_fn(segment_id_end)))
  validate_overlapping_event_segments(inferred)
  return(inferred)
}

assign_to_event_segments <- function(sensor_data, time_segments){
  sensor_data <- sensor_data %>% 
    group_by(local_timezone) %>% 
    nest() %>% 
    mutate(inferred_time_segments = map(local_timezone, infer_event_segments, time_segments),
           data = map2(data, inferred_time_segments, assign_rows_to_segments)) %>% 
    select(-inferred_time_segments) %>% 
    unnest(data) %>% 
    arrange(timestamp) %>%
    ungroup()
}