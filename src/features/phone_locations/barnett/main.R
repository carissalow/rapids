source("renv/activate.R")
library("dplyr", warn.conflicts = F)
library("stringr")
library("lubridate")
library("purrr")

create_empty_file <- function(requested_features){
  return(data.frame(local_segment= character(), 
                    hometime= numeric(), 
                    disttravelled= numeric(), 
                    rog= numeric(), 
                    maxdiam= numeric(), 
                    maxhomedist= numeric(), 
                    siglocsvisited= numeric(), 
                    avgflightlen= numeric(), 
                    stdflightlen= numeric(), 
                    avgflightdur= numeric(), 
                    stdflightdur= numeric(), 
                    probpause= numeric(), 
                    siglocentropy= numeric(), 
                    minsmissing= numeric(), 
                    circdnrtn= numeric(), 
                    wkenddayrtn= numeric(),
                    minutes_data_used= numeric()
  ) %>% select(all_of(requested_features)))
}

summarise_multiday_segments <- function(segments, features){
  features <- features %>% mutate(local_date=ymd(local_date))
  segments <- segments %>% extract(col = local_segment, 
                                   into = c ("local_segment_start_datetime", "local_segment_end_datetime"),
                                   ".*#(.*) .*,(.*) .*", 
                                   remove = FALSE) %>% 
    mutate(local_segment_start_datetime = ymd(local_segment_start_datetime),
           local_segment_end_datetime = ymd(local_segment_end_datetime)) %>% 
    group_by(local_segment) %>% 
    nest() %>% 
    mutate(data = map(data, function(nested_data, nested_features){
      
      summary <- nested_features %>% filter(local_date >= nested_data$local_segment_start_datetime &
                                          local_date <= nested_data$local_segment_end_datetime)
      if(nrow(summary) > 0)
        summary <-  summary %>% 
                 summarise(across(c(hometime, disttravelled, siglocsvisited, minutes_data_used), sum),
                           across(c(maxdiam, maxhomedist), max),
                           across(c(rog, avgflightlen, stdflightlen, avgflightdur, stdflightdur, probpause, siglocentropy, circdnrtn, wkenddayrtn, minsmissing), mean))
      return(summary)
      
    }, features)) %>% 
    unnest(cols = everything()) %>%
    ungroup()
  return(segments)
}

barnett_features <- function(sensor_data_files, time_segment, params){
  location_features <- NULL
  daily_features <- read.csv(sensor_data_files[["barnett_daily"]], stringsAsFactors = FALSE)
  location <- read.csv(sensor_data_files[["sensor_data"]], stringsAsFactors = FALSE)
  minutes_data_used <- params[["MINUTES_DATA_USED"]]
  available_features <- c("hometime","disttravelled","rog","maxdiam", "maxhomedist","siglocsvisited","avgflightlen", "stdflightlen",
                          "avgflightdur","stdflightdur", "probpause","siglocentropy", "circdnrtn","wkenddayrtn")
  requested_features <- intersect(unlist(params["FEATURES"], use.names = F), available_features)
  requested_features <- c("local_segment", requested_features)
  if(minutes_data_used)
    requested_features <- c(requested_features, "minutes_data_used")
  
  if (nrow(location) > 0 & nrow(daily_features) > 0){
    location <- location %>% filter_data_by_segment(time_segment)
    datetime_start_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2} 00:00:00"
    datetime_end_regex = "[0-9]{4}[\\-|\\/][0-9]{2}[\\-|\\/][0-9]{2} 23:59:59"
    location <- location %>% mutate(is_daily = str_detect(local_segment, paste0(time_segment, "#", datetime_start_regex, ",", datetime_end_regex))) 
    if(nrow(location) == 0 || !all(location$is_daily)){
      message(paste("Barnett's location features cannot be computed for data or time segmentes that do not span entire days (00:00:00 to 23:59:59). Skipping ", time_segment))
      location_features <- create_empty_file(requested_features)  
    } else {
      location_dates_segments <- location %>% select(local_segment) %>% distinct(local_segment, .keep_all = TRUE)
      features <- summarise_multiday_segments(location_dates_segments, daily_features)
      location_features <- features %>% select(all_of(requested_features))
    } 
  } else
    location_features <- create_empty_file(requested_features)  
  return(location_features)
}