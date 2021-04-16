source("renv/activate.R")
library(tidyverse)
library(lubridate)
library(glue)

create_empty_dataframe <- function(episode_type){
  integer_columns <- c("countepisode{episode_type}", "starttimefirstepisode{episode_type}", "endtimefirstepisode{episode_type}", "starttimelastepisode{episode_type}", "endtimelastepisode{episode_type}", "starttimelongestepisode{episode_type}", "endtimelongestepisode{episode_type}")
  integer_columns <- sapply(integer_columns, function(x) glue(x), simplify = TRUE, USE.NAMES = FALSE)
  double_columns <- c()
  for(col in c("duration", "calories", "mets"))
    for(fun in c("sum", "mean", "min","max","sd"))
      double_columns <- c(double_columns, glue("{fun}{col}episode{episode_type}"))
  
  as_tibble(c(sapply(integer_columns, function(x) integer()), sapply(double_columns, function(x) numeric())))
}

longest <- function(duration, time){
  position_longest <- min(which(duration == max(duration)))
  time[position_longest]
}

episode_type_features <- function(data, episode_type, episode_id_column){
  if(nrow(data) == 0)
    return(create_empty_dataframe(episode_type))
  
  data %>%
    group_by(across(all_of(episode_id_column))) %>% 
    summarise(duration = (max(timestamp) - min(timestamp)) / 60000 + 1,
              mets = sum(mets),
              calories = sum(value),
              start_time = min(time_since_ref),
              end_time = max(time_since_ref) + 1) %>% 
    summarise("countepisode{episode_type}" := n(), 
              "starttimefirstepisode{episode_type}" := first(start_time),
              "endtimefirstepisode{episode_type}" := first(end_time),
              "starttimelastepisode{episode_type}" := last(start_time),
              "endtimelastepisode{episode_type}" := last(end_time),
              "starttimelongestepisode{episode_type}" := longest(duration, start_time),
              "endtimelongestepisode{episode_type}" := longest(duration, end_time),
              across(duration, list(sum=sum, avg=mean, min=min,max=max,std=sd), .names = "{.fn}{.col}episode{episode_type}"),
              across(calories, list(sum=sum, avg=mean, min=min,max=max,std=sd), .names = "{.fn}{.col}episode{episode_type}"),
              across(mets, list(sum=sum, avg=mean, min=min,max=max,std=sd), .names = "{.fn}{.col}episode{episode_type}"))
}

rapids_features <- function(sensor_data_files, time_segment, provider){
    calories <- read_csv(snakemake@input[["sensor_data"]], 
                          col_types = cols_only(level="i", mets="d", value="d", local_date_time="T",assigned_segments="c", timestamp="d"))# %>%
    MET_THRESHOLD <- provider[["EPISODE_MET_THRESHOLD"]]
    MVPA_LABELS <- provider[["EPISODE_MVPA_CATEGORIES"]]
    FITBIT_LEVELS <- c("sedentary", "lightlyactive", "fairlyactive", "veryactive")
    MVPA_LEVELS <- which(FITBIT_LEVELS %in% MVPA_LABELS) - 1
    EPISODE_TIME_THRESHOLD <- provider[["EPISODE_TIME_THRESHOLD"]]
    EPISODE_REFERENCE_TIME <- provider[["EPISODE_REFERENCE_TIME"]]
    REQUESTED_EPISODES <-  provider[["EPISODE_TYPE"]]
    REQUESTED_FEATURES <- provider[["FEATURES"]]

    calories <- calories %>% filter_data_by_segment(time_segment) 

    if(nrow(calories) == 0)
      return(bind_cols(lapply(REQUESTED_EPISODES, function(episode_type) episode_type_features(calories, episode_type, ""))) %>% 
              add_column(local_segment = character(), .before = 1) %>%
              select(starts_with(c("local_segment", REQUESTED_FEATURES))))

    calories <- calories %>% 
      extract(timestamps_segment, regex = "(\\d*),", into = c("segment_start_ts"), remove = TRUE, convert = TRUE) %>%
      arrange(timestamp) %>% 
      mutate(consecutive = c(0,diff(timestamp) / 60000),
            level_diff = c(0, diff(level)),
            mvpa_diff = c(1, diff(if_else(level %in% MVPA_LEVELS, 1, 0))),
            met_diff = c(1, diff(if_else(mets >= MET_THRESHOLD, 1, 0))),
            level_episode_id = cumsum(consecutive > EPISODE_TIME_THRESHOLD | level_diff != 0),
            mvpa_episode_id = cumsum(consecutive > EPISODE_TIME_THRESHOLD | mvpa_diff != 0),
            met_episode_id = cumsum(consecutive > EPISODE_TIME_THRESHOLD | met_diff != 0),
            time_since_ref = case_when(EPISODE_REFERENCE_TIME == "MIDNIGHT" ~ ((hour(local_date_time) *3600) + (minute(local_date_time) * 60) + second(local_date_time))/60,
                                        EPISODE_REFERENCE_TIME == "START_OF_THE_SEGMENT" ~ (timestamp - segment_start_ts) / 60000)
            ) %>% 
      select(-consecutive, -level_diff, -mvpa_diff, -met_diff) %>% 
      group_by(local_segment) %>%
      nest() %>%
      mutate(sedentary = map(data, ~ episode_type_features(.x %>% filter(level == 0) , "sedentary", "level_episode_id")),
            lightlyactive = map(data, ~ episode_type_features(.x %>% filter(level == 1) , "lightlyactive", "level_episode_id")),
            fairlyactive = map(data, ~ episode_type_features(.x %>% filter(level == 2) , "fairlyactive", "level_episode_id")),
            veryactive = map(data, ~ episode_type_features(.x %>% filter(level == 3) , "veryactive", "level_episode_id")),
            mvpa = map(data, ~ episode_type_features(.x %>% filter(level >= 2) , "mvpa", "mvpa_episode_id")),
            lowmet = map(data, ~ episode_type_features(.x %>% filter(mets < MET_THRESHOLD) , "lowmet", "met_episode_id")),
            highmet = map(data, ~ episode_type_features(.x %>% filter(mets >= MET_THRESHOLD) , "highmet", "met_episode_id"))
            ) %>% 
      ungroup() %>% 
      select(all_of(c("local_segment", REQUESTED_EPISODES))) %>% 
      unnest(everything(), keep_empty=TRUE) %>% 
      select(starts_with(c("local_segment", REQUESTED_FEATURES)))
}