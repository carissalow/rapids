source("packrat/init.R")

library(dplyr)
library(tidyr)
library(stringr)

screen <- read.csv(snakemake@input[["screen"]])
participant_info <- snakemake@input[["participant_info"]]
platform <- readLines(participant_info, n=2)[[2]]

# Screen States
# Android: https://github.com/denzilferreira/aware-client/blob/78ccc22f0f822f8421bef9b1a73d36e71b8aa85b/aware-core/src/main/java/com/aware/Screen.java
# iOS: https://github.com/tetujin/aware-client-ios-v2/blob/master/Pods/AWAREFramework/AWAREFramework/Classes/Sensors/Screen/Screen.m#L120
# 0. OFF (Not existent for iOS due to API limitations)
# 1. ON (Not existent for iOS due to API limitations)
# 2. LOCKED
# 3. UNLOCKED

swap_screen_status <- function(data, status1, status2, time_buffer){
  # 800L is an auxiliary state to swap
  data %>% mutate(screen_status = case_when(screen_status == status2 & lag(screen_status) == status1 & (timestamp - lag(timestamp)) < time_buffer ~ 800L,
                                            TRUE ~ screen_status),
                  screen_status = case_when(screen_status == status1 & lead(screen_status) == 800L ~ status2,
                                            TRUE ~ screen_status),
                  screen_status = ifelse(screen_status == 800L, status1, screen_status))
}

get_ios_screen_episodes <- function(screen){
  episodes <- screen %>%
    # only keep consecutive pairs of 3,2 events
    filter( (screen_status == 3 & lead(screen_status) == 2) | (screen_status == 2 & lag(screen_status) == 3) ) %>%
    # in iOS and after our filtering, screen episodes should end with a LOCK event (2)
    mutate(episode_id = ifelse(screen_status == 2, 1:n(), NA_integer_)) %>%
    fill(episode_id, .direction = "updown") %>%
    group_by(episode_id) %>%
    summarise(episode = "unlock",
              screen_sequence = toString(screen_status),
              time_diff = (last(timestamp) - first(timestamp)) / 1000,
              local_start_date_time = first(local_date_time),
              local_end_date_time = last(local_date_time),
              local_start_date = first(local_date),
              local_end_date = last(local_date),
              local_start_day_segment = first(local_day_segment),
              local_end_day_segment = last(local_day_segment))
}

get_android_screen_episodes <- function(screen){
  episodes <- screen %>% 
    # filter out UNLOCK events (2) that come within 50 milliseconds of an ON or OFF event
    filter(!(screen_status == 2 & lag(screen_status) == 1 & timestamp - lag(timestamp) < 50)) %>% 
    filter(!(screen_status == 2 & lag(screen_status) == 0 & timestamp - lag(timestamp) < 50)) %>% 
    # in Android and after our filtering, screen episodes should end with a OFF event (0)
    mutate(episode_id = ifelse(screen_status == 0, 1:n(), NA_integer_)) %>% 
    fill(episode_id, .direction = "updown") %>% 
    group_by(episode_id)  %>% 
    # Rarely, UNLOCK events (3) get logged just before ON events (1). If this happens within 800ms, swap them
    swap_screen_status(3L, 1L, 800) %>% 
    # to be consistent with iOS we get rid off events (and thus sequences) starting with an ON (1) event
    filter(screen_status != 1) %>%
    summarise(episode = "unlock",
              screen_sequence = toString(screen_status),
              time_diff = (last(timestamp) - first(timestamp)) / 1000,
              local_start_date_time = first(local_date_time),
              local_end_date_time = last(local_date_time),
              local_start_date = first(local_date),
              local_end_date = last(local_date),
              local_start_day_segment = first(local_day_segment),
              local_end_day_segment = last(local_day_segment)) %>% 
    filter(str_detect(screen_sequence,
                      paste0("^(",
                            paste(c(3), collapse = "|"), # Filter sequences that start with 3 (UNLOCK) AND
                            ").*(",
                            paste(c(0), collapse = "|"), # Filter sequences that end with 0 (OFF)
                            ")$")))
}

if(nrow(screen) < 1){
  episodes <- data.frame(episode = character(), 
                                time_diff = numeric(),
                                local_start_date_time = character(),
                                local_end_date_time = character(),
                                local_start_date = character(),
                                local_end_date = character(),
                                local_start_day_segment = character(),
                                local_end_day_segment = character())
} else if(platform == "ios"){
  episodes <- get_ios_screen_episodes(screen)
} else if(platform == "android"){
  episodes <- get_android_screen_episodes(screen)
} else {
  print(paste0("The platform (second line) in ", participant_info, " should be android or ios"))
}

write.csv(episodes, snakemake@output[[1]], row.names = FALSE)
