source("renv/activate.R")

library("dplyr", warn.conflicts = F)
library(tidyr)
library(stringr)

screen <- read.csv(snakemake@input[["screen"]])
participant_info <- snakemake@input[["participant_info"]]

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

get_screen_episodes <- function(screen){  
  # Aware Android logs LOCK events after turning the screen ON or OFF but we filter them out to simplify this analysis. 
  # The code below only process UNLOCK to OFF episodes, but it's possible to modify it for ON to OFF (see line 61) or ON to UNLOCK episodes.

  episodes <- screen %>% 
    # Relevant for Android. Remove LOCK events (2) that come within 50 milliseconds of an ON (1) or OFF (0) event
    filter(!(screen_status == 2 & lag(screen_status) == 1 & timestamp - lag(timestamp) < 50)) %>% 
    filter(!(screen_status == 2 & lag(screen_status) == 0 & timestamp - lag(timestamp) < 50)) %>% 
    # After our filtering, screen episodes should end with a OFF event (0)
    mutate(episode_id = ifelse(screen_status == 0, 1:n(), NA_integer_)) %>% 
    fill(episode_id, .direction = "updown") %>% 
    group_by(episode_id)  %>% 
    # Relevant for Android. Rarely, UNLOCK events (3) get logged just before ON events (1). If this happens within 800ms, swap them
    swap_screen_status(3L, 1L, 800) %>% 
    # Relevant for Android. To be consistent with iOS we remove events (and thus sequences) starting with an ON (1) event
    filter(screen_status != 1) %>%
    # Only keep consecutive 3,0 pairs (UNLOCK, OFF)
    filter( (screen_status == 3 & lead(screen_status) == 0) | (screen_status == 0 & lag(screen_status) == 3) ) %>%
    summarise(episode = "unlock",
              device_id = first(device_id),
              screen_sequence = toString(screen_status),
              start_timestamp = first(timestamp),
              end_timestamp = last(timestamp)) %>% 
    filter(str_detect(screen_sequence,
                      paste0("^(",
                            paste(c(3), collapse = "|"), # Filter sequences that start with 3 (UNLOCK) AND
                            ").*(",
                            paste(c(0), collapse = "|"), # Filter sequences that end with 0 (OFF)
                            ")$")))
}

if(nrow(screen) < 2){
  episodes <- data.frame(device_id = character(),
                                episode = character(), 
                                screen_sequence = character(),
                                start_timestamp = character(),
                                end_timestamp = character())
} else {
  episodes <- get_screen_episodes(screen)
}

write.csv(episodes, snakemake@output[[1]], row.names = FALSE)
