source("packrat/init.R")

library("tidyverse")

screen <- read.csv(snakemake@input[[1]])

if(nrow(screen) > 0){
  unlock_episodes <-
    screen %>%
    # in iOS there are unlock (3) events on the same second, discard them
    distinct(screen_status, utc_date_time, .keep_all = TRUE) %>%
    # in Android we discard on and off events (0,1) for now (iOS does not collect them)
    filter(screen_status == 2 | screen_status == 3) %>% 
    # create groups of consecutive unlock/lock (3/2) events
    mutate(screen_episode = cumsum(c(1, head(screen_status, -1) == 2 & tail(screen_status, -1) == 3))) %>% 
    group_by(screen_episode) %>%
    # in Android there are multiple consecutive unlock/lock events so we keep the closest pair
    # this happens because ACTION_SCREEN_OFF and ON are "sent when the device becomes non-interactive 
    # which may have nothing to do with the screen turning off" see:
    # https://developer.android.com/reference/android/content/Intent.html#ACTION_SCREEN_OFF
    filter((screen_status == 2 & screen_status != lag(screen_status, default="1")) |
            (screen_status == 3 & screen_status != lead(screen_status, default="1"))) %>%
    filter(n() == 2) %>%
    summarize(episode = "unlock",
              time_diff = (last(timestamp) - first(timestamp)) / (1000 * 60),
              local_start_date_time = first(local_date_time),
              local_end_date_time = last(local_date_time),
              local_start_date = first(local_date),
              local_end_date = last(local_date),
              local_start_day_segment = first(local_day_segment),
              local_end_day_segment = last(local_day_segment)) %>%
    select(-screen_episode)
} else {
  unlock_episodes <- data.frame(episode = character(), 
                            time_diff = numeric(),
                            local_start_date_time = character(),
                            local_end_date_time = character(),
                            local_start_date = character(),
                            local_end_date = character(),
                            local_start_day_segment = character(),
                            local_end_day_segment = character())
}

write.csv(unlock_episodes, snakemake@output[[1]], row.names = FALSE)
