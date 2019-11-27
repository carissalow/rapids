source("packrat/init.R")

library("tidyverse")

gar <- read.csv(snakemake@input[[1]])

if(nrow(gar) > 0){
  activity_episodes <-
    gar %>%
    mutate(activity_episode = cumsum(c(1, head(activity_type, -1) != tail(activity_type, -1)))) %>% 
    group_by(activity_episode) %>%
    filter(n() > 1) %>%
    summarize(episode = first(activity_name),
              time_diff = (last(timestamp) - first(timestamp)) / (1000 * 60),
              local_start_date_time = first(local_date_time),
              local_end_date_time = last(local_date_time),
              local_start_date = first(local_date),
              local_end_date = last(local_date),
              local_start_day_segment = first(local_day_segment),
              local_end_day_segment = last(local_day_segment)) %>%
    select(-activity_episode)
} else {
  activity_episodes <- data.frame(episode = character(), 
                            time_diff = numeric(),
                            local_start_date_time = character(),
                            local_end_date_time = character(),
                            local_start_date = character(),
                            local_end_date = character(),
                            local_start_day_segment = character(),
                            local_end_day_segment = character())
}

write.csv(activity_episodes, snakemake@output[[1]], row.names = FALSE)
