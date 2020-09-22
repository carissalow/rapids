source("renv/activate.R")
library("dplyr")

battery <- read.csv(snakemake@input[[1]])

if(nrow(battery) > 0){
  # TODO expose this in the config file
  threshold_between_rows = 30
  battery_episodes <- battery %>% 
    filter(battery_status >= 2 ) %>% # discard unknown states
    mutate(start_timestamp = timestamp,
          end_timestamp = lead(start_timestamp) - 1,
          time_diff = (end_timestamp - start_timestamp) / 1000 / 60,
          time_diff = if_else(time_diff > threshold_between_rows, threshold_between_rows, time_diff),
          episode_id = 1:n()) %>%
    select(episode_id, start_timestamp, end_timestamp, battery_level)
} else {
  battery_episodes <- data.frame(episode_id = numeric(), 
                            start_timestamp = numeric(),
                            end_timestamp = character(),
                            battery_level = character())
}

write.csv(battery_episodes, snakemake@output[[1]], row.names = FALSE)
