source("renv/activate.R")
library("dplyr", warn.conflicts = F)

battery <- read.csv(snakemake@input[[1]])
episode_threshold_between_rows <- snakemake@params[["episode_threshold_between_rows"]]

if(nrow(battery) > 0){
  episode_threshold_between_rows = episode_threshold_between_rows * 60000

  battery_episodes <- battery %>% 
    filter(battery_status >= 2 ) %>% # discard unknown states
    mutate(start_timestamp = timestamp, # a battery level starts as soon as is logged
          end_timestamp = lead(timestamp) - 1, # a battery level ends as soon as a new one is logged
          time_diff = (end_timestamp - start_timestamp),
          # we assume the current level existed until the next row only if that row is logged within [episode_threshold_between_rows] minutes
          end_timestamp = if_else(is.na(time_diff) | time_diff > (episode_threshold_between_rows), start_timestamp + (episode_threshold_between_rows), end_timestamp)) %>%
    mutate(time_diff = c(1, diff(start_timestamp)),
          level_diff = c(1, diff(battery_level)),
          status_diff = c(1, diff(battery_status)),
          episode_id = cumsum(level_diff != 0 | status_diff != 0 | time_diff > (episode_threshold_between_rows))) %>%
    group_by(episode_id) %>%
    summarise(device_id = first(device_id), battery_level = first(battery_level), battery_status = first(battery_status), start_timestamp=first(start_timestamp), end_timestamp = last(end_timestamp))
} else {
  battery_episodes <- data.frame(device_id = character(),
                            episode_id = numeric(), 
                            start_timestamp = numeric(),
                            end_timestamp = character(),
                            battery_level = character(),
                            battery_status = character())
}

write.csv(battery_episodes, snakemake@output[[1]], row.names = FALSE)
