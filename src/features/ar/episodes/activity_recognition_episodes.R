source("renv/activate.R")
library("dplyr")

activity_recognition <- read.csv(snakemake@input[[1]])

if(nrow(activity_recognition) > 0){
  threshold_between_rows = 5
  ar_episodes <- activity_recognition %>% 
    mutate(start_timestamp = timestamp,
          end_timestamp = lead(start_timestamp) - 1,
          time_diff = (end_timestamp - start_timestamp) / 1000 / 60,
          time_diff = if_else(time_diff > threshold_between_rows, threshold_between_rows, time_diff),
          episode_id = 1:n()) %>%
    select(episode_id, start_timestamp, end_timestamp, activity_type)

} else {
  ar_episodes <- data.frame(start_timestamp = numeric(), 
                            end_timestamp = numeric(),
                            episode_id = numeric())
}

write.csv(ar_episodes, snakemake@output[[1]], row.names = FALSE)
