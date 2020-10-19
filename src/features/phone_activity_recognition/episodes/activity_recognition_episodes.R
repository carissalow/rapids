source("renv/activate.R")
library("dplyr")

activity_recognition <- read.csv(snakemake@input[[1]])

if(nrow(activity_recognition) > 0){
  # TODO expose this in the config file
  threshold_between_rows = 5 * 60000

  ar_episodes <- activity_recognition %>% 
  mutate(start_timestamp = timestamp, # a battery level starts as soon as is logged
         time_diff = (lead(timestamp) - start_timestamp), # lead diff
         # we assume the current activity existed until the next row only if that row is logged within [threshold_between_rows] minutes
         end_timestamp = if_else(is.na(time_diff) | time_diff > (threshold_between_rows), start_timestamp + (threshold_between_rows), lead(timestamp) - 1), 
         time_diff = c(1, diff(start_timestamp)), # lag diff
         type_diff = c(1, diff(activity_type)),
         episode_id = cumsum(type_diff != 0 | time_diff > (threshold_between_rows))) %>% 
  group_by(episode_id) %>%
  summarise(activity_name = first(activity_name), activity_type = first(activity_type), start_timestamp=first(start_timestamp), end_timestamp = last(end_timestamp))

} else {
  ar_episodes <- data.frame(start_timestamp = numeric(), 
                            end_timestamp = numeric(),
                            episode_id = numeric(),
                            activity_type = numeric(),
                            activity_name = character())
}

write.csv(ar_episodes, snakemake@output[[1]], row.names = FALSE)
