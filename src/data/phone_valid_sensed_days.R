source("renv/activate.R")
library("dplyr")
library("tidyr")

phone_sensed_bins <- read.csv(snakemake@input[["phone_sensed_bins"]])
min_valid_hours_per_day <- snakemake@params[["min_valid_hours_per_day"]]
min_valid_bins_per_hour <- snakemake@params[["min_valid_bins_per_hour"]]
output_file <- snakemake@output[[1]]

phone_valid_sensed_days <- phone_sensed_bins %>% 
  pivot_longer(cols = -local_date, names_to = c("hour", "bin"), names_sep = "_") %>% 
  filter(value > 0) %>%
  group_by(local_date, hour) %>%
  summarise(valid_bins = n()) %>% 
  filter(valid_bins >= min_valid_bins_per_hour) %>% 
  group_by(local_date) %>% 
  summarise(valid_sensed_hours = n()) %>% 
  mutate(is_valid_sensed_day = ifelse(valid_sensed_hours >= min_valid_hours_per_day, TRUE, FALSE))

write.csv(phone_valid_sensed_days, output_file, row.names = FALSE)
