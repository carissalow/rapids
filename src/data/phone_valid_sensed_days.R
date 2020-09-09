source("renv/activate.R")
library("dplyr")
library("tidyr")

phone_sensed_bins <- read.csv(snakemake@input[["phone_sensed_bins"]])
min_valid_hours_per_day <- as.integer(snakemake@params[["min_valid_hours_per_day"]])
min_valid_bins_per_hour <- as.integer(snakemake@params[["min_valid_bins_per_hour"]])
output_file <- snakemake@output[[1]]

phone_valid_sensed_days <- phone_sensed_bins %>% 
  pivot_longer(cols = -local_date, names_to = c("hour", "bin"), names_sep = "_") %>% 
  group_by(local_date, hour) %>%
  summarise(valid_bins = sum(value > 0)) %>% 
  group_by(local_date) %>% 
  summarise(valid_sensed_hours = sum(valid_bins >= min_valid_bins_per_hour)) %>% 
  mutate(is_valid_sensed_day = ifelse(valid_sensed_hours >= min_valid_hours_per_day, TRUE, FALSE))

write.csv(phone_valid_sensed_days, output_file, row.names = FALSE)
