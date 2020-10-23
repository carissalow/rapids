source("renv/activate.R")

library("dplyr", warn.conflicts = F)
library(tidyr)

all_sensors <- snakemake@input[["all_sensors"]]
bin_size <- snakemake@params[["bin_size"]]
output_file <- snakemake@output[[1]]

# Load all sensors and extract timestamps
all_sensor_data <- data.frame(timestamp = c())
for(sensor in all_sensors){
  sensor_data <- read.csv(sensor, stringsAsFactors = F) %>% 
    select(local_date, local_hour, local_minute) %>% 
    mutate(sensor = basename(sensor))
  all_sensor_data <- rbind(all_sensor_data, sensor_data)
}

phone_sensed_bins <- all_sensor_data %>% 
  mutate(bin = (local_minute %/% bin_size) * bin_size) %>% # bin rows into bin_size-minute bins
  group_by(local_date, local_hour, bin) %>% 
  summarise(sensor_count = n_distinct(sensor)) %>%
  ungroup() %>% 
  complete(nesting(local_date), 
           local_hour = seq(0, 23, 1), 
           bin = seq(0, (59 %/% bin_size) * bin_size, bin_size), 
           fill = list(sensor_count=0)) %>% 
  pivot_wider(names_from = c(local_hour, bin), values_from = sensor_count)

write.csv(phone_sensed_bins, output_file, row.names = FALSE)

