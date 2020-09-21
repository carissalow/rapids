source("renv/activate.R")

library(dplyr)
library(tidyr)
library(lubridate)

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

if(nrow(all_sensor_data) == 0){
  bins = seq(0, 59, by = bin_size)
  hours = seq(0, 23, 1)
  write.csv(crossing(hours, bins) %>% unite("hour_bin",hours, bins, sep = "_") %>% mutate(value = NA, local_date = NA) %>% pivot_wider(names_from = hour_bin, values_from=value) %>% head(0), output_file, row.names = FALSE)
} else{
  phone_sensed_bins <- all_sensor_data %>% 
    mutate(bin = (local_minute %/% bin_size) * bin_size) %>% # bin rows into bin_size-minute bins
    group_by(local_date, local_hour, bin) %>% 
    summarise(sensor_count = n_distinct(sensor)) %>%
    ungroup() %>% 
    mutate(local_date = lubridate::ymd(local_date)) %>% 
    complete(local_date = seq.Date(min(local_date), max(local_date), by="day"), 
            fill = list(local_hour = 0, bin = 0, sensor_count = 0)) %>% 
    complete(nesting(local_date), 
            local_hour = seq(0, 23, 1), 
            bin = seq(0, 59, bin_size), 
            fill = list(sensor_count=0)) %>% 
    pivot_wider(names_from = c(local_hour, bin), values_from = sensor_count)

  write.csv(phone_sensed_bins, output_file, row.names = FALSE)
}
