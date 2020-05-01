source("renv/activate.R")

library(dplyr)

all_sensors <- snakemake@input[["all_sensors"]]
bin_size <- snakemake@params[["bin_size"]]
min_valid_hours <- snakemake@params[["min_valid_hours"]]
min_bins_per_hour <- snakemake@params[["min_bins_per_hour"]]
output_file <- snakemake@output[[1]]

# Load all sensors and extract timestamps
all_sensor_data <- data.frame(timestamp = c())
for(sensor in all_sensors){
    sensor_data <- read.csv(sensor, stringsAsFactors = F) %>% select(local_date, local_hour, local_minute)
    all_sensor_data <- rbind(all_sensor_data, sensor_data)
}

phone_valid_sensed_days <- all_sensor_data %>% 
    mutate(bin = (local_minute %/% bin_size) * bin_size) %>% # bin rows into bin_size-minute bins
    group_by(local_date, local_hour, bin) %>% 
    summarise(minute_period = first(bin)) %>% #filter repeated bins (if rows were logged within bin_size minutes)
    ungroup() %>% 
    group_by(local_date, local_hour) %>% 
    summarise(bins = n()) %>% # Count how many bins there are per hour
    ungroup() %>% 
    filter(bins >= min_bins_per_hour) %>% # Discard those hours where there were fewer than min_bins_per_hour
    group_by(local_date) %>% 
    summarise(valid_hours = n()) %>% # Count how many valid hours each day has
    filter(valid_hours >= min_valid_hours) # Discard those days where there were fewer than min_valid_hours

write.csv(phone_valid_sensed_days, output_file, row.names = FALSE)
