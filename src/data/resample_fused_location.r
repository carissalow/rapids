source("packrat/init.R")


library(dplyr)
library(readr)
library(tidyr)

bin_size <- snakemake@params[["bin_size"]]
timezone <- snakemake@params[["timezone"]]
consecutive_threshold <- snakemake@params[["consecutive_threshold"]]
time_since_valid_location <- snakemake@params[["time_since_valid_location"]]

locations <- read_csv(snakemake@input[["locations"]], col_types = cols()) %>% filter(provider == "fused")
phone_sensed_bins  <- read_csv(snakemake@input[["phone_sensed_bins"]], col_types = cols(local_date = col_character()))

if(nrow(locations) > 0){
    sensed_minute_bins <- phone_sensed_bins %>% 
        pivot_longer(-local_date, names_to = c("hour", "bin"), names_ptypes = list(hour = integer(), bin = integer()), names_sep = "_", values_to = "sensor_count") %>% 
        complete(nesting(local_date, hour), bin = seq(0, 59,1)) %>% 
        fill(sensor_count) %>% 
        mutate(timestamp = as.numeric(as.POSIXct(paste0(local_date, " ", hour,":", bin,":00"), format = "%Y-%m-%d %H:%M:%S", tz = timezone)) * 1000 ) %>%
        filter(sensor_count > 0) %>% 
        select(timestamp)

    resampled_locations <- locations %>%
        bind_rows(sensed_minute_bins) %>% 
        arrange(timestamp) %>% 
        # We group and therefore, fill in, missing rows that appear after a valid fused location record and exist
        # within consecutive_threshold minutes from each other
        mutate(consecutive_time_diff = c(1, diff(timestamp)),
            resample_group = cumsum(!is.na(double_longitude) | consecutive_time_diff > (1000 * 60 * consecutive_threshold))) %>% 
        group_by(resample_group) %>% 
        # drop rows that are logged after time_since_valid_location hours from the last valid fused location
        filter((timestamp - first(timestamp) < (1000 * 60 * 60 * time_since_valid_location))) %>% 
        fill(-timestamp, -resample_group) %>% 
        select(-consecutive_time_diff) %>% 
        drop_na(double_longitude, double_latitude, accuracy)

    write.csv(resampled_locations,snakemake@output[[1]], row.names = F)
} else {
    write.csv(locations,snakemake@output[[1]], row.names = F)
}
