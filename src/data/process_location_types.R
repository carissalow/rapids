source("renv/activate.R")
library(dplyr)
library(readr)
library(tidyr)

source("src/data/assign_to_day_segment.R")

bin_size <- snakemake@params[["bin_size"]]
timezone <- snakemake@params[["timezone"]]
consecutive_threshold <- snakemake@params[["consecutive_threshold"]]
time_since_valid_location <- snakemake@params[["time_since_valid_location"]]
location_to_used <- snakemake@params[["time_since_valocation_to_usedlid_location"]]
day_segments <- read.csv(snakemake@input[["day_segments"]])
day_segments_type <- snakemake@params[["day_segments_type"]]

phone_sensed_bins  <- read_csv(snakemake@input[["phone_sensed_bins"]], col_types = cols(local_date = col_character()))
locations <- read_csv(snakemake@input[["locations"]], col_types = cols()) %>% filter(provider == "fused") %>% 
            filter(double_latitude != 0 & double_longitude != 0)

locations_to_use <- snakemake@params["locations_to_use"]
if(!locations_to_use %in% c("ALL", "FUSED_RESAMPLED", "GPS")){
    print("Unkown location filter, provide one of the following three: ALL, GPS, or FUSED_RESAMPLED")
    quit(save = "no", status = 1, runLast = FALSE)
  }


if(locations_to_use == "ALL"){
    processed_locations <- locations
} else if(locations_to_use == "GPS"){
    processed_locations <- locations %>% filter(provider == "gps")
} else if(locations_to_use == "FUSED_RESAMPLED"){
    locations <- locations %>% filter(provider == "fused")
    if(nrow(locations) > 0){
        sensed_minute_bins <- phone_sensed_bins %>% 
            pivot_longer(-local_date, names_to = c("hour", "bin"), names_sep = "_", values_to = "sensor_count") %>% 
            mutate(hour = as.integer(hour), bin = as.integer(bin)) %>% 
            complete(nesting(local_date, hour), bin = seq(0, 59,1)) %>% 
            fill(sensor_count) %>% 
            mutate(timestamp = as.numeric(as.POSIXct(paste0(local_date, " ", hour,":", bin,":00"), format = "%Y-%m-%d %H:%M:%S", tz = timezone)) * 1000 ) %>%
            filter(sensor_count > 0) %>% 
            select(timestamp)

        resampled_locations <- locations %>%
            select(-assigned_segments) %>% 
            bind_rows(sensed_minute_bins) %>%
            mutate(provider = replace_na(provider, "resampled"))  %>% 
            arrange(timestamp) %>% 
            # We group and therefore, fill in, missing rows that appear after a valid fused location record and exist
            # within consecutive_threshold minutes from each other
            mutate(consecutive_time_diff = c(1, diff(timestamp)),
                resample_group = cumsum(!is.na(double_longitude) | consecutive_time_diff > (1000 * 60 * consecutive_threshold))) %>% 
            group_by(resample_group) %>% 
            # drop rows that are logged after time_since_valid_location minutes from the last valid fused location
            filter((timestamp - first(timestamp) < (1000 * 60 * time_since_valid_location))) %>% 
            fill(-timestamp, -resample_group) %>% 
            select(-consecutive_time_diff) %>% 
            drop_na(double_longitude, double_latitude, accuracy) %>% 
            # Add local date_time
            mutate(utc_date_time = as.POSIXct(timestamp/1000, origin="1970-01-01", tz="UTC"),
                local_date_time = format(utc_date_time, tz = timezone, usetz = F)) %>% 
            separate(local_date_time, c("local_date","local_time"), "\\s", remove = FALSE) %>% 
            separate(local_time, c("local_hour", "local_minute"), ":", remove = FALSE, extra = "drop") %>%
            mutate(local_hour = as.numeric(local_hour),
                            local_minute = as.numeric(local_minute)) %>%
            # Delete resampled rows that exist in the same minute as other original (fused) rows
            group_by(local_date, local_hour, local_minute) %>%
            mutate(n = n()) %>%
            filter(n == 1 | (n > 1 & provider == "fused")) %>% 
            select(-n) %>% 
            ungroup()
        processed_locations <- assign_to_day_segment(resampled_locations, day_segments, day_segments_type, timezone)
    } else {
        processed_locations <- locations
    }
}
write.csv(processed_locations,snakemake@output[[1]], row.names = F)
