source("renv/activate.R")
library("dplyr", warn.conflicts = F)
library(readr)
library(tidyr)

consecutive_threshold <- snakemake@params[["consecutive_threshold"]]
time_since_valid_location <- snakemake@params[["time_since_valid_location"]]
locations_to_use <- snakemake@params[["locations_to_use"]]
accuracy_limit <- snakemake@params[["accuracy_limit"]]

locations <- read.csv(snakemake@input[["locations"]]) %>% 
            filter(double_latitude != 0 & double_longitude != 0 & accuracy < accuracy_limit) %>% 
            drop_na(double_longitude, double_latitude) %>% 
            group_by(timestamp) %>% # keep only the row with the best accuracy if two or more have the same timestamp
            filter(accuracy == min(accuracy, na.rm=TRUE)) %>%  
            filter(row_number()==1) %>% 
            ungroup()

if(!locations_to_use %in% c("ALL", "FUSED_RESAMPLED", "GPS", "ALL_RESAMPLED")){
    print("Unkown location filter, provide one of the following three: ALL, GPS, ALL_RESAMPLED, or FUSED_RESAMPLED")
    quit(save = "no", status = 1, runLast = FALSE)
  }

# keep the location row that has the best (lowest) accuracy if more than 1 row was logged within any 1 second
if(locations_to_use %in% c("FUSED_RESAMPLED", "ALL_RESAMPLED"))
    locations <- locations %>% drop_na(double_longitude, double_latitude) %>% 
        mutate(minute_bin = timestamp %/% 1001) %>% 
        group_by(minute_bin) %>% 
        slice(which.min(accuracy)) %>% 
        ungroup() %>% 
        select(-minute_bin)

if(locations_to_use == "ALL"){
    processed_locations <- locations
} else if(locations_to_use == "GPS"){
    processed_locations <- locations %>% filter(provider == "gps")
} else if(locations_to_use %in% c("FUSED_RESAMPLED", "ALL_RESAMPLED")){
    if (locations_to_use == "FUSED_RESAMPLED"){
        locations <- locations %>% filter(provider == "fused")
        providers_to_keep = c("fused")
    } else if(locations_to_use == "ALL_RESAMPLED"){
        providers_to_keep = c("fused", "gps", "network")
    }

    if(nrow(locations) > 0){
        phone_sensed_timestamps  <- read_csv(snakemake@input[["phone_sensed_timestamps"]], col_types = cols_only(timestamp = col_double()))

        processed_locations <- locations %>%
            distinct(timestamp, .keep_all = TRUE) %>% 
            bind_rows(phone_sensed_timestamps) %>%
            arrange(timestamp) %>%
            # We group and therefore, fill in, missing rows that appear after a valid fused location record and exist
            # within consecutive_threshold minutes from each other
            mutate(consecutive_time_diff = c(1, diff(timestamp)),
                    resample_group = cumsum(!is.na(double_longitude) | consecutive_time_diff > (1000 * 60 * consecutive_threshold))) %>%
            group_by(resample_group) %>%
            # Filter those rows that are further away than time_since_valid_location since the last fused location
            mutate(time_from_fused = timestamp - first(timestamp)) %>%
            filter(provider %in% providers_to_keep | (time_from_fused < (1000 * 60 * time_since_valid_location))) %>%
            select(-consecutive_time_diff, -time_from_fused) %>% 
            # Summarise the period to resample for
            summarise(across(timestamp, max, .names = "limit"), across(everything(), first)) %>% 
            # the limit will be equal to the next timestamp-1 or the last binded timestamp (limit) plus the consecutive_threshold buffer
            # you can think of consecutive_threshold as the period a location row is valid for
            mutate(limit = pmin(lead(timestamp, default = 9999999999999) - 1, limit + (1000 * 60 * consecutive_threshold)),
                    n_resample = (limit - timestamp)%/%60001,
                    n_resample = n_resample + 1) %>% 
            drop_na(double_longitude, double_latitude) %>%
            uncount(weights = n_resample, .id = "id") %>% 
            mutate(provider = if_else(id > 1, "resampled", provider),
                    id = id -1,
                    timestamp = timestamp + (id * 60000)) %>% 
            ungroup() %>% 
            select(-resample_group, -limit, -id)
    } else {
        processed_locations <- locations
    }
}
write.csv(processed_locations,snakemake@output[[1]], row.names = F)
