source("packrat/init.R")

library("tidyverse")

input <- read.csv(snakemake@input[[1]])
sensor_output <- snakemake@output[[1]]
timezone_periods <- snakemake@params[["timezone_periods"]]
fixed_timezone <- snakemake@params[["fixed_timezone"]]

if(!is.null(timezone_periods)){
    timezones <- read.csv(timezones)
    tz_starts <- timezones$start
    output <- input %>% 
                mutate(timezone = findInterval(timestamp / 1000, tz_starts), # Set an interval ID based on timezones' start column
                        timezone = ifelse(timezone == 0, 1, timezone), # Correct the first timezone ID
                        timezone = recode(timezone, !!! timezones$timezone), # Swap IDs for text labels
                        timezone = as.character(timezone)) %>%
                rowwise() %>%
                mutate(utc_date_time = as.POSIXct(timestamp/1000, origin="1970-01-01", tz="UTC"),
                        local_date_time = format(utc_date_time, tz = timezone, usetz = F),
                        local_date = as.Date(local_date_time),
                        local_time = strsplit(local_date_time, " ")[[1]][2],
                        local_hour = as.numeric(strsplit(local_time, ":")[[1]][1]),
                        day_segment = case_when(local_hour %in% 0:5 ~ "night",
                                                local_hour %in% 6:11 ~ "morning",
                                                local_hour %in% 12:17 ~ "afternoon",
                                                local_hour %in% 18:23 ~ "evening"))

    write.csv(output, sensor_output)
} else if(!is.null(fixed_timezone)){
    output <- input %>% 
                mutate(utc_date_time = as.POSIXct(timestamp/1000, origin="1970-01-01", tz="UTC"),
                        local_date_time = format(utc_date_time, tz = fixed_timezone, usetz = F),
                        local_date = as.Date(local_date_time),
                        local_time = strsplit(local_date_time, " ")[[1]][2],
                        local_hour = as.numeric(strsplit(local_time, ":")[[1]][1]),
                        local_day_segment = case_when(local_hour %in% 0:5 ~ "night",
                                                local_hour %in% 6:11 ~ "morning",
                                                local_hour %in% 12:17 ~ "afternoon",
                                                local_hour %in% 18:23 ~ "evening"))

    write.csv(output, sensor_output)
}
