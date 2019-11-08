source("packrat/init.R")

library("tidyverse")

battery <- read.csv(snakemake@input[[1]])

consumption <- battery %>%
    mutate(group = ifelse(lag(battery_status) != battery_status, 1, 0) %>% coalesce(0),
            group_id = cumsum(group) + 1) %>%
    filter(battery_status == 2 || battery_status == 3) %>%
    group_by(group_id) %>%
    summarize(battery_diff = first(battery_level) - last(battery_level),
            time_diff = (last(timestamp) - first(timestamp)) / (1000 * 60 * 60),
            local_start_date_time = first(local_date_time),
            local_end_date_time = last(local_date_time),
            local_start_date = first(local_date),
            local_end_date = last(local_date)) %>%
    select(-group_id) %>%
    filter(time_diff > 0.1) # Avoids including quick cycles

write.csv(consumption, snakemake@output[[1]], row.names = FALSE)
