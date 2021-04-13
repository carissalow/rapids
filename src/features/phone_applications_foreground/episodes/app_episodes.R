source("renv/activate.R")

library("dplyr", warn.conflicts = F)
options(scipen=999)

screen_ep <- read.csv(snakemake@input[["screen"]]) %>% mutate(timestamp = end_timestamp)

app <- read.csv(snakemake@input[["app"]]) %>% 
            select(c('timestamp', 'device_id', 'package_name', 'application_name', 'is_system_app', 'genre'))

if (nrow(screen_ep) > 0 & nrow(app) > 0){
    
    screen_ep <- screen_ep %>% select(-c('episode_id','episode','screen_sequence','device_id'))
    joined_dt <- full_join(app, screen_ep, by = "timestamp") %>% 
                 arrange("timestamp") %>%  
                 mutate(start_timestamp = timestamp, end_timestamp = lead(timestamp)) %>% 
                 filter(!is.na(application_name)) %>% 
                 select(-c('timestamp')) %>% head(-1)
    write.csv(joined_dt, snakemake@output[[1]], row.names = FALSE)
}