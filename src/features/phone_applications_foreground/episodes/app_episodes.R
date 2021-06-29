source("renv/activate.R")

library("dplyr", warn.conflicts = F)
library("tidyverse")
options(scipen=999)

screen_ep <- read_csv(snakemake@input[["screen"]], col_types = cols_only(start_timestamp = col_double(), end_timestamp = col_double())) %>% 
             mutate(timestamp = end_timestamp)

app <- read_csv(snakemake@input[["app"]], col_types = cols_only(timestamp = col_double(), device_id = col_character(), package_name = col_character(), application_name = col_character(), is_system_app = col_integer(), genre = col_character()))

if (nrow(screen_ep) > 0 & nrow(app) > 0){
        
    joined_dt <- full_join(app, screen_ep, by = "timestamp") %>% 
                 arrange(timestamp) %>%  
                 mutate(start_timestamp = timestamp, end_timestamp = lead(timestamp)) %>% 
                 filter(!is.na(application_name)) %>% 
                 select(-c('timestamp')) %>% head(-1) %>% 
                 # converting the duration from milliseconds to minutes
                 mutate(duration = (end_timestamp - start_timestamp)/(1000*60))
    write.csv(joined_dt, snakemake@output[[1]], row.names = FALSE)
} else {
    empty <- tibble(device_id = character(), package_name = character(), application_name = character(), is_system_app = integer(), genre = character(), start_timestamp = double(), end_timestamp = double(), duration = double())
    write.csv(empty, snakemake@output[[1]], row.names = FALSE)
}