source("renv/activate.R")
library("dplyr", warn.conflicts = F)
library(readr)
library(tidyr)
library(purrr)
options(scipen=999)

all_sensors = snakemake@input[["all_sensors"]]

sensor_timestamps <- tibble(files = all_sensors) %>% 
  mutate(timestamps = map(files,~ read_csv(.,col_types = cols_only(timestamp = col_double(), device_id = col_character()))),
         sensor = row_number(),
         files = NULL) %>% 
  unnest(timestamps) %>%
  mutate(timestamp = (timestamp %/% 1000) * 1000) %>% 
  distinct(timestamp, .keep_all = TRUE) %>%
  arrange(timestamp)

write.csv(sensor_timestamps, snakemake@output[[1]], row.names = FALSE)