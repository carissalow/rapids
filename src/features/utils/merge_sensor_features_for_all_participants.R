source("renv/activate.R")

library(tidyr)
library(purrr)
library("dplyr", warn.conflicts = F)
library(stringr)

feature_files  <- snakemake@input[["feature_files"]]


features_of_all_participants <- tibble(filename = feature_files) %>% # create a data frame
  mutate(file_contents = map(filename, ~ read.csv(., stringsAsFactors = F, colClasses = c(local_segment = "character", local_segment_label = "character", local_segment_start_datetime="character", local_segment_end_datetime="character"))),
         pid = str_match(filename, ".*/(.*)/all_sensor_features.csv")[,2]) %>%
  unnest(cols = c(file_contents)) %>%
  select(-filename)

write.csv(features_of_all_participants, snakemake@output[[1]], row.names = FALSE)