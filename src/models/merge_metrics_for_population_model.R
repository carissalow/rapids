source("packrat/init.R")

library(tidyr)
library(purrr)
library(dplyr)
library(stringr)

metric_files  <- snakemake@input[["metric_files"]]

metrics_of_all_participants <- data_frame(filename = metric_files) %>% # create a data frame
  mutate(file_contents = map(filename, ~ read.csv(., stringsAsFactors = F, colClasses = c(local_date = "character"))),
         pid = str_match(filename, ".*/(p[0-9]{2})/.*")[,2]) %>%
  unnest() %>%
  select(-filename)

write.csv(metrics_of_all_participants, snakemake@output[[1]], row.names = FALSE)