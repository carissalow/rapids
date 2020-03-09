source("packrat/init.R")

library(tidyr)
library(purrr)
library(dplyr)

metric_files  <- snakemake@input[["metric_files"]]

metrics_of_single_participant <- metric_files %>%
  map(read.csv, stringsAsFactors = F, colClasses = c(local_date = "character")) %>%
  reduce(full_join, by="local_date")

write.csv(metrics_of_single_participant, snakemake@output[[1]], row.names = FALSE)