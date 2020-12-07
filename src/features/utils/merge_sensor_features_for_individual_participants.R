source("renv/activate.R")

library(tidyr)
library(purrr)
library("dplyr", warn.conflicts = F)
library("methods")
library("mgm")
library("qgraph")
library("dplyr", warn.conflicts = F)
library("scales")
library("ggplot2")
library("purrr")
library("tidyr")
library("reshape2")

feature_files  <- snakemake@input[["feature_files"]]

features_for_individual_model <- feature_files %>%
  map(read.csv, stringsAsFactors = F, colClasses = c(local_segment = "character", local_segment_label = "character", local_segment_start_datetime="character", local_segment_end_datetime="character")) %>%
  reduce(full_join, by=c("local_segment","local_segment_label","local_segment_start_datetime","local_segment_end_datetime"))

write.csv(features_for_individual_model, snakemake@output[[1]], row.names = FALSE)
