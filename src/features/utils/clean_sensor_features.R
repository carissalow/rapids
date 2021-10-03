source("renv/activate.R")
library(tidyr)
library("dplyr", warn.conflicts = F)
library(tidyverse)
library(caret)
library(corrr)


clean_features <- read.csv(snakemake@input[[1]])
cols_nan_threshold <- as.numeric(snakemake@params[["cols_nan_threshold"]])
drop_zero_variance_columns <- as.logical(snakemake@params[["cols_var_threshold"]])
rows_nan_threshold <- as.numeric(snakemake@params[["rows_nan_threshold"]])
data_yielded_hours_ratio_threshold <- as.numeric(snakemake@params[["data_yielded_hours_ratio_threshold"]])
corr_valid_pairs_threshold <- as.numeric(snakemake@params[["corr_valid_pairs_threshold"]])
corr_threshold <- as.numeric(snakemake@params[["corr_threshold"]])

# drop rows with the value of "phone_data_yield_rapids_ratiovalidyieldedhours" column less than data_yielded_hours_ratio_threshold
clean_features <- clean_features %>% 
  filter(phone_data_yield_rapids_ratiovalidyieldedhours > data_yielded_hours_ratio_threshold)

# drop columns with a percentage of NA values above cols_nan_threshold
if(nrow(clean_features))
    clean_features <- clean_features %>% select_if(~ sum(is.na(.)) / length(.) <= cols_nan_threshold )

if(drop_zero_variance_columns)
  clean_features <- clean_features %>% select_if(grepl("pid|local_segment|local_segment_label|local_segment_start_datetime|local_segment_end_datetime",names(.)) | sapply(., n_distinct, na.rm = T) > 1)

# drop highly correlated features
features_for_corr <- clean_features %>% 
  select_if(is.numeric) %>% 
  select_if(sapply(., n_distinct, na.rm = T) > 1)

valid_pairs <- crossprod(!is.na(features_for_corr)) >= corr_valid_pairs_threshold * nrow(features_for_corr)

highly_correlated_features <- features_for_corr %>% 
  correlate(use = "pairwise.complete.obs", method = "spearman") %>% 
  column_to_rownames(., var = "term") %>% 
  as.matrix() %>% 
  replace(!valid_pairs | is.na(.), 0) %>% 
  findCorrelation(., cutoff = corr_threshold, verbose = F, names = T)

clean_features <- clean_features[, !names(clean_features) %in% highly_correlated_features]

# drop rows with a percentage of NA values above rows_nan_threshold
clean_features <- clean_features %>% 
  mutate(percentage_na =  rowSums(is.na(.)) / ncol(.)) %>% 
  filter(percentage_na < rows_nan_threshold) %>% 
  select(-percentage_na)

write.csv(clean_features, snakemake@output[[1]], row.names = FALSE)
