source("renv/activate.R")
library(tidyr)
library(dplyr)

filter_participant_without_enough_days <- function(clean_features, days_before_threshold, days_after_threshold){
  clean_features$day_type <- ifelse(clean_features$day_idx < 0, -1, ifelse(clean_features$day_idx > 0, 1, 0))
  if("pid" %in% colnames(clean_features)){
    clean_features <- clean_features %>% 
      group_by(pid) %>% 
      add_count(pid, day_type) # this adds a new column "n"
  } else {
    clean_features <- clean_features %>% add_count(day_type < 0)
  }

  # Only keep participants with enough days before surgery and after discharge
  clean_features <- clean_features %>% 
    mutate(count_before = ifelse(day_type == -1, n, NA), # before surgery
          count_after = ifelse(day_type == 1, n, NA)) %>%  # after discharge
    fill(count_before, .direction = "downup") %>% 
    fill(count_after, .direction = "downup") %>% 
    filter(count_before >= days_before_threshold & count_after >= days_after_threshold) %>% 
    select(-n, -count_before, -count_after, -day_type) %>% 
    ungroup()

  return(clean_features)
}

clean_features <- read.csv(snakemake@input[[1]])
cols_nan_threshold <- as.numeric(snakemake@params[["cols_nan_threshold"]])
drop_zero_variance_columns <- as.logical(snakemake@params[["cols_var_threshold"]])
rows_nan_threshold <- as.numeric(snakemake@params[["rows_nan_threshold"]])
days_before_threshold <- as.numeric(snakemake@params[["days_before_threshold"]])
days_after_threshold <- as.numeric(snakemake@params[["days_after_threshold"]])
features_exclude_day_idx <- as.logical(snakemake@params[["features_exclude_day_idx"]])


# We have to do this before and after dropping rows, that's why is duplicated
clean_features <- filter_participant_without_enough_days(clean_features, days_before_threshold, days_after_threshold)

# drop columns with a percentage of NA values above cols_nan_threshold
if(nrow(clean_features))
    clean_features <- clean_features %>% select_if(~ sum(is.na(.)) / length(.) <= cols_nan_threshold )

if(drop_zero_variance_columns)
  clean_features <- clean_features %>% select_if(grepl("pid|local_date",names(.)) | sapply(., n_distinct, na.rm = T) > 1)

# drop rows with a percentage of NA values above rows_nan_threshold
clean_features <- clean_features %>% 
  mutate(percentage_na =  rowSums(is.na(.)) / ncol(.)) %>% 
  filter(percentage_na < rows_nan_threshold) %>% 
  select(-percentage_na)

if(nrow(clean_features) != 0){
  clean_features <- filter_participant_without_enough_days(clean_features, days_before_threshold, days_after_threshold)
  
  # include "day_idx" as features or not
  if(features_exclude_day_idx)
    clean_features <- clean_features %>% select(-day_idx)
}

write.csv(clean_features, snakemake@output[[1]], row.names = FALSE)
