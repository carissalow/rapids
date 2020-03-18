source("packrat/init.R")
library(tidyr)
library(dplyr)

filter_participant_without_enough_days <- function(clean_metrics, participants_day_threshold){
  if("pid" %in% colnames(clean_metrics))
    clean_metrics <- clean_metrics %>% group_by(pid)
  
  clean_metrics <- clean_metrics %>% 
    filter(n() >= participants_day_threshold) %>% 
    ungroup()
  
  return(clean_metrics)
}

clean_metrics <- read.csv(snakemake@input[[1]])
cols_nan_threshold <- snakemake@params[["cols_nan_threshold"]]
drop_zero_variance_columns <- snakemake@params[["cols_var_threshold"]]
rows_nan_threshold <- snakemake@params[["rows_nan_threshold"]]
participants_day_threshold <- snakemake@params[["participants_day_threshold"]]

# We have to do this before and after dropping rows, that's why is duplicated
clean_metrics <- filter_participant_without_enough_days(clean_metrics, participants_day_threshold)

# drop columns with a percentage of NA values above cols_nan_threshold
if(nrow(clean_metrics))
    clean_metrics <- clean_metrics %>% select_if(~ sum(is.na(.)) / length(.) <= cols_nan_threshold )

if(drop_zero_variance_columns)
  clean_metrics <- clean_metrics %>% select_if(grepl("pid|local_date",names(.)) | sapply(., n_distinct) > 1)

# drop rows with a percentage of NA values above rows_nan_threshold
clean_metrics <- clean_metrics %>% 
  mutate(percentage_na =  rowSums(is.na(.)) / ncol(.)) %>% 
  filter(percentage_na < rows_nan_threshold) %>% 
  select(-percentage_na)

clean_metrics <- filter_participant_without_enough_days(clean_metrics, participants_day_threshold)

write.csv(clean_metrics, snakemake@output[[1]], row.names = FALSE)
