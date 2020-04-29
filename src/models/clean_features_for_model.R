source("packrat/init.R")
library(tidyr)
library(dplyr)

filter_participant_without_enough_days <- function(clean_features, days_before_threshold, days_after_threshold){
  if("pid" %in% colnames(clean_features)){
    clean_features <- clean_features %>% 
      group_by(pid) %>% 
      add_count(pid, day_type) # this adds a new column "n"
  } else {
    clean_features <- clean_features %>% add_count(day_type)
  }

  # Only keep participants with enough days before surgery and after discharge
  clean_features <- clean_features %>% 
    mutate(count_before = ifelse(day_type == -1, n, NA), # before surgery
          count_after = ifelse(day_type == 1, n, NA)) %>%  # after discharge
    fill(count_before, .direction = "downup") %>% 
    fill(count_after, .direction = "downup") %>% 
    filter(count_before >= days_before_threshold & count_after >= days_after_threshold) %>% 
    select(-n, -count_before, -count_after) %>% 
    ungroup()

  return(clean_features)
}

clean_features <- read.csv(snakemake@input[[1]])
cols_nan_threshold <- snakemake@params[["cols_nan_threshold"]]
drop_zero_variance_columns <- snakemake@params[["cols_var_threshold"]]
rows_nan_threshold <- snakemake@params[["rows_nan_threshold"]]
days_before_threshold <- snakemake@params[["days_before_threshold"]]
days_after_threshold <- snakemake@params[["days_after_threshold"]]


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
  clean_features <- clean_features %>% select(-day_type)
}

write.csv(clean_features, snakemake@output[[1]], row.names = FALSE)
