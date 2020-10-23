source("renv/activate.R")
library("tibble")
library("dplyr", warn.conflicts = F)
library("tidyr")
library("tibble")
options(scipen=999)

# Using mostly indeixng instead of tidyr because is faster
resampled_episodes <- read.csv(snakemake@input[[1]]) 
resampled_episodes["n_resamples"] <- 1 + (resampled_episodes["end_timestamp"] - resampled_episodes["start_timestamp"]) %/% 60001
resampled_episodes <- resampled_episodes %>% uncount(n_resamples, .id = "nrow")

resampled_episodes["nrow"] <- (resampled_episodes["nrow"] - 1) * 60000 
resampled_episodes["start_timestamp"] <- resampled_episodes["start_timestamp"] + resampled_episodes["nrow"]
# Use +59999 because each resampled minute should not overlap with each other
resampled_episodes["end_timestamp"] <- pmin(resampled_episodes["start_timestamp"] + 59999, resampled_episodes["end_timestamp"])
resampled_episodes <- resampled_episodes %>% select(-nrow)
resampled_episodes <- resampled_episodes %>% uncount(2, .id = "end_flag") 

resampled_episodes <- resampled_episodes %>% add_column(timestamp = NA_real_)
if(nrow(resampled_episodes) > 0){
  resampled_episodes[resampled_episodes$end_flag ==1, "timestamp"] = resampled_episodes[resampled_episodes$end_flag ==1, "start_timestamp"]
  resampled_episodes[resampled_episodes$end_flag ==2, "timestamp"] = resampled_episodes[resampled_episodes$end_flag ==2, "end_timestamp"]
}
resampled_episodes <- resampled_episodes %>% select(-end_flag)

write.csv(resampled_episodes, snakemake@output[[1]], row.names = FALSE)
