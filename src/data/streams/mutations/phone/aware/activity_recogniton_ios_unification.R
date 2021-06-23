source("renv/activate.R")
library("dplyr", warn.conflicts = F)
library(stringr)

clean_ios_activity_column <- function(ios_gar){
    ios_gar <- ios_gar %>%
        mutate(activities = str_replace_all(activities, pattern = '("|\\[|\\])', replacement = ""))

    existent_multiple_activities <- ios_gar %>%
        filter(str_detect(activities, ",")) %>% 
        group_by(activities) %>%
        summarise(mutiple_activities = unique(activities), .groups = "drop_last") %>% 
        pull(mutiple_activities)

    known_multiple_activities <- c("stationary,automotive")
    unkown_multiple_actvities <- setdiff(existent_multiple_activities, known_multiple_activities)
    if(length(unkown_multiple_actvities) > 0){
        stop(paste0("There are unkwown combinations of ios activities, you need to implement the decision of the ones to keep: ", unkown_multiple_actvities))
    }

    ios_gar <- ios_gar %>%
        mutate(activities = str_replace_all(activities, pattern = "stationary,automotive", replacement = "automotive"))
    
    return(ios_gar)
}

unify_ios_activity_recognition <- function(ios_gar){
    # discard rows where activities column is blank
    ios_gar <- ios_gar %>% filter(!is.na(activities) & activities != "" )
    # clean "activities" column of ios_gar
    ios_gar <- clean_ios_activity_column(ios_gar)

    # make it compatible with android version: generate "activity_name" and "activity_type" columns
    ios_gar  <-  ios_gar %>% 
        mutate(activity_name = case_when(activities == "automotive" ~ "in_vehicle",
                                         activities == "cycling" ~ "on_bicycle",
                                         activities == "walking" ~ "walking",
                                         activities == "running" ~ "running",
                                         activities == "stationary" ~ "still",
                                         activities == "unknown" ~ "unknown"),
               activity_type = case_when(activities == "automotive" ~ 0,
                                         activities == "cycling" ~ 1,
                                         activities == "walking" ~ 7,
                                         activities == "running" ~ 8,
                                         activities == "stationary" ~ 3,
                                         activities == "unknown" ~ 4),
                confidence = case_when(confidence == 0 ~ 0,
                                       confidence == 1 ~ 50,
                                       confidence == 2 ~ 100)
                                       ) %>% 
        select(-activities)
    
    return(ios_gar)
}

main <- function(data, stream_parameters){
    return(unify_ios_activity_recognition(data))
}