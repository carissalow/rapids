source("renv/activate.R")
library("dplyr", warn.conflicts = F)

unify_ios_conversation <- function(conversation){
    if(nrow(conversation) > 0){
        duration_check <- conversation %>% 
            select(double_convo_start, double_convo_end) %>% 
            mutate(start_is_seconds = double_convo_start <= 9999999999,
                end_is_seconds = double_convo_end <= 9999999999) # Values smaller than 9999999999 are in seconds instead of milliseconds
        start_end_in_seconds = sum(duration_check$start_is_seconds) + sum(duration_check$end_is_seconds)

        if(start_end_in_seconds > 0) # convert seconds to milliseconds
            conversation <- conversation %>% mutate(double_convo_start = double_convo_start * 1000, double_convo_end = double_convo_end * 1000)
    }
    return(conversation) 
}
main <- function(data, stream_parameters){
    return(unify_ios_conversation(data))
}