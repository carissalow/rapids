source("renv/activate.R")
library("dplyr", warn.conflicts = F)


unify_ios_screen <- function(ios_screen){
    # In Android we only process UNLOCK to OFF episodes. In iOS we only process UNLOCK to LOCKED episodes,
    # thus, we replace LOCKED with OFF episodes (2 to 0) so we can use Android's code for iOS
    ios_screen <- ios_screen %>% 
        # only keep consecutive pairs of 3,2 events
        filter( (screen_status == 3 & lead(screen_status) == 2) | (screen_status == 2 & lag(screen_status) == 3) ) %>%
        mutate(screen_status = replace(screen_status, screen_status == 2, 0))
    return(ios_screen)
}

main <- function(data, stream_parameters){
    return(unify_ios_screen(data))
}