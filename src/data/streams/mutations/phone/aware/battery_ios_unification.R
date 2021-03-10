source("renv/activate.R")
library("dplyr", warn.conflicts = F)

unify_ios_battery <- function(ios_battery){
    # We only need to unify battery data for iOS client V1. V2 does it out-of-the-box
    # V1 will not have rows where battery_status is equal to 4
    if(nrow(ios_battery %>% filter(battery_status == 4)) == 0)
        ios_battery <- ios_battery %>%
            mutate(battery_status = replace(battery_status, battery_status == 3, 5),
                battery_status = replace(battery_status, battery_status == 1, 3))
    return(ios_battery)
}

main <- function(data, stream_parameters){
    return(unify_ios_battery(data))
}