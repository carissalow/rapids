source("renv/activate.R")
# Load Ian Barnett's code. Taken from https://scholar.harvard.edu/ibarnett/software/gpsmobility
file.sources = list.files(c("src/features/location_barnett"), pattern="*.R$", full.names=TRUE, ignore.case=TRUE)
sapply(file.sources,source,.GlobalEnv)

library(dplyr)

write_empty_file <- function(file_path, requested_features){
  write.csv(data.frame(local_date= character(), 
                        location_barnett_hometime= numeric(), 
                        location_barnett_disttravelled= numeric(), 
                        location_barnett_rog= numeric(), 
                        location_barnett_maxdiam= numeric(), 
                        location_barnett_maxhomedist= numeric(), 
                        location_barnett_siglocsvisited= numeric(), 
                        location_barnett_avgflightlen= numeric(), 
                        location_barnett_stdflightlen= numeric(), 
                        location_barnett_avgflightdur= numeric(), 
                        location_barnett_stdflightdur= numeric(), 
                        location_barnett_probpause= numeric(), 
                        location_barnett_siglocentropy= numeric(), 
                        location_barnett_minsmissing= numeric(), 
                        location_barnett_circdnrtn= numeric(), 
                        location_barnett_wkenddayrtn= numeric(),
                        minutes_data_used= numeric()
                      ) %>% select(requested_features), file_path, row.names = F)
}

location <- read.csv(snakemake@input[["locations"]], stringsAsFactors = F) 
# The choice between RESAMPLE_FUSED and the original location data happens at the rule level in the function
# optional_location_input in features.snakefile
locations_to_use <- snakemake@params[["locations_to_use"]]
accuracy_limit <- snakemake@params[["accuracy_limit"]]
timezone <- snakemake@params[["timezone"]]
minutes_data_used <- snakemake@params[["minutes_data_used"]]
requested_features <- intersect(unlist(snakemake@params["features"], use.names = F), 
                                c("hometime","disttravelled","rog","maxdiam","maxhomedist","siglocsvisited","avgflightlen","stdflightlen","avgflightdur","stdflightdur","probpause","siglocentropy","minsmissing","circdnrtn","wkenddayrtn"))
requested_features <- c("local_date", paste("location_barnett", requested_features, sep = "_"))
if(minutes_data_used)
  requested_features <- c(requested_features, "minutes_data_used")

if(!locations_to_use %in% c("ALL_EXCEPT_FUSED", "RESAMPLE_FUSED", "ALL")){
  print("Unkown filter, provide one of the following three: ALL, ALL_EXCEPT_FUSED, or RESAMPLE_FUSED")
  quit(save = "no", status = 1, runLast = FALSE)
}

 # excludes fused and resample
if(locations_to_use == "ALL_EXCEPT_FUSED")
  location <- location %>% filter(provider == "gps")

# Remove 0,0 location coordinates
location <- location %>% filter(double_latitude != 0 & double_longitude != 0)

# Excludes datasets with less than 24 hours of data
if(max(location$timestamp) - min(location$timestamp) < 86400000)
  location <- head(location, 0)

if (nrow(location) > 1){

    # Count how many minutes of data we use to get location features
    # Some minutes have multiple fused  rows
    location_minutes_used <- location %>% 
      group_by(local_date, local_hour) %>% 
      summarise(n_minutes = n_distinct(local_minute)) %>% 
      group_by(local_date) %>% 
      summarise(minutes_data_used = sum(n_minutes)) %>% 
      select(local_date, minutes_data_used)

    location <- location %>%
      select(timestamp, latitude = double_latitude, longitude = double_longitude, altitude = double_altitude, accuracy)

    outputMobility <- MobilityFeatures(location, ACCURACY_LIM = accuracy_limit, tz = timezone)

    if(is.null(outputMobility)){
      write_empty_file(snakemake@output[[1]], requested_features)
    } else{
      # Copy index (dates) as a column 
      features <- cbind(rownames(outputMobility$featavg), outputMobility$featavg)
      features <- as.data.frame(features)
      features[-1] <- lapply(lapply(features[-1], as.character), as.numeric)
      colnames(features)=c("local_date",tolower(paste("location_barnett", colnames(outputMobility$featavg), sep = "_")))
      # Add the minute count column
      features <- left_join(features, location_minutes_used, by = "local_date")
      write.csv(features %>% select(requested_features), snakemake@output[[1]], row.names = F)
    }
    
} else {
    write_empty_file(snakemake@output[[1]], requested_features)    
}
