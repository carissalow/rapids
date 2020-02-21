source("packrat/init.R")

library(dplyr)

write_empty_file <- function(file_path, metrics_to_include){
  write.csv(data.frame(local_date= character(), 
                        hometime= numeric(), 
                        disttravelled= numeric(), 
                        rog= numeric(), 
                        maxdiam= numeric(), 
                        maxhomedist= numeric(), 
                        siglocsvisited= numeric(), 
                        avgflightlen= numeric(), 
                        stdflightlen= numeric(), 
                        avgflightdur= numeric(), 
                        stdflightdur= numeric(), 
                        probpause= numeric(), 
                        siglocentropy= numeric(), 
                        minsmissing= numeric(), 
                        circdnrtn= numeric(), 
                        wkenddayrtn= numeric()
                      ) %>% select(metrics_to_include), file_path, row.names = F)
}

# Load Ian Barnett's code. Taken from https://scholar.harvard.edu/ibarnett/software/gpsmobility
file.sources = list.files(c("src/features/location_barnett"), pattern="*.R$", full.names=TRUE, ignore.case=TRUE)
sapply(file.sources,source,.GlobalEnv)

locations_to_use <- snakemake@params[["locations_to_use"]]
accuracy_limit <- snakemake@params[["accuracy_limit"]]
timezone <- snakemake@params[["timezone"]]
metrics_to_include <- intersect(unlist(snakemake@params["metrics"], use.names = F), 
                                c("hometime","disttravelled","rog","maxdiam","maxhomedist","siglocsvisited","avgflightlen","stdflightlen","avgflightdur","stdflightdur","probpause","siglocentropy","minsmissing","circdnrtn","wkenddayrtn"))

# By deafult we use all raw locations: fused without resampling and not fused (gps, network)
location <- read.csv(snakemake@input[["raw"]], stringsAsFactors = F) %>%
  select(timestamp, latitude = double_latitude, longitude = double_longitude, altitude = double_altitude, accuracy)

if(locations_to_use == "ALL_EXCEPT_FUSED"){
  location <- location %>% filter(provider != "fused")
} else if (locations_to_use == "RESAMPLE_FUSED"){
  location <- read.csv(snakemake@input[["fused"]], stringsAsFactors = F) %>%
    select(timestamp, latitude = double_latitude, longitude = double_longitude, altitude = double_altitude, accuracy)
} else if (locations_to_use != "ALL"){
  print("Unkown filter, provide one of the following three: ALL, ALL_EXCEPT_FUSED, or RESAMPLE_FUSED")
  quit(save = "no", status = 1, runLast = FALSE)
}

if (nrow(location) > 1){
    features <- MobilityFeatures(location, ACCURACY_LIM = accuracy_limit, tz = timezone)
    if(is.null(features)){
      write_empty_file(snakemake@output[[1]], metrics_to_include)
    } else{
      # Copy index (dates) as a column 
      outmatrix <- cbind(rownames(features$featavg), features$featavg)
      outmatrix <- as.data.frame(outmatrix)
      outmatrix[-1] <- lapply(lapply(outmatrix[-1], as.character), as.numeric)
      colnames(outmatrix)=c("local_date",tolower(colnames(features$featavg)))
      write.csv(outmatrix %>% select(metrics_to_include), snakemake@output[[1]], row.names = F)
    }
    
} else {
    write_empty_file(snakemake@output[[1]], metrics_to_include)    
}