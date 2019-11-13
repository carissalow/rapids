source("packrat/init.R")

library(dplyr)

write_empty_file <- function(file_path){
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
                      ), file_path, row.names = F)
}

# Load Ian Barnett's code. Taken from https://scholar.harvard.edu/ibarnett/software/gpsmobility
file.sources = list.files(c("src/features/location_barnett"), pattern="*.R$", full.names=TRUE, ignore.case=TRUE)
sapply(file.sources,source,.GlobalEnv)

accuracy_limit <- snakemake@params[["accuracy_limit"]]
timezone <- snakemake@params[["timezone"]]

location <- read.csv(snakemake@input[[1]], stringsAsFactors = F) %>%
  select(timestamp, latitude = double_latitude, longitude = double_longitude, altitude = double_altitude, accuracy)

if (nrow(location) > 1){
    features <- MobilityFeatures(location, ACCURACY_LIM = accuracy_limit, tz = timezone)
    if(is.null(features)){
      write_empty_file(snakemake@output[[1]])
    } else{
      # Copy index (dates) as a column 
      outmatrix <- cbind(rownames(features$featavg), features$featavg)
      outmatrix <- as.data.frame(outmatrix)
      outmatrix[-1] <- lapply(lapply(outmatrix[-1], as.character), as.numeric)
      colnames(outmatrix)=c("local_date",tolower(colnames(features$featavg)))
      write.csv(outmatrix,snakemake@output[[1]], row.names = F)
    }
    
} else {
    write_empty_file(snakemake@output[[1]])    
}