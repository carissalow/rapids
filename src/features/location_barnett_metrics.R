library(dplyr)

# Load Ian Barnett's code. Taken from https://scholar.harvard.edu/ibarnett/software/gpsmobility
file.sources = list.files(c("src/features/location_barnett"), pattern="*.R$", full.names=TRUE, ignore.case=TRUE)
sapply(file.sources,source,.GlobalEnv)

accuracy_limit <- snakemake@params[["accuracy_limit"]]
timezone <- snakemake@params[["timezone"]]

location <- read.csv(snakemake@input[[1]], stringsAsFactors = F) %>%
  select(timestamp, latitude = double_latitude, longitude = double_longitude, altitude = double_altitude, accuracy)

if (nrow(location) > 0){
    features <- MobilityFeatures(location, ACCURACY_LIM = accuracy_limit, tz = timezone)

    # Copy index (dates) as a column 
    outmatrix = cbind(rownames(features$featavg), features$featavg)
    colnames(outmatrix)=c("local_date",tolower(colnames(features$featavg)))
    write.csv(outmatrix,snakemake@output[[1]], row.names = F)
} else {
    write.csv(data.frame(),snakemake@output[[1]])
}