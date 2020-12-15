source("renv/activate.R")

library("tidyr")
library("dplyr", warn.conflicts = F)

empatica_files <- snakemake@input[["input_files"]]
empatica_data <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), c("timestamp"))


for(file in empatica_files){
    data <- read.csv(file)
    if(! ("timestamp" %in% colnames(data)))
        stop(paste("This file does not have a timestamp column, something might have gone wrong while unzipping it:", file))
    empatica_data <- merge(empatica_data, data, all = TRUE)
}

write.csv(empatica_data, snakemake@output[[1]], row.names = FALSE)