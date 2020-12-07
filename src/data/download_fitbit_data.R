source("renv/activate.R")
library(RMariaDB)
library("dplyr", warn.conflicts = F)
library(readr)
library(stringr)
library(yaml)


participant_file <- snakemake@input[["participant_file"]]
input_file <- snakemake@input[["input_file"]]
data_configuration <- snakemake@params[["data_configuration"]]
source <- data_configuration$SOURCE
sensor <- snakemake@params[["sensor"]]
table <- snakemake@params[["table"]]
sensor_file <- snakemake@output[[1]]

participant <- read_yaml(participant_file)
if(! "FITBIT" %in% names(participant)){
  stop(paste("The following participant file does not have a FITBIT section, create one manually or automatically (see the docs):", participant_file))
}
device_ids <- participant$FITBIT$DEVICE_IDS
unified_device_id <- tail(device_ids, 1)
# As opposed to phone data, we dont' filter by date here because data can still be in JSON format, we need to parse it first

if(source$TYPE == "DATABASE"){
  dbEngine <- dbConnect(MariaDB(), default.file = "./.env", group = source$DATABASE_GROUP)
  query <- paste0("SELECT * FROM ", table, " WHERE ",source$DEVICE_ID_COLUMN," IN ('", paste0(device_ids, collapse = "','"), "')")
  sensor_data <- dbGetQuery(dbEngine, query)
  dbDisconnect(dbEngine)
} else if(source$TYPE == "FILES"){
  sensor_data <- read_csv_chunked(input_file, callback = DataFrameCallback$new(function(x, pos) subset(x,x[[source$DEVICE_ID_COLUMN]] %in% device_ids)), progress = T, chunk_size = 50000)
  if(is.null(sensor_data)) # emtpy file
    sensor_data <- read.csv(input_file)
}

sensor_data <- sensor_data %>%
  rename(device_id = source$DEVICE_ID_COLUMN) %>% 
  mutate(device_id = unified_device_id) # Unify device_id

if("HIDDEN" %in% names(data_configuration) && data_configuration$HIDDEN$SINGLE_FITBIT_TABLE == TRUE) # For MoSHI use, we didn't split fitbit sensors into different tables
  sensor_data <- sensor_data %>% filter(fitbit_data_type == str_split(sensor, "_", simplify = TRUE)[[2]])

# Droping duplicates on all columns except for _id or id
sensor_data <- sensor_data %>% distinct(!!!syms(setdiff(names(sensor_data), c("_id", "id"))))

write_csv(sensor_data, sensor_file)