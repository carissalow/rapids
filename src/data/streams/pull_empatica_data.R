source("renv/activate.R")

library(yaml)
library(dplyr)
library(readr)
# we use reticulate but only load it if we are going to use it to minimize the case when old RAPIDS deployments need to update ther renv
mutate_data <- function(scripts, data){
  for(script in scripts){
    if(grepl("\\.(R)$", script)){
      myEnv <- new.env()    
      source(script, local=myEnv)
      attach(myEnv, name="sourced_scripts_rapids")
      if(exists("main", myEnv)){
        message(paste("Applying mutation script", script))
        data <- main(data)
      } else{
        stop(paste0("The following mutation script does not have main function: ", script))
      }
      # rm(list = ls(envir = myEnv), envir = myEnv, inherits = FALSE)
      detach("sourced_scripts_rapids")
    } else{ # python
      library(reticulate)
      module <- gsub(pattern = "\\.py$", "", basename(script))
      script_functions <- import_from_path(module, path = dirname(script))
      if(py_has_attr(script_functions, "main")){
        message(paste("Applying mutation script", script))
        data <- script_functions$main(data)
      } else{
        stop(paste0("The following mutation script does not have a main function: ", script))
      }
    }
  }

  return(data)
}

rename_columns <- function(name_maps, data){
  for(name in names(name_maps))
    data <- data %>% rename(!!tolower(name) := name_maps[[name]])
  return(data)
}

validate_expected_columns_mapping <- function(schema, rapids_schema, sensor, rapids_schema_file, stream_format){
  columns <- names(schema[[sensor]][["COLUMN_MAPPINGS"]])
  columns <- columns[(columns != "FLAG_AS_EXTRA")]
  rapids_columns <- rapids_schema[[sensor]]

  if(is.null(rapids_columns))
    stop(paste(sensor, " columns are not listed in RAPIDS' column specification. If you are adding support for a new phone sensor, add any mandatory columns in ", rapids_schema_file))
  if(length(setdiff(rapids_columns, columns)) > 0)
    stop(paste(sensor," mappings are missing one or more mandatory columns. The missing column mappings are for ", paste(setdiff(rapids_columns, columns), collapse=","),"in", stream_format, " (the mappings are case sensitive)"))
  if(length(setdiff(columns, rapids_columns)) > 0)
    stop(paste(sensor," mappings have one or more columns than required, add them as FLAG_AS_EXTRA instead. The extra column mappings are for ", paste(setdiff(columns, rapids_columns), collapse=","),"in", stream_format, " (the mappings are case sensitive)"))
}

load_container_script <- function(stream_container){
  language <- if_else(endsWith(tolower(stream_container), "py"), "python", "r")
  if(language == "python"){
    library(reticulate)
    container <- import_from_path(gsub(pattern = "\\.py$", "", basename(stream_container)), path = dirname(stream_container))
    if(!py_has_attr(container, "pull_data"))
      stop(paste0("The following container.py script does not have a pull_data function: ", stream_container))
    return(container$pull_data)
  } else if(language == "r"){
    source(stream_container)
    if(!exists("pull_data"))
      stop(paste0("The following container.R script does not have a pull_data function: ", stream_container))
    return(pull_data)
  }
}

pull_empatica_data_main <- function(){
  participant_file <- snakemake@input[["participant_file"]]
  stream_format <- snakemake@input[["stream_format"]]
  rapids_schema_file <- snakemake@input[["rapids_schema_file"]]
  stream_container <- snakemake@input[["stream_container"]]
  data_configuration <- snakemake@params[["data_configuration"]]
  pid <- snakemake@params[["pid"]]
  table <- snakemake@params[["table"]]
  sensor <- toupper(snakemake@params[["sensor"]])
  output_data_file <- snakemake@output[[1]]


  participant_data <- read_yaml(participant_file)
  stream_schema <- read_yaml(stream_format)
  rapids_schema <- read_yaml(rapids_schema_file)
  devices <- participant_data$EMPATICA$DEVICE_IDS
  if(length(devices) == 0)
    devices <- c(pid)
  validate_expected_columns_mapping(stream_schema, rapids_schema, sensor, rapids_schema_file, stream_format)
  expected_columns <- tolower(names(stream_schema[[sensor]][["COLUMN_MAPPINGS"]]))
  expected_columns <- expected_columns[(expected_columns != "flag_extra")]
  participant_data <- setNames(data.frame(matrix(ncol = length(expected_columns), nrow = 0)), expected_columns)

  pull_data_container <- load_container_script(stream_container)

  for(idx in seq_along(devices)){ #TODO remove length    
    device <- devices[idx]
    message(paste0("\nProcessing ", sensor, " for ", device))

    columns_to_download <- stream_schema[[sensor]][["COLUMN_MAPPINGS"]]
    columns_to_download <- columns_to_download[(columns_to_download != "FLAG_TO_MUTATE")]
    data <- pull_data_container(data_configuration, device, sensor, columns_to_download)

    columns_to_rename <- stream_schema[[sensor]][["COLUMN_MAPPINGS"]]
    columns_to_rename <- (columns_to_rename[(columns_to_rename != "FLAG_TO_MUTATE" & names(columns_to_rename) != "FLAG_AS_EXTRA")])
    renamed_data <- rename_columns(columns_to_rename, data)

    mutation_scripts <- stream_schema[[sensor]][["MUTATION_SCRIPTS"]]
    mutated_data <- mutate_data(mutation_scripts, renamed_data)

    if(length(setdiff(expected_columns, colnames(mutated_data))) > 0)
      stop(paste("The mutated data for", device, "is missing these columns expected by RAPIDS: [", paste(setdiff(expected_columns, colnames(mutated_data)), collapse=","),"]. One ore more mutation scripts in [", sensor,"][",toupper(device_os), "]","[MUTATION_SCRIPTS] are removing or not adding these columns"))
    participant_data <- rbind(participant_data, mutated_data)
      
  }
  
  write_csv(participant_data, output_data_file)
}

pull_empatica_data_main()