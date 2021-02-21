source("renv/activate.R")

library(RMariaDB)
library(stringr)
library(purrr)
library(readr)
library("dplyr", warn.conflicts = F)

config <- snakemake@params[["config"]]
group <- config$SOURCE$DATABASE_GROUP
timezone <- config$SOURCE$TIMEZONE
phone_device_id_column = config$PHONE_SECTION$DEVICE_ID_COLUMN
fitbit_device_id_column = config$FITBIT_SECTION$DEVICE_ID_COLUMN
add_phone_section = config$PHONE_SECTION$ADD
add_fitbit_section = config$FITBIT_SECTION$ADD
add_empatica_section = config$EMPATICA_SECTION$ADD
phone_ignored = config$PHONE_SECTION$IGNORED_DEVICE_IDS
fitbit_ignored = config$FITBIT_SECTION$IGNORED_DEVICE_IDS

rmysql.settingsfile <- "./.env"

if(config$SOURCE$TYPE == "AWARE_DEVICE_TABLE"){
  database <- dbConnect(MariaDB(), default.file = rmysql.settingsfile, group = group)
  if(config$FITBIT_SECTION$ADD == TRUE){
    query <- paste("SELECT",phone_device_id_column, ",",fitbit_device_id_column," as _temp_fitbit_id, brand, label, timestamp FROM aware_device order by timestamp asc")
    fitbit_device_id_column <- "_temp_fitbit_id"
  }
  else 
    query <- paste("SELECT ",phone_device_id_column,", brand, label, timestamp FROM aware_device order by timestamp asc")
  participants <- dbGetQuery(database, query)
  dbDisconnect(database)
  participants <- participants %>% 
    mutate(pid = if_else(row_number()<10, paste0("p","0",row_number()), paste0("p", row_number())),
           platform = if_else(brand == "iPhone", "ios", "android"), brand = NULL,
           label = iconv(if_else(label == "", "EMPTY_LABEL", label), from = "UTF-8", to = "UTF-8", sub=''),
           start_date = format(as.POSIXct(timestamp / 1000, origin = "1970-01-01", tz = timezone), "%Y-%m-%d"),
           end_date = format(Sys.Date(), "%Y-%m-%d"),
           !!phone_device_id_column := if_else(!!rlang::sym(phone_device_id_column) %in% phone_ignored, NA_character_, !!rlang::sym(phone_device_id_column)),
           !!fitbit_device_id_column := if_else(!!rlang::sym(fitbit_device_id_column) %in% fitbit_ignored, NA_character_, !!rlang::sym(fitbit_device_id_column)))

} else if(config$SOURCE$TYPE == "CSV_FILE"){
  participants <- read_csv(config$SOURCE$CSV_FILE_PATH, col_types=cols_only(device_id="c",pid="c",label="c",platform="c",
                            start_date=col_date(format = "%Y-%m-%d"),end_date=col_date(format = "%Y-%m-%d"),fitbit_id="c")) %>% 
                            mutate(start_date = as.character(start_date), end_date = as.character(end_date)) # we read as date to validate format
  participants <- participants %>% 
  mutate(!!phone_device_id_column := str_replace(!!rlang::sym(phone_device_id_column), ";",","),
         platform = str_replace(platform, ";",","),
         !!phone_device_id_column := if_else(!!rlang::sym(phone_device_id_column) %in% phone_ignored, NA_character_, !!rlang::sym(phone_device_id_column)),
         !!fitbit_device_id_column := if_else(!!rlang::sym(fitbit_device_id_column) %in% fitbit_ignored, NA_character_, !!rlang::sym(fitbit_device_id_column)))
}

dir.create(file.path("./data/external/participant_files/"))

participants %>%
  pwalk(function(add_phone_section, add_fitbit_section, phone_device_id_column, fitbit_device_id_column, ...) {
    empty_phone <- c("PHONE:", "  DEVICE_IDS:", "  PLATFORMS:","  LABEL:", "  START_DATE:", "  END_DATE:")
    empty_fitbit <- c("FITBIT:", "  DEVICE_IDS:", "  LABEL:", "  START_DATE:", "  END_DATE:")
    empty_empatica <- c("EMPATICA:", "  LABEL:", "  START_DATE:", "  END_DATE:")
    row <- tibble(...)
    lines <- c()
    start_date = if_else(is.na(row$start_date), "", row$start_date)
    end_date = if_else(is.na(row$end_date), "", row$end_date)

    if(add_phone_section == TRUE && !is.na(row[phone_device_id_column])){
      lines <- append(lines, c("PHONE:", paste0("  DEVICE_IDS: [",row[phone_device_id_column],"]"), paste0("  PLATFORMS: [",row$platform,"]"),
                               paste("  LABEL:",row$label), paste("  START_DATE:", start_date), paste("  END_DATE:", end_date)))
    }else
      lines <- append(lines, empty_phone)
    
    if(add_fitbit_section == TRUE && !is.na(row[fitbit_device_id_column])){
      lines <- append(lines, c("FITBIT:", paste0("  DEVICE_IDS: [",row[fitbit_device_id_column],"]"),
                               paste("  LABEL:",row$label), paste("  START_DATE:", start_date), paste("  END_DATE:", end_date)))
    } else
      lines <- append(lines, empty_fitbit)

    if(add_empatica_section == TRUE){
      lines <- append(lines, c("EMPATICA:",
                               paste("  LABEL:",row$label), paste("  START_DATE:", start_date), paste("  END_DATE:", end_date)))
    } else
      lines <- append(lines, empty_empatica)
    
    file_connection <- file(paste0("./data/external/participant_files/", row$pid, ".yaml"))
    writeLines(lines, file_connection)
    close(file_connection)

  }, add_phone_section, add_fitbit_section, phone_device_id_column, fitbit_device_id_column)

file_lines <-readLines("./config.yaml")
for (i in 1:length(file_lines)){
  if(startsWith(file_lines[i], "PIDS:")){
    file_lines[i] <- paste0("PIDS: [", paste(participants$pid, collapse = ", "), "]")
  }
}
writeLines(file_lines, con = "./config.yaml") 