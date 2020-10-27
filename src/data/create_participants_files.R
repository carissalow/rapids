source("renv/activate.R")

library(RMySQL)
library(stringr)
library(purrr)
library(readr)
library("dplyr", warn.conflicts = F)

config <- snakemake@params[["config"]]
group <- config$SOURCE$DATABASE_GROUP
timezone <- config$SOURCE$TIMEZONE
phone_device_id_column = config$PHONE_SECTION$DEVICE_ID_COLUMN
fitbit_device_id_column = config$FITBIT_SECTION$DEVICE_ID_COLUMN
add_fitbit_section = config$PHONE_SECTION$ADD
add_phone_section = config$FITBIT_SECTION$ADD
phone_ignored = config$PHONE_SECTION$IGNORED_DEVICE_IDS
fitbit_ignored = config$FITBIT_SECTION$IGNORED_DEVICE_IDS

rmysql.settingsfile <- "./.env"

if(config$SOURCE$TYPE == "AWARE_DEVICE_TABLE"){
  database <- dbConnect(MySQL(), default.file = rmysql.settingsfile, group = group)
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
                            start_date=col_date(format = "%Y-%m-%d"),end_date=col_date(format = "%Y-%m-%d"),fitbit_id="c"))
  participants <- participants %>% 
  mutate(!!phone_device_id_column := str_replace(!!rlang::sym(phone_device_id_column), ";",","),
         platform = str_replace(platform, ";",","),
         !!phone_device_id_column := if_else(!!rlang::sym(phone_device_id_column) %in% phone_ignored, NA_character_, !!rlang::sym(phone_device_id_column)),
         !!fitbit_device_id_column := if_else(!!rlang::sym(fitbit_device_id_column) %in% fitbit_ignored, NA_character_, !!rlang::sym(fitbit_device_id_column)))
}

participants %>%
  pwalk(function(add_phone_section, add_fitbit_section, phone_device_id_column, fitbit_device_id_column, ...) {
    empty_phone <- c("PHONE:", "  DEVICE_IDS:", "  PLATFORMS:","  LABEL:", "  START_DATE:", "  END_DATE:")
    empty_fitbit <- c("FITBIT:", "  DEVICE_IDS:", "  LABEL:", "  START_DATE:", "  END_DATE:")
    row <- tibble(...)
    lines <- c()

    if(add_phone_section == TRUE && !is.na(row[phone_device_id_column])){
      lines <- append(lines, c("PHONE:", paste0("  DEVICE_IDS: [",row[phone_device_id_column],"]"), paste0("  PLATFORMS: [",row$platform,"]"),
                               paste("  LABEL:",row$label), paste("  START_DATE:", row$start_date), paste("  END_DATE:", row$end_date)))
    }else
      lines <- append(lines, empty_phone)
    
    if(add_fitbit_section == TRUE && !is.na(row[fitbit_device_id_column])){
      lines <- append(lines, c("FITBIT:", paste0("  DEVICE_IDS: [",row[fitbit_device_id_column],"]"),
                               paste("  LABEL:",row$label), paste("  START_DATE:", row$start_date), paste("  END_DATE:", row$end_date)))
    } else
      lines <- append(lines, empty_fitbit)
    
    file_connection <- file(paste0("./data/external/participant_files/", row$pid, ".yaml"))
    writeLines(lines, file_connection)
    close(file_connection)

  }, add_phone_section, add_fitbit_section, phone_device_id_column, fitbit_device_id_column)