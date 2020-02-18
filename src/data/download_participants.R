source("packrat/init.R")

library(RMySQL)

group <- snakemake@params[["group"]]
ignored_device_ids <- snakemake@params[["ignored_device_ids"]]
timezone <- snakemake@params[["timezone"]]
rmysql.settingsfile <- "./.env"

stopDB <- dbConnect(MySQL(), default.file = rmysql.settingsfile, group = group)
query <- "SELECT device_id, brand, label, timestamp FROM aware_device order by timestamp asc"
participants <- dbGetQuery(stopDB, query)
pids <- c()

end_date <- format(Sys.Date(), "%Y/%m/%d")

for(id in 1:nrow(participants)){
    device_id <- participants$device_id[[id]]
    brand <- ifelse(participants$brand[[id]] == "iPhone", "ios", "android")
    label <- ifelse(participants$label[[id]] == "", "EMPTY_LABEL", participants$label[[id]])
    start_date <- format(as.POSIXct(participants$timestamp[[id]] / 1000, origin = "1970-01-01", tz = timezone), "%Y/%m/%d")
    if(!(device_id %in% ignored_device_ids)){
        pid <- paste0("p", ifelse(id < 10, paste0("0", id), id))
        pids <- append(pids, pid)
        file_connection <- file(paste0("./data/external/", pid))
        writeLines(c(device_id, brand, label, paste0(start_date, ",", end_date)), file_connection)
        close(file_connection)
    }
}

file_lines <-readLines("./config.yaml")
for (i in 1:length(file_lines)){
  if(startsWith(file_lines[i], "PIDS:")){
    file_lines[i] <- paste0("PIDS: [", paste(pids, collapse = ", "), "]")
  }
}
writeLines(file_lines, con = "./config.yaml")