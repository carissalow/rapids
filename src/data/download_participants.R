source("packrat/init.R")

library(RMySQL)

group <- snakemake@params[["group"]]
ignored_device_ids <- snakemake@params[["ignored_device_ids"]]
rmysql.settingsfile <- "./.env"

stopDB <- dbConnect(MySQL(), default.file = rmysql.settingsfile, group = group)
query <- paste0("SELECT device_id, brand FROM aware_device order by timestamp asc")
participants <- dbGetQuery(stopDB, query)
pids <- c()

for(id in 1:nrow(participants)){
    device_id <- participants$device_id[[id]]
    brand <- ifelse(participants$brand[[id]] == "iPhone", "ios", "android")
    if(!(device_id %in% ignored_device_ids)){
        pid <- paste0("p", ifelse(id < 10, paste0("0", id), id))
        pids <- append(pids, pid)
        file_connection <- file(paste0("./data/external/", pid))
        writeLines(c(device_id, brand), file_connection)
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