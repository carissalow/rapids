source("renv/activate.R")
library("dplyr")

if(!is.null(snakemake@input[["visible_access_points"]]) && is.null(snakemake@input[["connected_access_points"]])){
  wifi_data <- read.csv(snakemake@input[["visible_access_points"]], stringsAsFactors = FALSE)
  wifi_data <- wifi_data %>% mutate(connected = 0)
} else if(is.null(snakemake@input[["visible_access_points"]]) && !is.null(snakemake@input[["connected_access_points"]])){
  wifi_data <- read.csv(snakemake@input[["connected_access_points"]], stringsAsFactors = FALSE)
  wifi_data <- wifi_data %>% mutate(connected = 1)
} else if(!is.null(snakemake@input[["visible_access_points"]]) && !is.null(snakemake@input[["connected_access_points"]])){
  visible_access_points <- read.csv(snakemake@input[["visible_access_points"]], stringsAsFactors = FALSE)
  visible_access_points <- visible_access_points %>% mutate(connected = 0)
  connected_access_points <- read.csv(snakemake@input[["connected_access_points"]], stringsAsFactors = FALSE)
  connected_access_points <- connected_access_points %>% mutate(connected = 1)
  wifi_data <- bind_rows(visible_access_points, connected_access_points) %>% arrange(timestamp)
}

write.csv(wifi_data, snakemake@output[[1]], row.names = FALSE)