source("renv/activate.R")

library(tidyr)
library("dplyr", warn.conflicts = F)
library(stringr)
library("rvest")

get_genre <- function(apps){
  urls = paste0("https://play.google.com/store/apps/details?id=", apps)
  destfiles = paste0(apps,".html")
  genres = vector("character", length(apps))

  for(i in seq_along(urls)){
    try_download <- try(download.file(urls[i], destfiles[i], quiet=TRUE), silent = T)
    if(is(try_download,"try-error") || (read_html(destfiles[i]) %>% html_nodes("title") %>% html_text()) == "Not Found"){
      genres[i] <- "unknown"
    }
    else{
      genres[i] <- read_html(destfiles[i]) %>% 
        html_nodes(xpath = '//*[@itemprop="genre"]') %>%
        html_text()
      genres[i]  <- tolower(str_remove_all(genres[i], "[\\s&]+")) # removes white spaces or ampersands
    }
    file.remove(destfiles[i])
  }
  return(data.frame(package_name = apps, genre = genres, stringsAsFactors = F))
}

apps <- read.csv(snakemake@input[[1]], stringsAsFactors = F)
genre_catalogue <- data.frame()
catalogue_source <- snakemake@params[["catalogue_source"]]
update_catalogue_file <- snakemake@params[["update_catalogue_file"]]
scrape_missing_genres <- snakemake@params[["scrape_missing_genres"]]
apps_with_genre <- data.frame(matrix(ncol=length(colnames(apps)) + 1,nrow=0, dimnames=list(NULL, c(colnames(apps), "genre"))))

if(nrow(apps) > 0){
  if(catalogue_source == "GOOGLE"){
    apps_with_genre <- apps %>% mutate(genre = NA_character_)
  } else if(catalogue_source == "FILE"){
    genre_catalogue <- read.csv(snakemake@params[["catalogue_file"]], colClasses = c("character", "character"))
    apps_with_genre <- left_join(apps, genre_catalogue, by = "package_name")
  }

  if(catalogue_source == "GOOGLE" || (catalogue_source == "FILE" && scrape_missing_genres)){
    apps_without_genre <- (apps_with_genre %>% filter(is.na(genre)) %>% distinct(package_name))$package_name
    updated_apps <- get_genre(apps_without_genre)
    apps_with_genre <- left_join(apps_with_genre, updated_apps, by = "package_name") %>%
      mutate(genre = coalesce(genre.x, genre.y)) %>%
      select(-genre.x, -genre.y)

    if(update_catalogue_file){
      genre_catalogue <- bind_rows(genre_catalogue, updated_apps) %>% distinct()
      write.csv(genre_catalogue, file = snakemake@params[["catalogue_file"]], row.names = FALSE)
    }
  }
}

write.csv(apps_with_genre, snakemake@output[[1]], row.names = FALSE)