# Download weather data for McDowell sites
# ER Zylstra
# 15 July 2025

# Using PRISM, since Daymet data for 2024 are not yet available

library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(here)
library(httr2)
library(terra)

# List files with formatted intensity data from McDowell
intensity_files <- list.files("npn-data",
                              pattern = "intensity-siteMCDO",
                              full.names = TRUE)

# Load NPN data ---------------------------------------------------------------#

sites <- str_sub(intensity_files, 24, 28)

for (site in sites) {
  filename <- intensity_files[grepl(site, intensity_files)]
  si1 <- read.csv(filename)
  
  # Remove any true duplicates
  si1 <- si1[!duplicated(si1),]
  
  if (site == sites[1]) {
    si <- si1
  } else {
    si <- rbind(si, si1)
  }
  rm(si1)
}  

# Grab lat/lon for each McDowell site
locs <- si %>%
  distinct(site_id, site_name, latitude, longitude, elevation_in_meters) %>%
  rename(id = site_id, 
         lat = latitude,
         lon = longitude,
         elev = elevation_in_meters)
rm(si)

# Using 2017-2024 data, so want to download weather data from Oct 2016 through
# Dec 2024
dates_list <- list()
dates_list[[1]] <- seq.Date("2016-10-01", "2016-12-31")
for(yy in 2017:2024) {
  dates_list[[yy-2015]] <- seq.Date(paste0(yy, "-01-01"), 
                                    paste0(yy, "-12-31"))
}

# Create folder to hold temporary prism data
prism_folder <- paste0(here::here(), "/prismtmp/")
dir.create(prism_folder)

# Base URL for PRISM Web Service (see https://www.prism.oregonstate.edu/downloads/)
url_base <- "https://services.nacse.org/prism/data/get/us/800m"


# Download daily data and extract values for McDowell sites
# for (yy in 1:length(dates_list)) {
for (yy in 2:2) {

  dates <- dates_list[[yy]]
  
  for (i in 1:length(dates)) {
    date_nodash <- str_remove_all(dates[i], "-")
    
    # Create folder/zip names
    ppt_folder <- paste0("prism_ppt_us_30s_", date_nodash)
    ppt_zip <- paste0(prism_folder, ppt_folder, ".zip")
    tmp_folder <- paste0("prism_tmean_us_30s_", date_nodash)
    tmp_zip <- paste0(prism_folder, tmp_folder, ".zip")
    
    # Download precip data (mm) and mean temperature data (degC)
    req <- request(url_base)
    req_ppt <- req %>%
      req_url_path_append("ppt", date_nodash)
    req_perform(req_ppt, path = ppt_zip)
    req_tmp <- req %>%
      req_url_path_append("tmean", date_nodash)
    req_perform(req_tmp, path = tmp_zip)
    
    # Unzip folders
    suppressWarnings(
      utils::unzip(ppt_zip, exdir = paste0(prism_folder, ppt_folder))
    )
    suppressWarnings(
      utils::unzip(tmp_zip, exdir = paste0(prism_folder, tmp_folder))
    )
    # Remove zips
    invisible(file.remove(ppt_zip))
    invisible(file.remove(tmp_zip))
  
    # Read files in as SpatRasters
    ppt_file <- paste0(prism_folder, ppt_folder, "/", ppt_folder, ".tif")
    ppt <- terra::rast(ppt_file)
    tmp_file <- paste0(prism_folder, tmp_folder, "/", tmp_folder, ".tif")
    tmp <- terra::rast(tmp_file)
    
    # Extract ppt and temperature values
    vals_ppt <- extract(ppt, as.matrix(locs[, c("lon", "lat")]), cells = TRUE)
    vals_tmp <- extract(tmp, as.matrix(locs[, c("lon", "lat")]), cells = TRUE)
    weather_df1 <- data.frame(cell = vals_ppt[, 1],
                              site_id = locs[, "id"],
                              date = date_nodash,
                              ppt = vals_ppt[, 2],
                              tmean = vals_tmp[, 2])
    
    if (yy == 1 & i == 1) {
      weather_df <- weather_df1
    } else {
      weather_df <- rbind(weather_df, weather_df1)
    }
  }
}

weather_2017 <- weather_df
weather_2017 <- weather_2017 %>%
  left_join(select(locs, id, site_name), by = c("site_id" = "id"))
write.csv(weather_2016, "weather-data/mcdo-2017.csv", row.names = FALSE)

# unlink(prism_folder, recursive = TRUE)

