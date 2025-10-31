# Canopy fullness, deciduous trees in GRSM
# ER Zylstra
# 31 October 2025

library(rnpn)
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(ordbetareg)
library(tidybayes)
library(terra)
# library(tidyterra) # Need this if I create a map

# Not sure I need these...
# library(broom)        # Convert model objects to data frames
# library(broom.mixed)  # Convert brms model objects to data frames


# Load shapefile with US state boundaries -------------------------------------#

states <- vect("states/cb_2017_us_state_500k.shp")

# Download data, basic formatting (if not done already) -----------------------#

grsm_data_file <- "npn-data/grsm-canopy.csv" #### If not too big

if(!file.exists(grsm_data_file)) {
  
  # Get species IDs for 6 deciduous trees
  tree_ids <- npn_species() %>%
    filter(common_name %in% c("American basswood",
                              "American beech",
                              "northern red oak",
                              "red maple",
                              "striped maple",
                              "sugar maple")) %>%
    pull(species_id)
  
  # Get IDs for sites in GRSM
  grsm_sites <- npn_stations() %>%
    filter(grepl("GRSM", station_name)) %>% 
    data.frame() %>%
    filter(latitude < 36) %>%
    pull(station_id)

  # Download data since 2016
  dl <- npn_download_status_data(
    request_source = "erinz",
    years = 2016:2025,
    species_ids = tree_ids,
    station_ids = grsm_sites,
    climate_data = FALSE) %>%
    data.frame()
  
  # Extract data for leaves phenophase ### Make sure that we don't need to trim whitespaces first#######
  df <- dl %>%
    filter(phenophase_description == "Flowers or flower buds")
  
  
  
  
  # There are some observations where the state field is missing. Will get state
  # for each site using shapefiles from the census bureau.
  state_fill <- df %>%
    distinct(site_id, latitude, longitude, state)
  state_fillv <- vect(state_fill, 
                      geom = c("longitude", "latitude"), 
                      crs = "epsg:4326")
  state_new <- terra::extract(states, state_fillv)
  state_fill <- cbind(state_fill, state_new = state_new$STUSPS)
  # Check how these compare to original assignments:
  count(state_fill, state, state_new)
  
  # Attach new state assignments
  df <- df %>%
    left_join(select(state_fill, site_id, state_new), by = "site_id") %>%
    select(-state) %>%
    rename(state = state_new)
  
  # Remove unnecessary columns
  df <- df %>%
    select(-c(update_datetime, genus, species, kingdom, phenophase_id,
              intensity_category_id, abundance_value)) 
  
  write.csv(df, saguaro_data_file, row.names = FALSE)
  rm(df, dl, state_fill, state_fillv, state_new, saguaro_id)
}

