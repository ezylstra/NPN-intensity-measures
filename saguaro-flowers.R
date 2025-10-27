# Floral abundance, saguaros
# ER Zylstra
# 27 October 2025

library(rnpn)
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(ordinal)
# library(maps)
library(terra)
# library(tidyterra)

# Load shapefile with US state boundaries -------------------------------------#

states <- vect("states/cb_2017_us_state_500k.shp")

# Download data, basic formatting (if not done already) -----------------------#

saguaro_data_file <- "npn-data/saguaro-flowers.csv"

if(!file.exists(saguaro_data_file)) {
  
  # Get species ID for saguaro
  saguaro_id <- npn_species() %>%
    filter(common_name == "saguaro") %>%
    pull(species_id)
  
  # Download most recent 10 years of data
  dl <- npn_download_status_data(
    request_source = "erinz",
    years = 2016:2025,
    species_ids = saguaro_id,
    climate_data = FALSE) %>%
    data.frame()
  
  # Extract data for flowers/flower buds
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
  
  # Attach new state assignments to the saguaro data and restrict sites to Arizona
  df <- df %>%
    left_join(select(state_fill, site_id, state_new), by = "site_id") %>%
    select(-state) %>%
    rename(state = state_new) %>%
    filter(!is.na(state), state == "AZ")
  
  # Remove unnecessary columns
  df <- df %>%
    select(-c(update_datetime, genus, species, kingdom, phenophase_id,
              intensity_category_id, abundance_value)) 
  
  write.csv(df, saguaro_data_file, row.names = FALSE)
  rm(df, dl, state_fill, state_fillv, state_new, saguaro_id)
}

# Load saguaro flowering data and simplify ------------------------------------#

df <- read.csv(saguaro_data_file)

# Remove unnecessary columns, rename others, and create year column
df <- df %>%
  select(-c(observation_id, species_id)) %>%
  rename(site = site_id, 
         lat = latitude,
         lon = longitude, 
         elev = elevation_in_meters, 
         id = individual_id,
         php = phenophase_description, 
         obsdate = observation_date,
         doy = day_of_year,
         status = phenophase_status) %>%
  mutate(yr = year(obsdate))

# Remove observations where status was unknown (-1)
# sum(df$status == -1)/nrow(df) * 100 # 0.17% of observations (n = 44)
df <- filter(df, status != -1)

# Create numeric column with ~midpoints for each intensity category (for easier 
# sorting and coding)
df <- df %>%
  mutate(intensity_midpoint = case_when(
    intensity_value == "Less than 3" ~ 1,
    intensity_value == "3 to 10" ~ 5,
    intensity_value == "11 to 100" ~ 50,
    intensity_value == "101 to 1,000" ~ 500,
    intensity_value == "More than 1,000" ~ 1001,
    .default = NA
  ))
  
# Occasionally there are two records for a plant in one day with 
# different intensity or status values. Keeping record that was in phase with 
# highest intensity value (or record that doesn't have NAs)
df <- df %>%
  group_by(id, obsdate) %>%
  mutate(n_obs = n()) %>%
  ungroup() %>%
  data.frame() %>%
  arrange(id, obsdate, desc(status), desc(intensity_midpoint)) %>%
  distinct(id, obsdate, .keep_all = TRUE) %>%
  dplyr::select(-n_obs)

# Check that there aren't any observations where the plant was out of phase 
# but an intensity value was reported.
# count(status, intensity_midpoint)

# Make intensity midpoints = 0 if status = 0
df <- df %>%
  mutate(intensity_midpoint = ifelse(status == 0, 0, intensity_midpoint))
# Now the only NAs left in the intensity_midpoint column occur when the status 
# is yes but no intensity value was provided 

# Filter data -----------------------------------------------------------------#

# Will use the same data to identify when flowering occurs. 
# Specifically, will identify the days-of-year that bracket the period when 80% 
# of the flowering observations occurred (flowering period), and the 
# days-of-year that bracket the period when 50% of the flowering occurred (peak
# flowering period). 

# Then:
# Exclude plant-year combos when <5 observations were made during the flowering
# period and <2 observations were made during the peak flowering period. 

# Identify when saguaros are typically flowering
flowering_periods <- df %>%
  filter(status == 1) %>%
  group_by(common_name) %>%
  summarize(nyrs = n_distinct(yr),
            nplants = n_distinct(id),
            nplantyrs = n_distinct(paste(yr, id)),
            mean = round(mean(doy)),
            q0.10 = quantile(doy, probs = 0.10),
            q0.25 = quantile(doy, probs = 0.25),
            q0.50 = quantile(doy, probs = 0.50),
            q0.75 = quantile(doy, probs = 0.75),
            q0.90 = quantile(doy, probs = 0.90),
            .groups = "keep") %>%
  mutate(across(mean:q0.90, round)) %>%
  mutate(length_50 = q0.75 - q0.25,
         length_80 = q0.90 - q0.10) %>%
  data.frame()
flowering_periods

q0.10 <- flowering_periods$q0.10[1]
q0.25 <- flowering_periods$q0.25[1]
q0.75 <- flowering_periods$q0.75[1]
q0.90 <- flowering_periods$q0.90[1]

# Summarize amount and quality of information for each plant and year
pl_yr <- df %>%
  mutate(in50 = ifelse(doy >= q0.25 & doy <= q0.75, 1, 0)) %>%
  mutate(in80 = ifelse(doy >= q0.10 & doy <= q0.90, 1, 0)) %>%
  group_by(id, site, yr) %>%
  summarize(nobs = n(),
            first_obs = min(doy),
            last_obs = max(doy),
            nobs50 = sum(in50),
            nobs80 = sum(in80),
            n_inphase = sum(status),
            n_intvalue = sum(!is.na(intensity_value)),
            prop_inphase = round(n_inphase / nobs, 2),
            prop_intvalue = round(n_intvalue / n_inphase, 2),
            max_intensity = ifelse(sum(is.na(intensity_midpoint)) == n(),
                                   NA, max(intensity_midpoint, na.rm = TRUE)),
            .groups = "keep") %>%
  data.frame()

# Identify which plant-yr combinations to filter out
pl_yr <- pl_yr %>%
  mutate(remove = case_when(
    is.na(max_intensity) ~ 1,
    is.na(nobs50) ~ 1,
    is.na(nobs80) ~ 1,
    nobs80 < 5 ~ 1,
    nobs50 < 2 ~ 1,
    .default = 0
  ))

df_filter <- df %>%
  left_join(select(pl_yr, id, yr, remove), by = c("id", "yr")) %>%
  filter(remove == 0) %>%
  select(-remove)







