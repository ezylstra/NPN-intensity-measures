# Canopy fullness, deciduous trees in GRSM
# ER Zylstra
# 3 November 2025

library(rnpn)
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(zoo)
library(ggplot2)
library(ordbetareg)
library(tidybayes)
library(terra)
# library(tidyterra) # Need this if I create a map

# Not sure I need these...
# library(broom)        # Convert model objects to data frames
# library(broom.mixed)  # Convert brms model objects to data frames

# Load shapefile with US state boundaries -------------------------------------#

# states <- vect("states/cb_2017_us_state_500k.shp")

# Download data, basic formatting (if not done already) -----------------------#

grsm_data_file <- "npn-data/grsm-canopy.csv"

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
  
  # Extract data for leaves phenophase
  df <- dl %>%
    filter(phenophase_description == "Leaves")

  # Remove unnecessary columns
  df <- df %>%
    select(-c(update_datetime, genus, species, kingdom, phenophase_id,
              intensity_category_id, abundance_value)) 
  
  write.csv(df, grsm_data_file, row.names = FALSE)
  rm(df, dl, states, grsm_sites, tree_ids)
}

# Load GRSM leaves data and simplify ------------------------------------------#

df <- read.csv(grsm_data_file)

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
# sum(df$status == -1)/nrow(df) * 100 # 0.12% of observations (n = 60)
df <- filter(df, status != -1)

# Create numeric column with ~midpoints for each intensity category (for easier 
# sorting and coding)
df <- df %>%
  mutate(intensity_midpoint = case_when(
    intensity_value == "Less than 5%" ~ 2,
    intensity_value == "5-24%" ~ 14,
    intensity_value == "25-49%" ~ 37,
    intensity_value == "50-74%" ~ 62,
    intensity_value == "75-94%" ~ 84,
    intensity_value == "95% or more" ~ 95,
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
# count(df, status, intensity_midpoint)

# Make intensity midpoints = 0 if status = 0
df <- df %>%
  mutate(intensity_midpoint = ifelse(status == 0, 0, intensity_midpoint))
# Now the only NAs left in the intensity_midpoint column occur when the status 
# is yes but no intensity value was provided 

# Filter data -----------------------------------------------------------------#

# Will focus on canopy closure in spring, so limiting observations to 
# the first half of the year (day 182 is June 30/July 1)
df <- filter(df, doy <= 182)

# Calculate intervals between consecutive observations for each tree, year
df <- df %>%
  arrange(common_name, id, yr, doy)

df$interval <- NA
for (i in 2:nrow(df)) {
  if (df$id[i] == df$id[i - 1] & df$yr[i] == df$yr[i - 1]) {
    df$interval[i] <- df$doy[i] - df$doy[i - 1]
  } else {
    df$interval[i] <- NA
  }
}

# Want to:
# remove plant-year combos with no observations in phase
# remove plant-year combos with no observations with intensity values
# remove plant-year combos with < 5 observations
# remove plant-year combos when max interval > 21 days

# Summarize amount and quality of information for each plant, year
pl_yr <- df %>%
  group_by(site, common_name, id, yr) %>%
  summarize(nobs = n(),
            first_obs = min(doy),
            last_obs = max(doy),
            mean_int = round(mean(interval, na.rm = TRUE), 2),
            max_int = ifelse(nobs == 1, NA, max(interval, na.rm = TRUE)),
            # Get number of observations in phase
            n_inphase = sum(status),
            # Get number of observations in phase with an intensity value
            n_intvalue = sum(!is.na(intensity_value)),
            # Calculate proportion of observations in phase
            prop_inphase = round(n_inphase / nobs, 2),
            # Calculate proportion of in-phase observations with intensity value
            prop_intvalue = round(n_intvalue / n_inphase, 2),
            .groups = "keep") %>%
  data.frame()

pl_yr <- pl_yr %>%
  mutate(remove = case_when(
    n_inphase == 0 ~ 1,
    n_intvalue == 0 ~ 1,
    nobs < 5 ~ 1,
    max_int > 21 ~ 1,
    .default = 0
  ))

df <- df %>%
  left_join(select(pl_yr, id, yr, remove), by = c("id", "yr")) %>%
  filter(remove == 0) %>%
  select(-remove)  

# Looks like there might be instances where an individual tree had high canopy
# values with a 0 value in-between. This isn't really possible, at least in such
# a short amount of time, so we'll try to filter these data out since they're
# likely due to observation error.

# Look for patterns with consecutive intensity values >80, then <50, then
# back up
dff <- df %>%
  mutate(intensity_cat = ifelse(intensity_midpoint > 80, 2,
                                ifelse(intensity_midpoint < 50, 0, 1))) %>%
  # Remove observations in phase with no intensity value
  filter(!is.na(intensity_midpoint)) %>%
  arrange(id, obsdate)

# Where do high-low-high patterns occur?
pattern <- c(2, 0, 2)
starts <- which(zoo::rollapply(dff$intensity_cat, 
                               length(pattern), 
                               function(x) all(x == pattern)))

# Identify anomalously low observations so we can filter them out
dff$anomalous <- 0
dff$anomalous[starts + 1] <- 1
# sum(dff$anomalous == 1) / nrow(dff) * 100 # 0.64% of observations

# Filter them out
dff <- dff %>%
  filter(anomalous == 0) %>%
  select(-c(anomalous, intensity_cat))

# Summarize remaining data by species, year
spp_yr <- dff %>% 
  group_by(common_name, yr) %>%
  summarize(n_trees = n_distinct(id),
            .groups = "keep") %>%
  data.frame()
spp_yr
sum(spp_yr$n_trees >= 5) / nrow(spp_yr)
# 42/58 species-years (72%) have >= 5 trees. 

# Summarize remaining data by species, across years
spp_yr %>%
  group_by(common_name) %>%
  summarize(n_yrs = n_distinct(yr),
            n_treeyrs = sum(n_trees),
            trees_per_yr_min = min(n_trees),
            trees_per_yr_max = max(n_trees),
            trees_per_yr_mn = round(mean(n_trees), 2)) %>%
  data.frame()
# American basswood averages only 3.5 trees per year; total of 35 tree-years.

# Summarize data by year
dff %>%
  group_by(yr) %>%
  summarize(n_spp = n_distinct(common_name),
            n_trees = n_distinct(id),
            n_sites = n_distinct(site)) %>%
  data.frame()
# All years had data for all 6 species except for 2020 (4 species)
# Number of trees much lower in 2020 than all other years (20 vs 60 in 2025
# and 98-167 in all other years).
# Only 5 sites in 2020, but 12-19 in all other years

# Create a dataset that excludes 2020 since we'll likely want to include year
# as a factor in some models
dff_no20 <- filter(dff, yr != 2020)

# Extract and summarize information about sites -------------------------------#

# Summarize data by site
grsm_sites <- dff %>%
  group_by(site, lat, lon, elev, state) %>%
  summarize(n_yrs = n_distinct(yr),
            yr_2020 = ifelse(2020 %in% yr, 1, 0),
            n_spp = n_distinct(common_name),
            n_trees = n_distinct(id),
            n_bass = n_distinct(id[common_name == "American basswood"]),
            n_beech = n_distinct(id[common_name == "American beech"]),
            n_oak = n_distinct(id[common_name == "northern red oak"]),
            n_red = n_distinct(id[common_name == "red maple"]),
            n_striped = n_distinct(id[common_name == "striped maple"]),
            n_sugar = n_distinct(id[common_name == "sugar maple"]),
            .groups = "keep") %>%
  data.frame()
  
# Write sites information to file to get PRISM weather data
# Downloading PRISM data through web interface, so we could get data at 800m 
# resolution without having to download daily rasters for CONUS, which takes a
# long time (https://prism.oregonstate.edu/explorer/bulk.php)
# write.table(select(grsm_sites, lat, lon, site), "weather-data/grsm-sites.csv",
#             sep = ",", row.names = FALSE, col.names = FALSE)

# Process and attach weather data ---------------------------------------------#

# Daily data
prism_files <- list.files("weather-data/grsm",
                          full.names = TRUE,
                          pattern = "stable|provisional")
for (i in 1:length(prism_files)) {
  prism1 <- read.csv(prism_files[i],
                     header = FALSE,
                     skip = 11,
                     col.names = c("site", "lon", "lat", "elev",
                                   "date", "ppt", "tmin", "tmean", "tmax"))
  if (i == 1) {
    weather <- prism1
  } else {
    weather <- rbind(weather, prism1)
  }
}
rm(prism1)

# Calculate daily GDD values, with base temp of 0 deg C (gdd = tmean - base)
base <- 0
weather <- weather %>%
  mutate(gdd = ifelse(tmean < base, 0, tmean - base)) %>%
  mutate(yr = year(date),
         doy = yday(date)) %>%
  filter(doy <= 200) %>%
  arrange(site, date)

# Calculate AGDD values
weather <- weather %>%
  group_by(site, yr) %>%
  mutate(agdd = cumsum(gdd)) %>%
  ungroup() %>%
  data.frame() %>%
  rename(obsdate = date)

# Add weather data to leaves datasets (and standardize potential variables)
dff <- dff %>%
  left_join(select(weather, site, obsdate, agdd), by = c("site", "obsdate")) %>%
  mutate(agdd_z = (agdd - mean(agdd)) / sd(agdd),
         doy_z = (doy - mean(doy)) / sd(doy),
         elev_z = (elev - mean(elev)) / sd(elev))
dff_no20 <- dff_no20 %>%
  left_join(select(weather, site, obsdate, agdd), by = c("site", "obsdate")) %>%
  mutate(agdd_z = (agdd - mean(agdd)) / sd(agdd),
         doy_z = (doy - mean(doy)) / sd(doy),
         elev_z = (elev - mean(elev)) / sd(elev))

# Prepare data for models -----------------------------------------------------#

# Calculate canopy proportion (considering 95% as full), and
# convert some variables to factors
dff <- dff %>%
  mutate(prop = intensity_midpoint/95,
         fyr = factor(yr),
         spp = factor(common_name),
         id = factor(id))
dff_no20 <- dff_no20 %>%
  mutate(prop = intensity_midpoint/95,
         fyr = factor(yr),
         spp = factor(common_name),
         id = factor(id))

# Model canopy development as a function of date ------------------------------#

# Year as a random effect (excluding 2020 data)
# Crossed random spp, year effects
start_time1 <- Sys.time()
m_no20_doy <- ordbetareg(
  prop ~ doy_z + (1 + doy_z|fyr) + (1 + doy_z|spp) + (1|id),
  data = dff_no20,
  control = list(adapt_delta = 0.99),
  iter = 4000, cores = 4, chains = 4,
  backend = "cmdstanr")
end_time1 <- Sys.time()
end_time1 - start_time1
save(m_no20_doy, file = "output/grsm-models/no20-doy.RData")
# load("output/grsm-models/no20-doy.RData")

# Year as a random effect (excluding 2020 data)
# Crossed random spp, year effects with an interaction
start_time2 <- Sys.time()
m_no20_doy_REint <- ordbetareg(
  prop ~ doy_z + (1 + doy_z|fyr) + (1 + doy_z|spp) + (1 + doy_z|fyr:spp) + (1|id),
  data = dff_no20,
  control = list(adapt_delta = 0.99),
  iter = 4000, cores = 4, chains = 4,
  backend = "cmdstanr")
end_time2 <- Sys.time()
end_time2 - start_time2
save(m_no20_doy_REint, file = "output/grsm-models/no20-doy-REint.RData")
# load("output/grsm-models/no20-doy-REint.RData")

# Year as a random effect (including 2020 data)
# Crossed random spp, year effects
start_time3 <- Sys.time()
m_allyr_doy <- ordbetareg(
  prop ~ doy_z + (1 + doy_z|fyr) + (1 + doy_z|spp) + (1|id),
  data = dff,
  control = list(adapt_delta = 0.99),
  iter = 4000, cores = 4, chains = 4,
  backend = "cmdstanr")
end_time3 <- Sys.time()
end_time3 - start_time3
save(m_allyr_doy, file = "output/grsm-models/allyr-doy.RData")
# load("output/grsm-models/allyr-doy.RData")

# Model canopy development as a function of AGDD ------------------------------#

# Species as a random effect (excluding 2020 data)
start_time4 <- Sys.time()
m_no20_gdd <- ordbetareg(
  prop ~ agdd_z + (1 + agdd_z|spp) + (1|id),
  data = dff_no20, 
  control = list(adapt_delta = 0.99),
  iter = 4000, cores = 4, chains = 4,
  backend = "cmdstanr")
end_time4 <- Sys.time()
end_time4 - start_time4
save(m_no20_gdd, file = "output/grsm-models/no20-gdd.RData")
# load("output/grsm-models/no20-gdd.RData")

# Species as a random effect (including 2020 data)
start_time5 <- Sys.time()
m_allyr_gdd <- ordbetareg(
  prop ~ agdd_z + (1 + agdd_z|spp) + (1|id),
  data = dff, 
  control = list(adapt_delta = 0.99),
  iter = 4000, cores = 4, chains = 4,
  backend = "cmdstanr")
end_time5 <- Sys.time()
end_time5 - start_time5
save(m_allyr_gdd, file = "output/grsm-models/allyr-gdd.RData")
# load("output/grsm-models/allyr-gdd.RData")

# Compare models, results -----------------------------------------------------#

load("output/grsm-models/no20-doy.RData")
load("output/grsm-models/no20-doy-REint.RData")
load("output/grsm-models/allyr-doy.RData")
load("output/grsm-models/no20-gdd.RData")
load("output/grsm-models/allyr-gdd.RData")

# Compare models
summary(m_no20_doy)
summary(m_no20_doy_REint)
summary(m_allyr_doy)
summary(m_no20_gdd)
summary(m_allyr_gdd)

# Compare models that excluded 2020 data
loo_doy <- loo(m_no20_doy, cores = 4)
loo_doy
loo_doy_int <- loo(m_no20_doy_REint, cores = 4)
loo_doy_int
loo_gdd <- loo(m_no20_gdd, cores = 4)
loo_gdd
loo_compare(loo_doy, loo_doy_int, loo_gdd)

# Compare models that included 2020 data
loo_allyr_doy <- loo(m_allyr_doy, cores = 4)
loo_allyr_doy
loo_allyr_gdd <- loo(m_allyr_gdd, cores = 4)
loo_allyr_gdd
loo_allyr_compare(loo_allyr_doy, loo_allyr_gdd)

# Visualize results for best models -------------------------------------------#
