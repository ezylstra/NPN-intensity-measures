# Exploring fruiting intensity data for oaks trees in midwest
# ER Zylstra

library(tidyverse)
library(leaflet)
library(elevatr)
library(sf)
library(terra)
library(ordinal)
library(lme4)

# Load data on white/red oak from Ohio DNR (for comparison) -------------------#
white_trees <- read.csv("oh-dnr-oak-data/percent-white-oak-trees-with-acorns.csv")
red_trees <- read.csv("oh-dnr-oak-data/percent-red-oak-trees-with-acorns.csv")
white_crowns <- read.csv("oh-dnr-oak-data/average-percent-of-white-crowns-with-acorns.csv")
red_crowns <- read.csv("oh-dnr-oak-data/average-percent-of-red-crowns-with-acorns.csv")

white_trees <- white_trees %>% mutate(common_name = "white oak")
red_trees <- red_trees %>% mutate(common_name = "northern red oak")
white_crowns <- white_crowns %>% mutate(common_name = "white oak")
red_crowns <- red_crowns %>% mutate(common_name = "northern red oak")

oh_trees <- rbind(white_trees, red_trees)
oh_crowns <- rbind(white_crowns, red_crowns)

# Load and format status-intensity data ---------------------------------------#

# Formatted file with just acorn data
ffile <- "npn-data/intensity-oaks/oak-acorns-2014-2025.csv"

if (file.exists(ffile)) {
  
  si <- read.csv(ffile)
  
} else {

  # List files with intensity data
  intensity_files <- list.files("npn-data/intensity-oaks",
                                full.names = TRUE)
  intensity_files <- intensity_files[grepl("intensity-oaks-", intensity_files)]
  
  for (i in 1:length(intensity_files)) {
    filename <- intensity_files[i]
    sitemp <- read.csv(filename)
    
    # Remove any true duplicates
    sitemp <- sitemp[!duplicated(sitemp),]
    
    if (i == 1) {
      si <- sitemp
    } else {
      si <- rbind(si, sitemp)
    }
    rm(sitemp)
  }  
  
  # Occasionally there are two records for a plant-phenophase in one day with 
  # different intensity or status values. Keeping record that was in phase with 
  # highest intensity value (or record that doesn't have NAs)
  si <- si %>%
    group_by(individual_id, phenophase_id, observation_date) %>%
    arrange(individual_id, observation_date, phenophase_id, 
            desc(phenophase_status), desc(intensity_midpoint)) %>%
    distinct(individual_id, phenophase_id, observation_date, .keep_all = TRUE) %>%
    data.frame()
  
  # Doublecheck that there's only one observation of each plant-phenophase per day
  if (nrow(si) != nrow(distinct(si, individual_id, 
                                phenophase_id, observation_date))) {
    warning("There is more than one observation of some plant phenophases at ",
            sitecap, " in a day")
  }
  
  # Remove any observations with phenophase_status = ?
  si <- si %>%
    filter(phenophase_status != -1)
  
  # Are there observations where the plant is out of phase but an intensity value
  # is reported? If so, change the status but add a column to note that the data
  # were amended.
  si <- si %>%
    mutate(amended_status = ifelse(phenophase_status == 0 & 
                                     !is.na(intensity_midpoint), 
                                   1, 0)) %>%
    mutate(phenophase_status = ifelse(phenophase_status == 0 & 
                                        !is.na(intensity_midpoint),
                                      1, phenophase_status))
  
  # Check that there's only one intensity category for each species-phenophase?
  spil <- si %>%
    filter(!is.na(intensity_category_id)) %>%
    distinct(common_name, phenophase_description, intensity_name,
             intensity_type, intensity_label) %>%
    filter(!is.na(intensity_name))
  count(spil, common_name, phenophase_description) %>%
    pull(n) %>% 
    max()
  # Value should be 1
  
  # Filter data by species, location --------------------------------------------#
  
  # Data for each species
  count(si, common_name)
  # Post oak has < 5000 records, while all other species have > 13000. 
  si <- si %>%
    filter(common_name != "post oak")
  
  sites <- si %>%
    group_by(site_id, site_name, state, 
             latitude, longitude, elevation_in_meters) %>%
    summarize(n_spp = n_distinct(common_name), 
              n_plants = n_distinct(individual_id),
              .groups = "keep") %>%
    rename(lat = latitude,
           lon = longitude,
           elev = elevation_in_meters) %>%
    data.frame()
  
  # Map
  leaflet(sites) %>% addTiles() %>%
    addCircleMarkers(lng = ~lon, lat = ~lat, radius = ~n_plants, 
                     fillOpacity = 0.6)
  
  # For now, limit to lat > 36.5 (about southern border of VA, KY)
  # and lon > (-95) (~western border of MO)
  sites <- sites %>% 
    filter(lat > 36.5 & lat < 49) %>%
    filter(lon > (-95) & lon < (-65))
  si <- si %>%
    filter(site_id %in% unique(sites$site_id))
  
  sites <- si %>%
    group_by(site_id, site_name, state, 
             latitude, longitude, elevation_in_meters) %>%
    summarize(n_plants = n_distinct(individual_id),
              black = ifelse("black oak" %in% common_name, 1, 0),
              bur = ifelse("bur oak" %in% common_name, 1, 0),
              red = ifelse("northern red oak" %in% common_name, 1, 0),
              pin = ifelse("pin oak" %in% common_name, 1, 0),
              swamp = ifelse("swamp white oak" %in% common_name, 1, 0),
              white = ifelse("white oak" %in% common_name, 1, 0),
              .groups = "keep") %>%
    rename(lat = latitude,
           lon = longitude,
           elev = elevation_in_meters) %>%
    data.frame()
  
  # Map by species
  leaflet(sites) %>% addTiles() %>%
    addCircles(lng = ~lon, lat = ~lat, 
               data = filter(sites, black == 1), 
               group = "black oak",
               color = "black", fillOpacity = 1, weight = 10) %>%
    addCircles(lng = ~lon, lat = ~lat, 
               data = filter(sites, bur == 1), 
               group = "bur oak",
               color = "purple", fillOpacity = 1, weight = 10) %>%
    addCircles(lng = ~lon, lat = ~lat, 
               data = filter(sites, red == 1), 
               group = "northern red oak",
               color = "red", fillOpacity = 1, weight = 10) %>%
    addCircles(lng = ~lon, lat = ~lat, 
               data = filter(sites, pin == 1), 
               group = "pin oak",
               color = "yellow", fillOpacity = 1, weight = 10) %>%
    addCircles(lng = ~lon, lat = ~lat, 
               data = filter(sites, swamp == 1), 
               group = "swamp white oak",
               color = "blue", fillOpacity = 1, weight = 10) %>%  
    addCircles(lng = ~lon, lat = ~lat, 
               data = filter(sites, white == 1), 
               group = "white oak",
               color = "green", fillOpacity = 1, weight = 10) %>%
    addLayersControl(overlayGroups = c("black oak",
                                       "bur oak",
                                       "northern red oak",
                                       "pin oak",
                                       "swamp white oak",
                                       "white oak"),
                     options = layersControlOptions(collapse = FALSE))
  
  # Some sites are missing elevation (and some may be wrong). Will grab elevations 
  # for all sites using the elevatr package
  elev_fill <- sites %>%
    distinct(site_id, lat, lon) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326)
  elev_fill <- get_elev_point(locations = elev_fill, src = "epqs")
  elev_fill <- data.frame(elev_fill) %>%
    mutate(elev_new = round(elevation)) %>%
    select(site_id, elev_new)
  elevations <- sites %>%
    left_join(elev_fill, by = "site_id")
  sites <- sites %>%
    left_join(select(elevations, site_id, elev_new), by = "site_id") %>% 
    select(-elev) %>%
    rename(elev = elev_new)
  # rm(elevations)
  
  # Some state assignments are missing (and a few look wrong).
  states <- terra::vect("states/cb_2017_us_state_500k.shp")
  state_fill <- sites %>%
    distinct(site_id, lat, lon)
  state_fillv <- vect(state_fill, 
                      geom = c("lon", "lat"), 
                      crs = "epsg:4326")
  state_new <- terra::extract(states, state_fillv)
  state_fill <- cbind(state_fill, state_new = state_new$STUSPS)
  sites <- sites %>%
    left_join(select(state_fill, site_id, state_new), by = "site_id") %>%
    data.frame() %>%
    mutate(state_best = ifelse(state %in% c("ON", "QC"), state, state_new)) %>%
    select(-c(state_new, state)) %>%
    rename(state = state_best)
  
  # Replace state, elevation in si dataframe with corrected values
  si <- si %>%
    select(-c(elevation_in_meters, state)) %>%
    left_join(select(sites, site_id, elev, state), by = "site_id")
  
  # Write status-intensity data to file
  write.csv(si, ffile, row.names = FALSE)
}

# Create or load weather data -------------------------------------------------#
# All data from PRISM (temps in degC and precip in mm)

# Extract site info
sites <- si %>%
  distinct(site_id, latitude, longitude) %>%
  rename(lon = longitude,
         lat = latitude)

# Specify location of PRISM data
prism_folder <- "C:/Users/erin/Documents/PRISMData/"

# Folder with weather data
weather_folder <- "weather-data/oaks/"

# Weather variables
weather_vars <- c("summer_temps_mn")

# Weather data
weather_files <- paste0(weather_folder, weather_vars, ".csv")

# Logical indicating whether to update weather data or not
update_weather <- FALSE

if (!all(file.exists(weather_files)) | update_weather) {
  
  # First, put site locations to a SpatVector, and convert to NAD83 to match
  # PRISM data
  sitesv <- vect(sites, geom = c("lon", "lat"), crs = "epsg:4326")
  # Convert to NAD83 to match PRISM data
  sitesv <- terra::project(sitesv, "epsg:4269")
  
  # Find years we have PRISM data for
  tmin_files <- list.files(prism_folder, pattern = "tmin-", full.names = TRUE)
  yrs <- str_sub(tmin_files, -8, -5) %>% as.numeric()

  for (yr in yrs) {
    
    # Mean summer temperatures ------------------------------------------------#
    # Get daily minimum temperatures in May-July
    rtmin <- readRDS(paste0(prism_folder, "tmin-", yr, ".rds"))
    rtmin_may <- rtmin[[str_sub(names(rtmin), 10, 11) == "05"]]
    tmin_may <- terra::extract(rtmin_may, sitesv, ID = FALSE)
    rtmin_jun <- rtmin[[str_sub(names(rtmin), 10, 11) == "06"]]
    tmin_jun <- terra::extract(rtmin_jun, sitesv, ID = FALSE)
    rtmin_jul <- rtmin[[str_sub(names(rtmin), 10, 11) == "07"]]
    tmin_jul <- terra::extract(rtmin_jul, sitesv, ID = FALSE)
    
    # Get daily maximum temperatures in May-July
    rtmax <- readRDS(paste0(prism_folder, "tmax-", yr, ".rds"))
    rtmax_may <- rtmax[[str_sub(names(rtmax), 10, 11) == "05"]]
    tmax_may <- terra::extract(rtmax_may, sitesv, ID = FALSE)
    rtmax_jun <- rtmax[[str_sub(names(rtmax), 10, 11) == "06"]]
    tmax_jun <- terra::extract(rtmax_jun, sitesv, ID = FALSE)
    rtmax_jul <- rtmax[[str_sub(names(rtmax), 10, 11) == "07"]]
    tmax_jul <- terra::extract(rtmax_jul, sitesv, ID = FALSE)
    
    # Calculate daily mean tempeartures
    tmn_daily_may <- (tmin_may + tmax_may)/2
    tmn_daily_jun <- (tmin_jun + tmax_jun)/2
    tmn_daily_jul <- (tmin_jul + tmax_jul)/2
    
    # Get monthly means
    tmn <- data.frame(site_id = sites$site_id,
                      season = yr,
                      may = apply(tmn_daily_may ,1, mean),
                      jun = apply(tmn_daily_jun ,1, mean),
                      jul = apply(tmn_daily_jul ,1, mean))
    rm(rtmin, rtmin_may, tmin_may, rtmin_jun, tmin_jun, rtmin_jul, tmin_jul,
       rtmax, rtmax_may, tmax_may, rtmax_jun, tmax_jun, rtmax_jul, tmax_jul,
       tmn_daily_may, tmn_daily_jun, tmn_daily_jul)

    # Combine data across years
    if (yr == min(yrs)) {
      tmn_all <- tmn
    } else {
      tmn_all <- rbind(tmn_all, tmn)
    }
  }
  
  # Save to file
  write.csv(tmn_all,
            "weather-data/oaks/summer_temps_mn.csv",
            row.names = FALSE)
}

# Load weather data
for (i in 1:length(weather_vars)) {
  assign(weather_vars[i], read.csv(weather_files[i]))
}

# Intensity categories for each phenophase ------------------------------------#
count(si, phenophase_description, intensity_name, 
      intensity_type, intensity_value, intensity_midpoint, intensity_label)
# For Fruits: No. fruits (largest bin = More than 10,000)
# For Ripe fruits: Ripe fruit %
# For Recent fruit or seed drop: No. fruit/seed drop (larges = More than 10,000)

# Just interested in number of fruit, so will remove Ripe Fruit observations
si <- si %>%
  filter(phenophase_description != "Ripe fruits")

# Clean up intensity data and remove NAs --------------------------------------#
# If status = 0, then set intensity value to 0

# Will eventually remove observations from the dataset that have status = 1 but
# no intensity value. However, I want to do a few things with the status data
# firstIf status = 1 and intensity not reported, remove observation from dataset

si <- si %>%
  mutate(intensity_midpoint = ifelse(phenophase_status == 0, 
                                     0, intensity_midpoint))

# check:
count(si, phenophase_description, phenophase_status, intensity_midpoint)

# Amount of data by species, year
si %>%
  group_by(common_name, yr) %>%
  summarize(n_fruit_plants = n_distinct(individual_id[phenophase_description == "Fruits"]),
            n_fruit_obs = sum(phenophase_description == "Fruits"),
            n_fruit_yes = sum(phenophase_status[phenophase_description == "Fruits"] == 1),
            n_drop_plants = n_distinct(individual_id[phenophase_description != "Fruits"]),
            n_drop_obs = sum(phenophase_description != "Fruits"),
            n_drop_yes = sum(phenophase_status[phenophase_description != "Fruits"] == 1),
            .groups = "keep") %>%
  mutate(prop_fruit_yes = round(n_fruit_yes/n_fruit_obs, 2),
         prop_drop_yes = round(n_drop_yes/n_drop_obs, 2)) %>%
  data.frame()
# Not too bad...

# When are acorns detected/dropped? -------------------------------------------#

# All species combined:
si %>% 
  select(phenophase_description, phenophase_status, day_of_year) %>%
  ggplot(aes(x = day_of_year)) +
  geom_histogram() +
  facet_grid(phenophase_status ~ phenophase_description)
# By species:
si %>%
  filter(phenophase_status == 1) %>%
  select(phenophase_description, common_name, day_of_year) %>%
  ggplot(aes(x = day_of_year)) +
  geom_histogram() +
  geom_vline(xintercept = 180, col = "steelblue3") +
  facet_grid(common_name ~ phenophase_description)
# Smaller seasonal window for fruit drop (almost all > doy 180 or 200)
# Swamp white oak has very little data, especially for fruit drop

# Prep data from OH DNR surveys (white, red oak only) for comparison ----------#
# Percent of trees with acorns (~status data)
trees_yr <- data.frame(
  yr = 2014:2024,
  oh_white = apply(white_trees[, 2:12], 2, mean, na.rm = TRUE),
  oh_red = apply(red_trees[, 2:12], 2, mean, na.rm = TRUE)) %>%
  mutate(across(oh_white:oh_red, ~round(., 2)))
cor(trees_yr[, 2:3])
trees_yr %>%
  pivot_longer(cols = oh_white:oh_red,
               names_to = "species",
               values_to = "percent") %>%
  ggplot(aes(x = yr, y = percent, group = species)) +
  geom_point(aes(color = species)) +
  geom_line(aes(color = species)) +
  scale_x_continuous(breaks = 2014:2024, labels = 2014:2024)
# Red and white oak in OH not doing the same thing each year

# Percent of crowns with acorns (~intensity data)
crowns_yr <- data.frame(
  yr = 2014:2024,
  oh_white = apply(white_crowns[, 2:12], 2, mean, na.rm = TRUE),
  oh_red = apply(red_crowns[, 2:12], 2, mean, na.rm = TRUE)) %>%
  mutate(across(oh_white:oh_red, ~round(., 2)))
cor(crowns_yr[, 2:3])
crowns_yr %>%
  pivot_longer(cols = oh_white:oh_red,
               names_to = "species",
               values_to = "percent") %>%
  ggplot(aes(x = yr, y = percent, group = species)) +
  geom_point(aes(color = species)) +
  geom_line(aes(color = species)) +
  scale_x_continuous(breaks = 2014:2024, labels = 2014:2024)
# Red and white oak in OH not doing the same thing each year (except for first
# few years [2014-2016])

# Look at proportion of NN trees with acorns each year ------------------------#
# This is just summarizing status data for fruits phenophase

fr <- si %>%
  filter(phenophase_description == "Fruits") %>%
  # Removing observations from early in year
  filter(day_of_year >= 100)

# Aggregate data for each plant-year
fr_py <- fr %>%
  group_by(common_name, individual_id, site_id, state, yr) %>%
  summarize(n_obs = n(),
            n_inphase = sum(phenophase_status),
            .groups = "keep") %>%
  data.frame() 

# Need to restrict summaries/analyses to plants that were observed a minimum
# number of times each year
sum(fr_py$n_obs >= 3)/nrow(fr_py) # 78% of plant-years had at least 3 obs
sum(fr_py$n_obs >= 5)/nrow(fr_py) # 66% of plant-years had at least 5 obs
min_obs <- 3
fr_py <- fr_py %>% filter(n_obs >= min_obs)

fr_yr <- fr_py %>%
  mutate(red = ifelse(common_name == "northern red oak", 1, 0),
         white = ifelse(common_name == "white oak", 1, 0)) %>%
  group_by(yr) %>%
  summarize(n_plants = n(),
            n_fruit = sum(n_inphase > 0),
            n_red = sum(red == 1),
            n_red_fruit = sum(n_inphase[red == 1] > 0),
            n_white = sum(white == 1),
            n_white_fruit = sum(n_inphase[white == 1] > 0)) %>%
  mutate(nn_all = round(n_fruit/n_plants * 100, 2),
         nn_red = round(n_red_fruit/n_red * 100, 2),
         nn_white = round(n_white_fruit /n_white * 100, 2)) %>%
  data.frame()

# Plot time series for each group of trees:
fr_yr %>%
  left_join(trees_yr, by = "yr") %>%
  data.frame() %>%
  filter(!is.na(oh_white)) %>%
  select(yr, nn_all, nn_red, nn_white, oh_white, oh_red) %>%
  pivot_longer(cols = nn_all:oh_red,
               names_to = "group",
               values_to = "percent") %>%
  ggplot(aes(x = yr, y = percent, group = group)) +
  geom_point(aes(color = group)) +
  geom_line(aes(color = group)) +
  scale_color_manual(values = c("darkred", "red", "coral", "darkblue", "lightblue3")) +
  scale_x_continuous(breaks = 2014:2024, labels = 2014:2024) +
  labs(y = "Percent of trees with acorns", x = "Year",
       title = "Comparing NN data (fruit) with OH DNR data") +
  theme_bw()

# Correlation coefficients:
fr_yr %>% 
  select(yr, contains("nn_")) %>%
  left_join(trees_yr, by = "yr") %>%
  filter(!is.na((oh_white))) %>%
  select(-yr) %>%
  cor() %>%
  round(2)
# White oaks at Ohio DNR site seem to be doing something different than all the 
# rest. Very weak positive correlation between red oak data (Ohio and NN).

# Look at proportion of NN trees that drop acorns each year -------------------#
# This is just summarizing status data for fruit drop phenophase

fd <- si %>%
  filter(phenophase_description == "Recent fruit or seed drop") %>%
  filter(common_name != "swamp white oak") %>%
  filter(day_of_year >= 180)

# Aggregate data for each plant-year
fd_py <- fd %>%
  group_by(common_name, individual_id, site_id, state, yr) %>%
  summarize(n_obs = n(),
            n_inphase = sum(phenophase_status),
            .groups = "keep") %>%
  data.frame() 

# Need to restrict summaries/analyses to plants that were observed a minimum
# number of times each year
sum(fd_py$n_obs >= 3)/nrow(fd_py) # 81% of plant-years had at least 3 obs
sum(fd_py$n_obs >= 5)/nrow(fd_py) # 73% of plant-years had at least 5 obs
min_obs <- 3
fd_py <- fd_py %>% filter(n_obs >= min_obs)

fd_yr <- fd_py %>%
  mutate(red = ifelse(common_name == "northern red oak", 1, 0),
         white = ifelse(common_name == "white oak", 1, 0)) %>%
  group_by(yr) %>%
  summarize(n_plants = n(),
            n_fruit = sum(n_inphase > 0),
            n_red = sum(red == 1),
            n_red_fruit = sum(n_inphase[red == 1] > 0),
            n_white = sum(white == 1),
            n_white_fruit = sum(n_inphase[white == 1] > 0)) %>%
  mutate(nn_all = round(n_fruit/n_plants * 100, 2),
         nn_red = round(n_red_fruit/n_red * 100, 2),
         nn_white = round(n_white_fruit /n_white * 100, 2)) %>%
  data.frame()

# Plot time series for each group of trees:
fd_yr %>%
  left_join(trees_yr, by = "yr") %>%
  data.frame() %>%
  filter(!is.na(oh_white)) %>%
  select(yr, nn_all, nn_red, nn_white, oh_white, oh_red) %>%
  pivot_longer(cols = nn_all:oh_red,
               names_to = "group",
               values_to = "percent") %>%
  ggplot(aes(x = yr, y = percent, group = group)) +
  geom_point(aes(color = group)) +
  geom_line(aes(color = group)) +
  scale_color_manual(values = c("darkred", "red", "coral", "darkblue", "lightblue3")) +
  scale_x_continuous(breaks = 2014:2024, labels = 2014:2024) +
  labs(y = "Percent of trees with acorns (OH) or dropped acorns (NN)", 
       x = "Year",
       title = "Comparing NN data (dropped fruit) with OH DNR data") +
  theme_bw()

# Correlation coefficients:
fd_yr %>% 
  select(yr, contains("nn_")) %>%
  left_join(trees_yr, by = "yr") %>%
  filter(!is.na((oh_white))) %>%
  select(-yr) %>%
  cor() %>%
  round(2)
# White oaks at Ohio DNR sites weakly correlated with NN white oaks. Patterns in 
# NN data (red oak or across all spp) somewhat correlated with Ohio DNR red 
# oaks, though mostly in 2017-2024.

# Annual summaries of intensity data ------------------------------------------#
# Remove observations where tree is in phase, but no intensity value is recorded
sii <- si %>%
  filter(!(phenophase_status == 1 & is.na(intensity_midpoint)))

# Fruit counts
  fri <- sii %>%
    filter(phenophase_description == "Fruits") %>%
    # Removing observations from early in year
    filter(day_of_year >= 100)
  
  # Aggregate data for each plant-year
  fri_py <- fri %>%
    group_by(common_name, individual_id, site_id, state, yr) %>%
    summarize(n_obs = n(),
              n_inphase = sum(phenophase_status),
              max_count = max(intensity_midpoint),
              .groups = "keep") %>%
    data.frame() 
  
  # Need to restrict summaries/analyses to plants that were observed a minimum
  # number of times each year
  sum(fri_py$n_obs >= 3)/nrow(fri_py) # 78% of plant-years had at least 3 obs
  sum(fri_py$n_obs >= 5)/nrow(fri_py) # 66% of plant-years had at least 5 obs
  min_obs <- 3
  fri_py <- fri_py %>% filter(n_obs >= min_obs)
  
  # Is there a correlation between max annual fruit counts and the number of
  # observations?
  ggplot(fri_py, aes(x = n_obs, y = log(max_count + 1))) +
    geom_point()
  cor(fri_py$n_obs, fri_py$max_count)
  # Weak (0.28)
  
  # How much variation among years for an individual?
  fri_p <- fri_py %>%
    group_by(common_name, individual_id, site_id, state) %>%
    summarize(n_yrs = n(),
              n_yrs_fruit = sum(n_inphase > 0),
              max_0 = sum(max_count == 0),
              max_10 = sum(max_count %in% c(1, 5)),
              max_50 = sum(max_count == 50),
              max_500 = sum(max_count == 500),
              max_5000 = sum(max_count == 5000),
              max_10001 = sum(max_count == 10001),
              .groups = "keep") %>%
    data.frame()
  fri_p %>%
    filter(n_yrs >= 3) %>%
    group_by(common_name) %>%
    summarize(n_plants = n(),
              # Number of plants that had acorns in at least one year
              n_fruit = sum(n_yrs_fruit > 0),
              # Number of plants that had acorns in at least one year but not every year
              n_vary = sum(n_yrs_fruit > 0 & n_yrs_fruit < n_yrs),
              # Number of plants that had annual max counts of 0 and annual max counts >= 500
              n_highlow = sum(max_0 > 0 & max_500 + max_5000 + max_10001 > 0),
              prop_highlow = round(n_highlow/n_fruit, 2),
              .groups = "keep") %>%
    data.frame()
  # All species except swamp white oaks have quite a bit of variation in max
  # acorn counts
  
# Dropped fruit counts
  fdi <- sii %>%
    filter(phenophase_description == "Recent fruit or seed drop") %>%
    filter(common_name != "swamp white oak") %>%
    filter(day_of_year >= 180)
  
  # Aggregate data for each plant-year
  fdi_py <- fdi %>%
    group_by(common_name, individual_id, site_id, latitude, longitude, state, yr) %>%
    summarize(n_obs = n(),
              n_inphase = sum(phenophase_status),
              max_count = max(intensity_midpoint),
              .groups = "keep") %>%
    data.frame() 
  
  # Need to restrict summaries/analyses to plants that were observed a minimum
  # number of times each year
  sum(fdi_py$n_obs >= 3)/nrow(fdi_py) # 80% of plant-years had at least 3 obs
  sum(fdi_py$n_obs >= 5)/nrow(fdi_py) # 72% of plant-years had at least 5 obs
  min_obs <- 3
  fdi_py <- fdi_py %>% filter(n_obs >= min_obs)
  
  # Is there a correlation between max annual fruit counts and the number of
  # observations?
  ggplot(fdi_py, aes(x = n_obs, y = log(max_count + 1))) +
    geom_point()
  cor(fdi_py$n_obs, fdi_py$max_count)
  # Weak (0.17)
  
  # How much variation among years for an individual?
  fdi_p <- fdi_py %>%
    group_by(common_name, individual_id, site_id, state) %>%
    summarize(n_yrs = n(),
              n_yrs_fruit = sum(n_inphase > 0),
              max_0 = sum(max_count == 0),
              max_10 = sum(max_count %in% c(1, 5)),
              max_50 = sum(max_count == 50),
              max_500 = sum(max_count == 500),
              max_5000 = sum(max_count == 5000),
              max_10001 = sum(max_count == 10001),
              .groups = "keep") %>%
    data.frame()
  fdi_p %>%
    filter(n_yrs >= 3) %>%
    group_by(common_name) %>%
    summarize(n_plants = n(),
              # Number of plants that had acorns in at least one year
              n_fruit = sum(n_yrs_fruit > 0),
              # Number of plants that had acorns in at least one year but not every year
              n_vary = sum(n_yrs_fruit > 0 & n_yrs_fruit < n_yrs),
              # Number of plants that had annual max counts of 0 and annual max counts >= 500
              n_highlow = sum(max_0 > 0 & max_500 + max_5000 + max_10001 > 0),
              prop_highlow = round(n_highlow/n_fruit, 2),
              .groups = "keep") %>%
    data.frame()
  # All species have quite a bit of variation in max acorn counts

  # Look more closely at spatial/temporal replication in dataset
  fdi_siteyr <- fdi_py %>%
    group_by(common_name, state, site_id, latitude, longitude, yr) %>%
    summarize(n_plants = n_distinct(individual_id), .groups = "keep") %>%
    data.frame()
  fdi_site <- fdi_siteyr %>%
    group_by(common_name, site_id, latitude, longitude) %>%
    summarize(n_yrs = n(),
              n_yrs_3min = sum(n_plants >= 3),
              .groups = "keep") %>%
    data.frame()
  
  # For each species, number of sites with x years where >= 3 trees monitored
  fdi_site %>%
    count(common_name, n_yrs_3min) %>%
    filter(n_yrs_3min > 0)
  
  # Map by species
  leaflet(fdi_site) %>% addTiles() %>%
    addCircles(lng = ~longitude, lat = ~latitude, 
               data = filter(fdi_site, common_name == "black oak", n_yrs_3min > 0), 
               group = "black oak",
               weight = ~n_yrs_3min*2, 
               color = "black", 
               fillOpacity = 1) %>%
    addCircles(lng = ~longitude, lat = ~latitude, 
               data = filter(fdi_site, common_name == "bur oak", n_yrs_3min > 0), 
               group = "bur oak",
               weight = ~n_yrs_3min*2, 
               color = "purple", 
               fillOpacity = 1)  %>%
    addCircles(lng = ~longitude, lat = ~latitude, 
               data = filter(fdi_site, common_name == "northern red oak", n_yrs_3min > 0), 
               group = "northern red oak",
               weight = ~n_yrs_3min*2, 
               color = "red", 
               fillOpacity = 1)  %>%
    addCircles(lng = ~longitude, lat = ~latitude, 
               data = filter(fdi_site, common_name == "pin oak", n_yrs_3min > 0), 
               group = "pin oak",
               weight = ~n_yrs_3min*2, 
               color = "blue", 
               fillOpacity = 1)  %>%
    addCircles(lng = ~longitude, lat = ~latitude, 
               data = filter(fdi_site, common_name == "white oak", n_yrs_3min > 0), 
               group = "white oak",
               weight = ~n_yrs_3min*2, 
               color = "green", 
               fillOpacity = 1)  %>%  
    addLayersControl(overlayGroups = c("black oak",
                                       "bur oak",
                                       "northern red oak",
                                       "pin oak",
                                       "white oak"),
                     options = layersControlOptions(collapse = FALSE))
  
  # Sites with a good number of years for black, northern red, and white oak
  # in the NYC area and north (to Albany)
  
  
# Not sure if stuff below is good or not.... 
  
  
# Regression models for max annual acorn count --------------------------------#
# (will eventually want to check if changing min number of obs/yr affects results)

# Remove swamp white oak data (not many observations with acorns)
fri_py <- fri_py %>%
  filter(common_name != "swamp white oak")
  
# Are there enough observations with acorns in 2025 (or is it too soon?)
fri_py %>%
  group_by(common_name, yr) %>%
  summarize(prop_acorns = sum(n_inphase > 0)/n()) %>%
  mutate(prop_acorns = round(prop_acorns, 2)) %>%
  data.frame()
# Proportion of plants where acorns detected doesn't seem lower in 2025, so 
# probably ok

# Create abundance categories
  # 0 (none), 1-50 (few), 500 (some), 5000-10001 (many)
  # May need to combine some and many
  acorns <- fri_py %>%
    mutate(abund4 = case_when(
      max_count == 0 ~ "none",
      max_count %in% 1:50 ~ "few",
      max_count == 500 ~ "some",
      max_count > 500 ~ "many",
    )) %>%
    mutate(abund4 = factor(abund4, 
                           levels = c("none", "few", "some", "many"),
                           ordered = TRUE)) %>%
    mutate(abund3 = case_when(
      max_count == 0 ~ "none",
      max_count %in% 1:50 ~ "few",
      max_count >= 500 ~ "many",
    )) %>%
    mutate(abund3 = factor(abund3, 
                           levels = c("none", "few", "many"),
                           ordered = TRUE))
  
  # Convert variables to factor
  acorns$fyr <- factor(acorns$yr)
  acorns$site <- factor(acorns$site_id)
  acorns$id <- factor(acorns$individual_id)

  # Where do we have data? (and append lat/lon to acorn data)
  sites_acorn <- sii %>%
    group_by(site_id, latitude, longitude) %>%
    summarize(n_plants = n_distinct(individual_id),
              black = ifelse("black oak" %in% common_name, 1, 0),
              bur = ifelse("bur oak" %in% common_name, 1, 0),
              red = ifelse("northern red oak" %in% common_name, 1, 0),
              pin = ifelse("pin oak" %in% common_name, 1, 0),
              white = ifelse("white oak" %in% common_name, 1, 0),
              .groups = "keep") %>%
    data.frame() %>%
    rename(lat = latitude,
           lon = longitude) %>%
    filter(site_id %in% acorns$site_id) %>%
    mutate(lat_z = (lat - mean(lat)) / sd(lat),
           lon_z = (lon - mean(lon)) / sd(lon))
  
  # Map by species
  leaflet(sites_acorn) %>% addTiles() %>%
    addCircles(lng = ~lon, lat = ~lat, 
               data = filter(sites_acorn, black == 1), 
               group = "black oak",
               color = "black", fillOpacity = 1, weight = 10) %>%
    addCircles(lng = ~lon, lat = ~lat, 
               data = filter(sites_acorn, bur == 1), 
               group = "bur oak",
               color = "purple", fillOpacity = 1, weight = 10) %>%
    addCircles(lng = ~lon, lat = ~lat, 
               data = filter(sites_acorn, red == 1), 
               group = "northern red oak",
               color = "red", fillOpacity = 1, weight = 10) %>%
    addCircles(lng = ~lon, lat = ~lat, 
               data = filter(sites_acorn, pin == 1), 
               group = "pin oak",
               color = "blue", fillOpacity = 1, weight = 10) %>%
    addCircles(lng = ~lon, lat = ~lat, 
               data = filter(sites_acorn, white == 1), 
               group = "white oak",
               color = "green", fillOpacity = 1, weight = 10) %>%
    # Add circle around Minneapolis/St. Paul
    addCircles(lng = -93.203542, lat = 44.962490, radius = 100000,
               weight = 2, color = "black") %>%
    # Add circle around NYC
    addCircles(lng = -74.00597, lat = 40.71427, radius = 100000,
               weight = 2, color = "black") %>%
    addLayersControl(overlayGroups = c("black oak",
                                       "bur oak",
                                       "northern red oak",
                                       "pin oak",
                                       "white oak"),
                     options = layersControlOptions(collapse = FALSE))
  # Note: we do have data for all species in the St. Paul MN area and the NYC 
  # areas...
  
  acorns <- acorns %>%
    left_join(select(sites_acorn, site_id, lat, lon, lat_z, lon_z), 
              by = "site_id")
  
  # Models with >2 abundance categories:
    # Wanted to run ML ordinal models (clmm or clmm2) with site and/or plant 
    # random effects but they were taking forever to run and eventually had 
    # convergence issues.
    # m_yr1 <- clmm(abund4 ~ fyr + (1|site), Hess = TRUE, nAGQ = 10,
    #               data = acorns)
    # m_yr1 <- clmm2(abund4 ~ fyr, random = site, Hess = TRUE, nAGQ = 10,
    #                data = acorns)
    
    # Tried a model with lat/lon and year as random effect:
    # m_ll <- clmm(abund4 ~ lat_z + I(lat_z^2) + lon_z + I(lon_z^2) + (1|fyr),
    #               Hess = TRUE, nAGQ = 10, data = acorns)
    # summary(m_ll)
    # Both quadratic effects signficant (max counts peak at intermediate lats,
    # and to a lesser degree, extreme longitudes)
    # Random year effects not that big/important

  # Maybe for this dataset, better to just create 2 categories (lots of acorns
  # versus few/none)
  acorns <- acorns %>%
    mutate(mast = ifelse(max_count >= 500, 1, 0))

  # Try a logistic regression model
  m_1 <- glmer(mast ~ fyr + lat_z + I(lat_z^2) + lon_z + I(lon_z^2) + (1|id), 
               data = acorns, 
               family = binomial, 
               control = glmerControl(optimizer = "bobyqa"),
               nAGQ = 10)
  summary(m_1)
  # Both lat/lon quadratic terms significant
  # Bunch of years had signficant positive effects: 2018-2019, 2023-2025
  # Plant REs seem important

  # Limit analyses to two areas (MSP and NYC)
  locs2 <- data.frame(loc = c("MSP", "NYC"),
                      lon = c(-93.203542, -74.00597),
                      lat = c(44.962490, 40.71247))
  locs2 <- vect(locs2, geom = c("lon", "lat"), crs = "epsg:4326")
  locs2_buffer <- terra::buffer(locs2, width = 100000)
  sites_acornv <- vect(sites_acorn, geom = c("lon", "lat"), crs = "epsg:4326")
  regions <- terra::extract(locs2_buffer, sites_acornv)
  sites_acorn$region <- regions$loc
  
  # check:
  # leaflet(sites_acorn) %>% addTiles() %>%
  #   addCircles(lng = ~lon, lat = ~lat, 
  #              data = filter(sites_acorn, region == "MSP"), 
  #              group = "MSP",
  #              color = "blue", fillOpacity = 1, weight = 10) %>%
  #   addCircles(lng = ~lon, lat = ~lat, 
  #              data = filter(sites_acorn, region == "NYC"), 
  #              group = "NYC",
  #              color = "purple", fillOpacity = 1, weight = 10) %>%
  #   addCircles(lng = ~lon, lat = ~lat, 
  #              data = filter(sites_acorn, is.na(region)), 
  #              group = "none",
  #              color = "red", fillOpacity = 1, weight = 10) %>%
  #   # Add circle around Minneapolis/St. Paul
  #   addCircles(lng = -93.203542, lat = 44.962490, radius = 100000,
  #              weight = 2, color = "black") %>%
  #   # Add circle around NYC
  #   addCircles(lng = -74.00597, lat = 40.71427, radius = 100000,
  #              weight = 2, color = "black") %>%
  #   addLayersControl(overlayGroups = c("MSP", "NYC", "none"),
  #                    options = layersControlOptions(collapse = FALSE))
    
  # Add region to acorn dataset
  acorns <- acorns %>%
    left_join(select(sites_acorn, site_id, region), by = "site_id")
  
  # Create new dataset for 2 regions:
  acorns2 <- acorns %>%
    filter(!is.na(region))
  # How much data do we have for each species?
  acorns2 %>%
    group_by(region, common_name) %>%
    summarize(nplants = n_distinct(individual_id), 
              .groups = "keep") %>%
    data.frame()
  # How much data do we have for each year?
  acorns2 %>%
    group_by(region, yr) %>%
    summarize(nplants = n_distinct(individual_id), 
              .groups = "keep") %>%
    data.frame()  
  
  acorns2$region <- factor(acorns2$region)
  table(acorns2$mast, acorns2$fyr, acorns2$region)
  
  # Try a logistic regression model for 2 regions (not accounting for spp):
  # Removing 2014 data because there were no trees with lots of acorns in MSP
  m_r1 <- glmer(mast ~ fyr*region + (1|id), 
                data = filter(acorns2, yr > 2014), 
                family = binomial,
                control = glmerControl(optimizer = "bobyqa"),
                nAGQ = 10)
  summary(m_r1)
  
  # Remove random effects:
  m_r2 <- glm(mast ~ fyr*region, 
              data = filter(acorns2, yr > 2014), 
              family = binomial)
  summary(m_r2)
  
  # Not sure that this model is very helpful -- don't we already know what years
  # were better than others? Would be more interesting to see whether the
  # probability a tree produces a lot of acorns is affected by weather...
  
  # Review Kelly et al. 2025 for ideas about weather?
  # Mean summer (June-July) temperatures 1 and 2 years prior
  
  # Load weather data
  summert <- summer_temps_mn %>%
    mutate(jjtemp_0 = (jun + jul)/2) %>%
    rename(yr = season)
  # Append same year temps to acorns data
  acorns2 <- acorns2 %>%
    left_join(select(summert, site_id, yr, jjtemp_0), by = c("site_id", "yr"))
  # Append temp data with 1- or 2-year lags
  summert_lag1 <- summert %>%
    mutate(yr = yr + 1) %>%
    rename(jjtemp_1 = jjtemp_0)
  summert_lag2 <- summert %>%
    mutate(yr = yr + 2) %>%
    rename(jjtemp_2 = jjtemp_0)
  acorns2 <- acorns2 %>%
    left_join(select(summert_lag1, site_id, yr, jjtemp_1),
              by = c("site_id", "yr")) %>%
    left_join(select(summert_lag2, site_id, yr, jjtemp_2),
              by = c("site_id", "yr")) %>%
    mutate(jjtemp_0_z = (jjtemp_0 - mean(jjtemp_0)) / sd(jjtemp_0),
           jjtemp_1_z = (jjtemp_1 - mean(jjtemp_1)) / sd(jjtemp_1),
           jjtemp_2_z = (jjtemp_2 - mean(jjtemp_2, na.rm = TRUE)) / sd(jjtemp_2, na.rm = TRUE))
  
  # Temps in same year, last year, 2 years ago
  m_rw012 <- glmer(mast ~ jjtemp_0_z + jjtemp_1_z + jjtemp_2_z + region + (1|id), 
                   data = filter(acorns2, yr > 2014), 
                   family = binomial,
                   control = glmerControl(optimizer = "bobyqa"),
                   nAGQ = 10)
  summary(m_rw012)
  # Positive effect of summer temperatures in same year
  # Negative (but smaller) effect of summer temperatures in previous year
  
  m_rw01 <- glmer(mast ~ jjtemp_0_z + jjtemp_1_z + region + (1|id), 
                   data = filter(acorns2, yr > 2014), 
                   family = binomial,
                   control = glmerControl(optimizer = "bobyqa"),
                   nAGQ = 10)
  summary(m_rw01)
  
  m_rw01i <- glmer(mast ~ jjtemp_0_z * jjtemp_1_z + region + (1|id), 
                   data = filter(acorns2, yr > 2014), 
                   family = binomial,
                   control = glmerControl(optimizer = "bobyqa"),
                   nAGQ = 10)
  summary(m_rw01i)
  # No support for an interaction between temperatures
  
  m_rw01ii <- glmer(mast ~ jjtemp_0_z*region + jjtemp_1_z*region + (1|id), 
                    data = filter(acorns2, yr > 2014), 
                    family = binomial,
                    control = glmerControl(optimizer = "bobyqa"),
                    nAGQ = 10)
  summary(m_rw01ii)
  # No support for interactions between region and temperatures
  
  # Species differences?
  acorns2$spp <- factor(acorns2$common_name)
  m_rws01 <- glmer(mast ~ jjtemp_0_z + jjtemp_1_z + spp + region + (1|id),
                   data = filter(acorns2, yr > 2014), 
                   family = binomial,
                   control = glmerControl(optimizer = "bobyqa"),
                   nAGQ = 10)
  summary(m_rws01)
  # No support for adding species to model (though pin oak had lower probability
  # of masting than other species)
  
  # What about about temporal autocorrelation? (if a tree produced a lot of
  # acorns in one year is it more/less likely to produce a lot of acorns the 
  # next year?)
  
  # Northern red oak in NYC:
  acorns2 %>%
    filter(region == "NYC") %>%
    filter(common_name == "northern red oak") %>%
    group_by(individual_id) %>%
    mutate(maxc = max(max_count),
           nyrs = n()) %>%
    filter(maxc > 0 & nyrs > 1) %>%
    ungroup() %>%
    ggplot(aes(x = yr, y = log(max_count + 0.1), group = individual_id)) +
    geom_jitter(aes(color = id), width = 0.1, height = 0.2) +
    # geom_line(aes(color = id)) +
    # facet_wrap(~id) +
    theme_bw() +
    theme(legend.position = "none")
