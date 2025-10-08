# Exploring fruiting intensity data for oaks trees in midwest
# ER Zylstra

library(tidyverse)
library(leaflet)
library(elevatr)
library(sf)
library(spData)
library(terra)
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

# List files with formatted intensity data
intensity_files <- list.files("npn-data/intensity-oaks",
                              full.names = TRUE)

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
rm(elevations)

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
# rest. Patterns in NN data (red oak particularly) look a little similar to 
# Ohio DNR red oaks, but almost always a lower percentage.

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

# Intensity data --------------------------------------------------------------#
# Remove observations where tree is in phase, but no intensity value is recorded
sii <- si %>%
  filter(!(phenophase_status == 1 & is.na(intensity_midpoint)))
