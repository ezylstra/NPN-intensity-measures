# Exploring acorn production in oaks
# ER Zylstra

library(tidyverse)
library(leaflet)
library(sf)
library(terra)
library(lme4)

# Load and format status-intensity data ---------------------------------------#

# Folder with oak data
oak_folder <- "npn-data/intensity-oaks/"

# Find species-specific files (don't have years/numbers in filename)
oakfiles <- list.files(oak_folder, full.names = TRUE)
oakfiles <- oakfiles[!grepl("[0-9]", oakfiles)]

for (i in 1:length(oakfiles)) {
  filename <- oakfiles[i]
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
sum(si$phenophase_status == 0 & !is.na(si$intensity_midpoint))
# No such records

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
  
# Remove unnecessary columns and rename, format -------------------------------#

si <- si %>%
  select(site_id, latitude, longitude, elevation_in_meters, state, 
         common_name, individual_id, phenophase_description, yr, 
         observation_date, day_of_year, phenophase_status, intensity_value, 
         intensity_type, intensity_midpoint, intensity_label) %>%
  rename(site = site_id, 
         lat = latitude,
         lon = longitude,
         elev = elevation_in_meters,
         spp = common_name,
         id = individual_id,
         php = phenophase_description, 
         obsdate = observation_date,
         doy = day_of_year,
         status = phenophase_status,
         intensity = intensity_value) %>%
  mutate(obsdate = ymd(obsdate))

# Fill in or correct state where needed ---------------------------------------#

# Some state assignments are missing (and a few look wrong).
states <- terra::vect("states/cb_2017_us_state_500k.shp")
state_fill <- si %>%
  distinct(site, lat, lon)
state_fillv <- vect(state_fill, 
                    geom = c("lon", "lat"), 
                    crs = "epsg:4326")
state_new <- terra::extract(states, state_fillv)
state_fill <- cbind(state_fill, state = state_new$STUSPS)

# Exclude any sites outside the lower 48 states
si <- si %>%
  select(-state) %>%
  left_join(select(state_fill, site, state), by = "site") %>%
  data.frame() %>%
  filter(!is.na(state)) %>%
  filter(state != "AK")

# Filter data by location -----------------------------------------------------#

# First, create indicator for CA or eastern species
si <- si %>% 
  mutate(region = ifelse(spp %in% c("valley oak", "blue oak"), "CA", "east"))

# Only keeping valley oak, blue oak trees in CA
si <- si %>%
  filter(!(spp %in% c("valley oak", "blue oak") & state != "CA"))
  
# Look at where the rest are located:
sites <- si %>%
  group_by(site, lat, lon) %>%
  summarize(n_plants = n_distinct(id),
            n_spp = n_distinct(spp),
            black = ifelse("black oak" %in% spp, 1, 0),
            bur = ifelse("bur oak" %in% spp, 1, 0),
            red = ifelse("northern red oak" %in% spp, 1, 0),
            white = ifelse("white oak" %in% spp, 1, 0),
            valley = ifelse("valley oak" %in% spp, 1, 0),
            blue = ifelse("blue oak" %in% spp, 1, 0),
            .groups = "keep") %>%
  data.frame()

# Map by species
# leaflet(sites) %>% addTiles() %>%
#   addCircles(lng = ~lon, lat = ~lat,
#              data = filter(sites, black == 1),
#              group = "black oak",
#              color = "black", fillOpacity = 1, weight = 10) %>%
#   addCircles(lng = ~lon, lat = ~lat,
#              data = filter(sites, bur == 1),
#              group = "bur oak",
#              color = "purple", fillOpacity = 1, weight = 10) %>%
#   addCircles(lng = ~lon, lat = ~lat,
#              data = filter(sites, red == 1),
#              group = "northern red oak",
#              color = "red", fillOpacity = 1, weight = 10) %>%
#   addCircles(lng = ~lon, lat = ~lat,
#              data = filter(sites, white == 1),
#              group = "white oak",
#              color = "green", fillOpacity = 1, weight = 10) %>%
#   addCircles(lng = ~lon, lat = ~lat,
#              data = filter(sites, valley == 1),
#              group = "valley oak",
#              color = "orange", fillOpacity = 1, weight = 10) %>%
#   addCircles(lng = ~lon, lat = ~lat,
#              data = filter(sites, blue == 1),
#              group = "blue oak",
#              color = "blue", fillOpacity = 1, weight = 10) %>%
#   addLayersControl(overlayGroups = c("black oak",
#                                      "bur oak",
#                                      "northern red oak",
#                                      "white oak",
#                                      "valley oak",
#                                      "blue oak"),
#                    options = layersControlOptions(collapse = FALSE))

# Remove any sites for 4 eastern species that are west of -100 deg
si <- si %>%
  filter(!(region == "east" & lon < (-100)))

# Remove sites dataframe since it now contains sites we're no longer using
rm(sites)

# Replace NA intensity midpoint values with 0s when status is 0 ---------------#

si <- si %>% 
  mutate(intensity_midpoint = ifelse(status == 0, 0, intensity_midpoint))

# Identify first yes date for each tree, phenophase, year ---------------------#
# For now, won't worry about whether the yes was preceeded by a no within X days

firstyes <- si %>%
  arrange(id, php, obsdate) %>%
  group_by(id, yr, php) %>%
  summarize(n_obs = n(),
            n_yes = sum(status),
            first_yes = ifelse(n_yes > 0, doy[1], NA),
            .groups = "keep") %>%
  data.frame()

# Explore white/bur oak data, by site -----------------------------------------#

wo <- si %>%
  filter(spp == "white oak") %>%
  filter(php %in% c("Fruits", "Recent fruit or seed drop"))
wo_sites <- wo %>%
  group_by(site, lat, lon, state) %>%
  summarize(n_trees = n_distinct(id),
            n_yrs = n_distinct(yr),
            n_treeyrs = n_distinct(paste0(id, yr)),
            .groups = "keep") %>%
  data.frame()

# Look for sites that monitored at least 5 trees in at least 5 years
wo_gsites <- wo_sites %>%
  filter(n_trees >= 5 & n_yrs >= 5 & n_treeyrs >= 25)
wo_gsites
leaflet(wo_gsites) %>% addTiles() %>%
  addCircleMarkers(lng = ~lon, lat = ~lat, radius = ~ n_treeyrs/5)
wo_grsm <- wo %>%
  filter(site %in% c(11896, 11999, 12002))
  
wo_grsm_fr <- wo_grsm %>%
  filter(php == "Fruits") %>%
  filter(!is.na(intensity_midpoint)) %>%
  arrange(id, obsdate) %>%
  group_by(id, yr, php) %>%
  summarize(n_obs = n(),
            first_obs = min(doy),
            last_obs = max(doy),
            n_obs_100300 = sum(doy >= 100 & doy <= 300),
            .groups = "keep") %>%
  # Identify trees observed at least 3 times between days 100 and 300
  mutate(keep_fr = ifelse(n_obs_100300 > 2, 1, 0)) %>%
  data.frame()

wo_grsm_dr <- wo_grsm %>%
  filter(php == "Recent fruit or seed drop") %>%
  filter(!is.na(intensity_midpoint)) %>%
  arrange(id, obsdate) %>%
  group_by(id, yr, php) %>%
  summarize(n_obs = n(),
            first_obs = min(doy),
            last_obs = max(doy),
            n_obs_180 = sum(doy >= 180),
            .groups = "keep") %>%
  # Identify trees observed at least 3 times on or after day 180
  mutate(keep_dr = ifelse(n_obs_180 > 2, 1, 0)) %>%
  data.frame()

# Attach include/exclude indicators to blue dataframe
wo_grsm <- wo_grsm %>%
  left_join(select(wo_grsm_fr, id, yr, php, keep_fr),
            by = c("id", "yr", "php")) %>%
  left_join(select(wo_grsm_dr, id, yr, php, keep_dr),
            by = c("id", "yr", "php"))

# What does annual acorn production look like (using dropped acorn counts)?
annual_drop <- wo_grsm %>% 
  filter(php == "Recent fruit or seed drop") %>%
  filter(!is.na(intensity_midpoint)) %>%
  filter(keep_dr == 1) %>%
  group_by(id, yr) %>%
  summarize(max_ct = max(intensity_midpoint)) %>%
  data.frame() %>%
  group_by(yr) %>%
  summarize(n_trees = n(),
            n_drop = sum(max_ct > 0),
            p_drop = round(n_drop/n_trees, 2),
            mn_ct = round(mean(max_ct[max_ct > 0]), 2),
            md_ct = median(max_ct[max_ct > 0])) %>%
  data.frame()
annual_drop

# What does annual acorn production look like (using acorn counts)?
annual_acorn <- wo_grsm %>% 
  filter(php == "Fruits") %>%
  filter(!is.na(intensity_midpoint)) %>%
  filter(keep_fr == 1) %>%
  group_by(id, yr) %>%
  summarize(max_ct = max(intensity_midpoint)) %>%
  data.frame() %>%
  group_by(yr) %>%
  summarize(n_trees = n(),
            n_acorn = sum(max_ct > 0),
            p_acorn = round(n_acorn/n_trees, 2),
            mn_ct = round(mean(max_ct[max_ct > 0]), 2),
            md_ct = median(max_ct[max_ct > 0])) %>%
  data.frame()
annual_acorn

# Bur oak instead?
bo <- si %>%
  filter(spp == "bur oak") %>%
  filter(php %in% c("Fruits", "Recent fruit or seed drop"))
bo_sites <- bo %>%
  group_by(site, lat, lon, state) %>%
  summarize(n_trees = n_distinct(id),
            n_yrs = n_distinct(yr),
            n_treeyrs = n_distinct(paste0(id, yr)),
            .groups = "keep") %>%
  data.frame()
leaflet(bo_sites) %>% addTiles() %>%
  addCircleMarkers(lng = ~lon, lat = ~lat, radius = ~ n_treeyrs/2)

# There just aren't that many trees monitored with sufficient acorn counts when 
# limiting things to one or a few sites. Try a region instead?

# Explore white/bur oak data, by region/city ----------------------------------#

wo <- si %>%
  filter(spp %in% c("white oak", "bur oak")) %>%
  filter(php %in% c("Fruits", "Recent fruit or seed drop"))

wo_sites <- wo %>%
  group_by(site, lat, lon, state) %>%
  summarize(n_trees = n_distinct(id),
            n_yrs = n_distinct(yr),
            n_treeyrs = n_distinct(paste0(id, yr)),
            .groups = "keep") %>%
  data.frame()

# Limit analyses to MSP, NYT, Boston, Raleigh/Durham area (circles with radius = 100 km)
locs <- data.frame(loc = c("MSP", "NYC", "BOS", "RDU"),
                    lon = c(-93.203542, -74.00597, -71.06548, -78.63853),
                    lat = c(44.962490, 40.71247, 42.36449, 35.78883))
locs <- vect(locs, geom = c("lon", "lat"), crs = "epsg:4326")
locs_buffer <- terra::buffer(locs, width = 100000)
wo_sitesv <- vect(wo_sites, geom = c("lon", "lat"), crs = "epsg:4326")
cities <- terra::extract(locs_buffer, wo_sitesv)
wo_sites$city <- cities$loc

# check:
leaflet(wo_sites) %>% addTiles() %>%
  addCircleMarkers(lng = ~lon, lat = ~lat, radius = ~n_treeyrs/2,
             data = filter(wo_sites, city == "MSP"),
             color = "blue", fillOpacity = 1, weight = 10) %>%
  addCircleMarkers(lng = ~lon, lat = ~lat, radius = ~n_treeyrs/2,
             data = filter(wo_sites, city == "NYC"),
             color = "blue", fillOpacity = 1, weight = 10) %>%
  addCircleMarkers(lng = ~lon, lat = ~lat, radius = ~n_treeyrs/2,
             data = filter(wo_sites, city == "BOS"),
             color = "blue", fillOpacity = 1, weight = 10) %>%
  addCircleMarkers(lng = ~lon, lat = ~lat, radius = ~n_treeyrs/2,
             data = filter(wo_sites, city == "RDU"),
             color = "blue", fillOpacity = 1, weight = 10) %>%
  addCircles(lng = ~lon, lat = ~lat,
             data = filter(wo_sites, is.na(city)),
             group = "none",
             color = "red", fillOpacity = 1, weight = 10) %>%
  # Add circle around Minneapolis/St. Paul
  addCircles(lng = -93.203542, lat = 44.962490, radius = 100000,
             weight = 2, color = "black") %>%
  # Add circle around NYC
  addCircles(lng = -74.00597, lat = 40.71427, radius = 100000,
             weight = 2, color = "black") %>%
  # Add circle around Boston
  addCircles(lng = -71.06548, lat = 42.36449, radius = 100000,
             weight = 2, color = "black") %>%
  # Add circle around Raleigh
  addCircles(lng = -78.63853, lat = 35.78883, radius = 100000,
             weight = 2, color = "black")

wo <- wo %>%
  left_join(select(wo_sites, site, city), by = "site")

table(wo$city, wo$spp)
# Best bet might be white oaks in Boston area

city1 <- "BOS"
wocity <- wo %>%
  filter(city == city1)

wocity_fr <- wocity %>%
  filter(php == "Fruits") %>%
  filter(!is.na(intensity_midpoint)) %>%
  arrange(id, obsdate) %>%
  group_by(spp, id, yr, php) %>%
  summarize(n_obs = n(),
            first_obs = min(doy),
            last_obs = max(doy),
            n_obs_100300 = sum(doy >= 100 & doy <= 300),
            .groups = "keep") %>%
  # Identify trees observed at least 3 times between days 100 and 300
  mutate(keep_fr = ifelse(n_obs_100300 > 2, 1, 0)) %>%
  data.frame()

wocity_dr <- wocity %>%
  filter(php == "Recent fruit or seed drop") %>%
  filter(!is.na(intensity_midpoint)) %>%
  arrange(id, obsdate) %>%
  group_by(id, yr, php) %>%
  summarize(n_obs = n(),
            first_obs = min(doy),
            last_obs = max(doy),
            n_obs_180 = sum(doy >= 180),
            .groups = "keep") %>%
  # Identify trees observed at least 3 times on or after day 180
  mutate(keep_dr = ifelse(n_obs_180 > 2, 1, 0)) %>%
  data.frame()

# Attach include/exclude indicators to blue dataframe
wocity <- wocity %>%
  left_join(select(wocity_fr, id, yr, php, keep_fr),
            by = c("id", "yr", "php")) %>%
  left_join(select(wocity_dr, id, yr, php, keep_dr),
            by = c("id", "yr", "php"))

# What does annual acorn production look like (using acorn counts)?
annual_acorn <- wocity %>% 
  filter(php == "Fruits") %>%
  filter(!is.na(intensity_midpoint)) %>%
  filter(keep_fr == 1) %>%
  group_by(spp, id, yr) %>%
  summarize(max_ct = max(intensity_midpoint)) %>%
  data.frame() %>%
  group_by(spp, yr) %>%
  summarize(n_trees = n(),
            n_acorn = sum(max_ct > 0),
            p_acorn = round(n_acorn/n_trees, 2),
            mn_ct = round(mean(max_ct[max_ct > 0]), 2),
            md_ct = median(max_ct[max_ct > 0])) %>%
  data.frame()

# What does annual acorn production look like (using dropped acorn counts)?
annual_drop <- wocity %>% 
  filter(php == "Recent fruit or seed drop") %>%
  filter(!is.na(intensity_midpoint)) %>%
  filter(keep_dr == 1) %>%
  group_by(spp, id, yr) %>%
  summarize(max_ct = max(intensity_midpoint)) %>%
  data.frame() %>%
  group_by(spp, yr) %>%
  summarize(n_trees = n(),
            n_drop = sum(max_ct > 0),
            p_drop = round(n_drop/n_trees, 2),
            mn_ct = round(mean(max_ct[max_ct > 0]), 2),
            md_ct = median(max_ct[max_ct > 0])) %>%
  data.frame()

annual_acorn
annual_drop

annual_acorn %>% filter(n_trees > 10)
annual_drop %>% filter(n_trees > 10)

wocity %>% 
  filter(php == "Recent fruit or seed drop") %>%
  filter(!is.na(intensity_midpoint)) %>%
  filter(keep_dr == 1) %>%
  group_by(spp, id, yr) %>%
  summarize(max_ct = max(intensity_midpoint)) %>%
  filter(yr %in% c(2014:2015, 2022:2025)) %>%
  filter(spp == "white oak") %>%
  ggplot(aes(x = log(max_ct + 1))) +
  geom_histogram() +
  facet_grid(yr~.)

# Explore white oak in MA -----------------------------------------------------#

whitema <- si %>%
  filter(spp == "white oak" & state == "MA")

# Want to look at annual acorn production via the max count of acorns or dropped 
# acorns for each tree and year. To reduce bias, probably want to remove trees
# that weren't observed during appropriate period 
whitema %>%
  filter(status == 1) %>%
  ggplot(aes(x = doy)) +
  geom_histogram() +
  facet_grid(php ~ .)

whitema_fr <- whitema %>%
  filter(php == "Fruits") %>%
  filter(!is.na(intensity_midpoint)) %>%
  arrange(id, obsdate) %>%
  group_by(id, yr, php) %>%
  summarize(n_obs = n(),
            first_obs = min(doy),
            last_obs = max(doy),
            n_obs_125300 = sum(doy >= 125 & doy <= 300),
            .groups = "keep") %>%
  # Identify trees observed at least 3 times between days 125 and 300
  mutate(keep_fr = ifelse(n_obs_125300 > 2, 1, 0)) %>%
  data.frame()

whitema_dr <- whitema %>%
  filter(php == "Recent fruit or seed drop") %>%
  filter(!is.na(intensity_midpoint)) %>%
  arrange(id, obsdate) %>%
  group_by(id, yr, php) %>%
  summarize(n_obs = n(),
            first_obs = min(doy),
            last_obs = max(doy),
            n_obs_180 = sum(doy >= 180),
            .groups = "keep") %>%
  # Identify trees observed at least 3 times on or after day 180
  mutate(keep_dr = ifelse(n_obs_180 > 2, 1, 0)) %>%
  data.frame()

# Attach include/exclude indicators to blue dataframe
whitema <- whitema %>%
  left_join(select(whitema_fr, id, yr, php, keep_fr),
            by = c("id", "yr", "php")) %>%
  left_join(select(whitema_dr, id, yr, php, keep_dr),
            by = c("id", "yr", "php"))

# What does annual acorn production look like (using dropped acorn counts)?
annual_drop <- whitema %>% 
  filter(php == "Recent fruit or seed drop") %>%
  filter(!is.na(intensity_midpoint)) %>%
  filter(keep_dr == 1) %>%
  group_by(id, yr) %>%
  summarize(max_ct = max(intensity_midpoint)) %>%
  data.frame() %>%
  group_by(yr) %>%
  summarize(n_trees = n(),
            n_drop = sum(max_ct > 0),
            p_drop = round(n_drop/n_trees, 2),
            mn_ct = round(mean(max_ct[max_ct > 0]), 2),
            md_ct = median(max_ct[max_ct > 0])) %>%
  data.frame()

annual_drop
cor(filter(annual_drop, n_drop > 0)) %>% round(2)
# Not great sample sizes, but %trees and counts not correlated

# What does annual acorn production look like (using acorn counts)?
annual_acorn <- whitema %>% 
  filter(php == "Fruits") %>%
  filter(!is.na(intensity_midpoint)) %>%
  filter(keep_fr == 1) %>%
  group_by(id, yr) %>%
  summarize(max_ct = max(intensity_midpoint)) %>%
  data.frame() %>%
  group_by(yr) %>%
  summarize(n_trees = n(),
            n_acorn = sum(max_ct > 0),
            p_acorn = round(n_acorn/n_trees, 2),
            mn_ct = round(mean(max_ct[max_ct > 0]), 2),
            md_ct = median(max_ct[max_ct > 0])) %>%
  data.frame()

annual_acorn
cor(annual_acorn) %>% round(2)


# Explore valley oak data -----------------------------------------------------#

valley <- si %>%
  filter(spp == "valley oak")

# Want to look at annual acorn production via the max count of acorns or dropped 
# acorns for each tree and year. To reduce bias, probably want to remove trees
# that weren't observed during appropriate period 
valley %>%
  filter(status == 1) %>%
  ggplot(aes(x = doy)) +
  geom_histogram() +
  facet_grid(php ~ .)

valley_fr <- valley %>%
  filter(php == "Fruits") %>%
  filter(!is.na(intensity_midpoint)) %>%
  arrange(id, obsdate) %>%
  group_by(id, yr, php) %>%
  summarize(n_obs = n(),
            first_obs = min(doy),
            last_obs = max(doy),
            n_obs_100300 = sum(doy >= 100 & doy <= 300),
            .groups = "keep") %>%
  # Identify trees observed at least 3 times between days 100 and 300
  mutate(keep_fr = ifelse(n_obs_100300 > 2, 1, 0)) %>%
  data.frame()

valley_dr <- valley %>%
  filter(php == "Recent fruit or seed drop") %>%
  filter(!is.na(intensity_midpoint)) %>%
  arrange(id, obsdate) %>%
  group_by(id, yr, php) %>%
  summarize(n_obs = n(),
            first_obs = min(doy),
            last_obs = max(doy),
            n_obs_180 = sum(doy >= 180),
            .groups = "keep") %>%
  # Identify trees observed at least 3 times on or after day 180
  mutate(keep_dr = ifelse(n_obs_180 > 2, 1, 0)) %>%
  data.frame()

# Attach include/exclude indicators to blue dataframe
valley <- valley %>%
  left_join(select(valley_fr, id, yr, php, keep_fr),
            by = c("id", "yr", "php")) %>%
  left_join(select(valley_dr, id, yr, php, keep_dr),
            by = c("id", "yr", "php"))

# What does annual acorn production look like (using dropped acorn counts)?
annual_drop <- valley %>% 
  filter(php == "Recent fruit or seed drop") %>%
  filter(!is.na(intensity_midpoint)) %>%
  filter(keep_dr == 1) %>%
  group_by(id, yr) %>%
  summarize(max_ct = max(intensity_midpoint)) %>%
  data.frame() %>%
  group_by(yr) %>%
  summarize(n_trees = n(),
            n_drop = sum(max_ct > 0),
            p_drop = round(n_drop/n_trees, 2),
            mn_ct = round(mean(max_ct[max_ct > 0]), 2),
            md_ct = median(max_ct[max_ct > 0])) %>%
  data.frame()

annual_drop
cor(annual_drop) %>% round(2)
# Bigger sample sizes for 2012-2019 (30+ trees each year) and strong correlation
# between counts and %trees with dropped acorns

# What does annual acorn production look like (using acorn counts)?
annual_acorn <- valley %>% 
  filter(php == "Fruits") %>%
  filter(!is.na(intensity_midpoint)) %>%
  filter(keep_fr == 1) %>%
  group_by(id, yr) %>%
  summarize(max_ct = max(intensity_midpoint)) %>%
  data.frame() %>%
  group_by(yr) %>%
  summarize(n_trees = n(),
            n_acorn = sum(max_ct > 0),
            p_acorn = round(n_acorn/n_trees, 2),
            mn_ct = round(mean(max_ct[max_ct > 0]), 2),
            md_ct = median(max_ct[max_ct > 0])) %>%
  data.frame()

annual_acorn
cor(annual_acorn) %>% round(2)
# Bigger sample sizes for 2012-2019 (30+ trees each year) and strong correlation
# between counts and %trees with acorns

# Explore blue oak data -------------------------------------------------------#

blue <- si %>%
  filter(spp == "blue oak")

# Want to look at annual acorn production via the max count of acorns or dropped 
# acorns for each tree and year. To reduce bias, probably want to remove trees
# that weren't observed during appropriate period 
blue %>%
  filter(status == 1) %>%
  ggplot(aes(x = doy)) +
  geom_histogram() +
  facet_grid(php ~ .)

blue_fr <- blue %>%
  filter(php == "Fruits") %>%
  filter(!is.na(intensity_midpoint)) %>%
  arrange(id, obsdate) %>%
  group_by(id, yr, php) %>%
  summarize(n_obs = n(),
            first_obs = min(doy),
            last_obs = max(doy),
            n_obs_100300 = sum(doy >= 100 & doy <= 300),
            .groups = "keep") %>%
  # Identify trees observed at least 3 times between days 100 and 300
  mutate(keep_fr = ifelse(n_obs_100300 > 2, 1, 0)) %>%
  data.frame()

blue_dr <- blue %>%
  filter(php == "Recent fruit or seed drop") %>%
  filter(!is.na(intensity_midpoint)) %>%
  arrange(id, obsdate) %>%
  group_by(id, yr, php) %>%
  summarize(n_obs = n(),
            first_obs = min(doy),
            last_obs = max(doy),
            n_obs_180 = sum(doy >= 180),
            .groups = "keep") %>%
  # Identify trees observed at least 3 times on or after day 180
  mutate(keep_dr = ifelse(n_obs_180 > 2, 1, 0)) %>%
  data.frame()

# Attach include/exclude indicators to blue dataframe
blue <- blue %>%
  left_join(select(blue_fr, id, yr, php, keep_fr),
            by = c("id", "yr", "php")) %>%
  left_join(select(blue_dr, id, yr, php, keep_dr),
            by = c("id", "yr", "php"))

# What does annual acorn production look like (using dropped acorn counts)?
annual_drop <- blue %>% 
  filter(php == "Recent fruit or seed drop") %>%
  filter(!is.na(intensity_midpoint)) %>%
  filter(keep_dr == 1) %>%
  group_by(id, yr) %>%
  summarize(max_ct = max(intensity_midpoint)) %>%
  data.frame() %>%
  group_by(yr) %>%
  summarize(n_trees = n(),
            n_drop = sum(max_ct > 0),
            p_drop = round(n_drop/n_trees, 2),
            mn_ct = round(mean(max_ct[max_ct > 0]), 2),
            md_ct = median(max_ct[max_ct > 0])) %>%
  data.frame()
annual_drop

cor(annual_drop) %>% round(2)
# Mean/median counts not correlated with % of trees that dropped acorns
# Mean counts but not % of trees that dropped acorns correlated with number of trees...
ggplot(annual_drop, aes(x = n_trees, y = mn_ct)) +
  geom_point()
# BUT this correlation driven entirely by 2022-2024 data, when way more trees
# were observed, all in the Bakersfield/Tehachapi area


# What does annual acorn production look like (using acorn counts)?
annual_acorn <- blue %>% 
  filter(php == "Fruits") %>%
  filter(!is.na(intensity_midpoint)) %>%
  filter(keep_fr == 1) %>%
  group_by(id, yr) %>%
  summarize(max_ct = max(intensity_midpoint)) %>%
  data.frame() %>%
  group_by(yr) %>%
  summarize(n_trees = n(),
            n_acorn = sum(max_ct > 0),
            p_acorn = round(n_acorn/n_trees, 2),
            mn_ct = round(mean(max_ct[max_ct > 0]), 2),
            md_ct = median(max_ct[max_ct > 0])) %>%
  data.frame()
annual_acorn

cor(annual_acorn) %>% round(2)
# Mean/median counts weakly correlated with % of trees that had acorns
# Mean counts but not % of trees that dropped acorns correlated with number of trees...
ggplot(annual_drop, aes(x = n_trees, y = mn_ct)) +
  geom_point()
# BUT this correlation driven entirely by 2022-2024 data, when way more trees
# were observed, all in the Bakersfield/Tehachapi area

# Variation among years for BLB:
blue %>%
  filter(status == 1 & php == "Breaking leaf buds") %>%
  ggplot(aes(x = doy)) +
  geom_histogram() + facet_wrap(~yr)
blue %>%
  filter(status == 1 & php == "Breaking leaf buds") %>%
  group_by(yr) %>%
  summarize(n = n(),
            mn = round(mean(doy), 2),
            md = median(doy),
            sd = round(sd(doy), 2),
            cv = round(sd / mn * 100, 2)) %>%
  data.frame()


# Summarize data by year
blue %>%
  group_by(yr) %>%
  summarize(n_sites = n_distinct(site),
            n_trees = n_distinct(id),
            n_drop = sum(php == "Recent fruit or seed drop"),
            n_drop1 = sum(status == 1 & php == "Recent fruit or seed drop"), 
            n_blb = sum(php == "Breaking leaf buds"),
            n_blb1 = sum(status == 1 & php == "Breaking leaf buds"),
            n_flower = sum(php == "Flowers or flower buds"),
            n_flower1 = sum(status == 1 & php == "Flowers or flower buds"),
            n_open = sum(php == "Open flowers"),
            n_open1 = sum(status == 1 & php == "Open flowers"),
            n_pollen = sum(php == "Pollen release (flowers)"),
            n_pollen1 = sum(status == 1 & php == "Pollen release (flowers)")) %>%
  data.frame() %>%
  mutate(p_drop = round(n_drop1/n_drop, 2),
         p_blb = round(n_blb1/n_blb, 2),
         p_flower = round(n_flower1/n_flower, 2),
         p_open = round(n_open1/n_open, 2), 
         p_pollen = round(n_pollen1/n_pollen, 2))
