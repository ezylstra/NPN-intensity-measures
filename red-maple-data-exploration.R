# Exploring masting events in red maple
# ER Zylstra

library(tidyverse)
library(leaflet)
library(sf)
library(terra)
library(lme4)

# Load and format status-intensity data ---------------------------------------#

# Folder with red maple data
rema_folder <- "npn-data/intensity-red-maple/"

# Find yr-specific files (don't have years/numbers in filename)
remafiles <- list.files(rema_folder, 
                        pattern = "intensity-rema-20",
                        full.names = TRUE)

for (i in 1:length(remafiles)) {
  filename <- remafiles[i]
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

# Limit to numeric intensity data categories ----------------------------------#

# First see if there are any pollen data
si %>%
  filter(php == "Pollen release (flowers)") %>%
  count(status, intensity)
# There are some, so keep for now

php_n <- si %>%
  filter(!is.na(intensity_type) & intensity_type != "percent") %>%
  distinct(php) %>%
  pull()

si_full <- si
si <- si_full %>%
  filter(php %in% php_n)

# Replace NA intensity midpoint values with 0s when status is 0 ---------------#

si <- si %>% 
  mutate(intensity_midpoint = ifelse(status == 0, 0, intensity_midpoint))

# Filter data by location -----------------------------------------------------#

# Map
sites <- si %>%
  group_by(site, state, lat, lon) %>%
  summarize(n_trees = n_distinct(id),
            n_yrs = n_distinct(yr),
            n_treeyrs = n_distinct(paste0(id, yr)),
            .groups = "keep") %>%
  data.frame()

leaflet(sites) %>% addTiles() %>%
    addCircleMarkers(lng = ~lon, lat = ~lat, radius = ~n_treeyrs/10,
                     color = "blue")
sites %>% arrange(desc(n_treeyrs)) %>% head(10)
sites %>% filter(n_trees >= 10) %>% arrange(desc(n_yrs)) %>% head(20)

# Look at intensity data for a few individual sites:
si %>%
  filter(site == 28745) %>%
  count(site, state, php, status, intensity_midpoint)

# site 18944 in MA: flowers, mostly with intensity data
# site 20120 in MA: flowers, fruits, drop, mostly with intensity values

# site 9342 in TN: flowers, few fruits, mostly with intensity values
# site 11897 in NC: flowers, fruits, mostly with intensity values

# site 25151 in LA: flowers, fruits, drop, mostly with intensity values
# site 28745 in MS: open flowers, fruits, few drop, mostly with intensity

# Explore data from paired site ----------------------------------------------# 
# Is there evidence of masting?

site_ids <- c(25151, 28745) 

sid <- si %>%
  filter(site %in% site_ids) %>%
  filter(!is.na(intensity_midpoint)) %>%
  filter(php != "Pollen release (flowers)")

# How much data (and yes reports) for each phenophase?
table(sid$php, sid$status)
# How do php observations break down by year?
table(sid$php, sid$yr)

# When do phenophases occur here?
sid %>%
  filter(status == 1) %>%
  ggplot(aes(x = doy)) +
  geom_histogram() +
  facet_grid(php ~ .)
# Might need to use water year for MS/LA site, but ignore for now

sid_plant <- sid %>%
  arrange(php, id, obsdate) %>%
  group_by(php, id, yr) %>%
  summarize(n_obs = n(),
            first_obs = min(doy),
            last_obs = max(doy),
            n_obs_180 = sum(doy <= 180),
            .groups = "keep") %>%
  # Identify trees observed at least 5 times before day 180
  mutate(keep = ifelse(n_obs_180 > 4, 1, 0)) %>%
  data.frame()

# Attach include/exclude indicators to dataframe
sid <- sid %>%
  left_join(select(sid_plant, id, yr, php, keep),
            by = c("id", "yr", "php"))

# What do intensity values look like on an annual basis?
sid_yr <- sid %>% 
  filter(keep == 1) %>%
  group_by(php, id, yr) %>%
  summarize(max_ct = max(intensity_midpoint), .groups = "keep") %>%
  data.frame() %>%
  group_by(php, yr) %>%
  summarize(n_trees = n(),
            n_inphase = sum(max_ct > 0),
            p_inphase = round(n_inphase/n_trees, 2),
            mn_ct = round(mean(max_ct[max_ct > 0]), 2),
            md_ct = median(max_ct[max_ct > 0]),
            .groups = "keep") %>%
  data.frame()
sid_yr

# The number of trees observed at least 5 times in first half of each year is
# small and the nubmer of trees with fruit/drop counts is really small.
# Think it might be better to look at patterns across multiple sites within
# a region.

# What if we identified a bunch of "regions" that contained multiple trees
# observed over multiple years. Then looked at the degree of synchrony in
# flower production as a function of distance or latitude (maybe more 
# seasonal environments are more likely to mast?)

# Map just sites with flower data:
sitesf <- si %>%
  filter(php == "Flowers or flower buds") %>%
  group_by(site, state, lat, lon) %>%
  summarize(n_trees = n_distinct(id),
            n_yrs = n_distinct(yr),
            n_treeyrs = n_distinct(paste0(id, yr)),
            .groups = "keep") %>%
  data.frame()
leaflet(sitesf) %>% addTiles() %>%
  addCircleMarkers(lng = ~lon, lat = ~lat, radius = ~n_treeyrs/10,
                   color = "blue")

# Potential "regions":
  # Gulf coast: LA, MS, AL
  # Tallahassee or FL generally?
  # GRSM area (TN and western NC)
  # Charlotte area
  # Raleigh area
  # Blacksburg/Roanoke
  # Washington DC
  # NYC
  # Boston
  # Portland, ME
  # Buffalo
  # Chicago
  # MSP
  # Colorado? (Denver + Greeley)

leaflet(sitesf) %>% addTiles() %>%
  addCircleMarkers(lng = ~lon, lat = ~lat, radius = ~n_treeyrs/10,
                   color = "blue") %>%
  # Southern LA, MS (includes NO, Pascagoula)
  addCircles(lng = -89.33069, lat = 30.31549, radius = 100000,
             weight = 2, color = "red") %>%
  # GRSM/Asheville NC
  addCircles(lng = -83.24039, lat = 35.63242, radius = 75000,
             weight = 2, color = "red") %>%
  # Tallahassee ??
  addCircles(lng = -84.27848, lat = 30.44395, radius = 50000,
             weight = 2, color = "red") %>%
  # Charlotte
  addCircles(lng = -81.06826, lat = 35.33918, radius = 50000,
             weight = 2, color = "red") %>%
  # Raleigh/Durham
  addCircles(lng = -78.77166, lat = 35.91345, radius = 50000,
             weight = 2, color = "red") %>%
  # Blacksburg/Roanoke
  addCircles(lng = -80.18785, lat = 37.24514, radius = 50000,
             weight = 2, color = "red") %>%
  # Washington DC/Baltimore
  addCircles(lng = -76.93418, lat = 39.03193, radius = 50000,
             weight = 2, color = "red") %>%
  # NYC
  addCircles(lng = -74.23687, lat = 40.62056, radius = 50000,
             weight = 2, color = "red") %>%
  # Buffalo
  addCircles(lng = -78.78174, lat = 42.94499, radius = 50000,
             weight = 2, color = "red") %>%
  # Boston
  addCircles(lng = -71.11074, lat = 42.36948, radius = 50000,
             weight = 2, color = "red") %>%
  # Portland, ME
  addCircles(lng = -70.27476, lat = 43.68543, radius = 50000,
             weight = 2, color = "red") %>%
  # Chicago
  addCircles(lng = -87.94979, lat = 41.87950, radius = 50000,
             weight = 2, color = "red") %>%
  # MSP
  addCircles(lng = -93.16071, lat = 45.00024, radius = 50000,
             weight = 2, color = "red") %>%
  # Colorado
  addCircles(lng = -104.86602, lat = 40.10633, radius = 50000,
            weight = 2, color = "red") 

# Assign "city" to sites
nola <- data.frame(city = "nola", lon = -89.33069, lat = 30.31549, rad = 100000)
grsm <- data.frame(city = "grsm", lon = -83.24039, lat = 35.63242, rad = 75000)
tall <- data.frame(city = "tall", lon = -84.27848, lat = 30.44395, rad = 50000)
char <- data.frame(city = "char", lon = -81.06826, lat = 35.33918, rad = 50000)
radu <- data.frame(city = "radu", lon = -78.77166, lat = 35.91345, rad = 50000)
blro <- data.frame(city = "blro", lon = -80.18785, lat = 37.24514, rad = 50000)
wash <- data.frame(city = "wash", lon = -76.93418, lat = 39.03193, rad = 50000)
nyc <- data.frame(city = "nyc", lon = -74.23687, lat = 40.62056, rad = 50000)
buff <- data.frame(city = "buff", lon = -78.78174, lat = 42.94499, rad = 50000)
bos <- data.frame(city = "bos", lon = -71.11074, lat = 42.36948, rad = 50000)
port <- data.frame(city = "port", lon = -70.27476, lat = 43.68543, rad = 50000)
chi <- data.frame(city = "chi", lon = -87.94979, lat = 41.87950, rad = 50000)
msp <- data.frame(city = "msp", lon = -93.16071, lat = 45.00024, rad = 50000)
colo <- data.frame(city = "colo", lon = -104.86602, lat = 40.10633, rad = 50000)

locs <- rbind(nola, grsm, tall, char, radu, blro, wash, nyc, buff, bos, port,
              chi, msp, colo)
locs <- vect(locs, geom = c("lon", "lat"), crs = "epsg:4326")
locs_buffer <- terra::buffer(locs, width = locs$rad)
sitesfv <- vect(sitesf, geom = c("lon", "lat"), crs = "epsg:4326")
cities <- terra::extract(locs_buffer, sitesfv)
sitesf$city <- cities$city

sitesf %>%
  filter(!is.na(city)) %>%
  group_by(city) %>%
  summarize(tree_yrs = sum(n_treeyrs)) %>%
  data.frame()

# Link city assignment to si data
si_city <- si %>%
  left_join(select(sitesf, site, city), by = "site") %>%
  filter(!is.na(city))

si_fl <- si_city %>%
  filter(php == "Flowers or flower buds") %>%
  filter(!is.na(intensity_midpoint)) %>%
  arrange(city, id, obsdate) %>%
  group_by(city, id, yr, php) %>%
  summarize(n_obs = n(),
            first_obs = min(doy),
            last_obs = max(doy),
            n_obs_180 = sum(doy <= 180),
            .groups = "keep") %>%
  # Identify trees observed at least 3 times in first half of year
  mutate(keep = ifelse(n_obs_180 > 2, 1, 0)) %>%
  data.frame()

# Attach include/exclude indicators to blue dataframe
si_city <- si_city %>%
  left_join(select(si_fl, id, yr, php, keep),
            by = c("id", "yr", "php"))

# What does annual flower production look like?
annual_fl <- si_city %>% 
  filter(php == "Flowers or flower buds") %>%
  filter(!is.na(intensity_midpoint)) %>%
  filter(!is.na(keep) & keep == 1) %>%
  group_by(city, id, yr) %>%
  summarize(max_ct = max(intensity_midpoint)) %>%
  data.frame() %>%
  group_by(city, yr) %>%
  summarize(n_trees = n(),
            n_fl = sum(max_ct > 0),
            p_fl = round(n_fl/n_trees, 2),
            mn_ct = round(mean(max_ct[max_ct > 0]), 2),
            md_ct = median(max_ct[max_ct > 0])) %>%
  data.frame()

# Exclude data for a city and year when < 4 trees had confirmed flowers
annual_fl <- annual_fl %>%
  filter(n_fl >= 4)

table(annual_fl$city, annual_fl$yr)

annual_flw <- annual_fl %>%
  pivot_wider(id_cols = yr,
              names_from = city,
              values_from = mn_ct) %>%
  data.frame() %>%
  arrange(yr)
annual_flw

annual_flw %>% 
  filter(yr %in% 2013:2025) %>% 
  select(grsm, port) %>% 
  cor() %>% round(2) 
# 0.75

annual_flw %>% 
  filter(yr %in% 2017:2025) %>% 
  select(grsm, msp, port, nola, nyc) %>% 
  cor() %>% round(2) 
# Highest is 0.50 (port, grsm); most negative

annual_flw %>% 
  filter(yr %in% 2019:2025) %>% 
  select(grsm, msp, bos, port, nola, nyc) %>% 
  cor() %>% round(2) 
# Highest is 0.45 (nyc, nola); many negative

# What does annual fruit production look like? (not filtering first)
annual_fr <- si_city %>% 
  filter(php == "Fruits") %>%
  filter(!is.na(intensity_midpoint)) %>%
  # filter(!is.na(keep) & keep == 1) %>%
  group_by(city, id, yr) %>%
  summarize(max_ct = max(intensity_midpoint)) %>%
  data.frame() %>%
  group_by(city, yr) %>%
  summarize(n_trees = n(),
            n_fr = sum(max_ct > 0),
            p_fr = round(n_fr/n_trees, 2),
            mn_ct = round(mean(max_ct[max_ct > 0]), 2),
            md_ct = median(max_ct[max_ct > 0])) %>%
  data.frame()

annual_fr <- annual_fr %>%
  filter(n_fr >= 4)

table(annual_fr$city, annual_fr$yr)

annual_frw <- annual_fr %>%
  pivot_wider(id_cols = yr,
              names_from = city,
              values_from = mn_ct) %>%
  data.frame() %>%
  arrange(yr)
annual_frw

annual_frw %>% 
  filter(yr %in% c(2013:2015, 2017:2021, 2023:2025)) %>% 
  select(grsm, port) %>% 
  cor() %>% round(2) 
# 0.57

annual_frw %>% 
  filter(yr %in% c(2017:2021, 2023:2025)) %>% 
  select(grsm, port, nola) %>% 
  cor() %>% round(2) 
# 0.70 for port, grsm
# -0.67 for nola, grsm
# -0.43 for nola, port

# This is going to be tough given that most cities don't have data in many years
# or they had few trees monitored in many of the years.

# Cities that might be worth looking into because they had more data:
# Flowers: Portland, ME (2012-2025); New Orleans (2017-2025); GRSM (2013-2025)
# Fruit: maybe same as above?
