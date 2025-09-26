# Exploring intensity data for eastern redbud
# ER Zylstra

library(tidyverse)
library(leaflet)
library(elevatr)
library(sf)
library(terra)

# Specify location of PRISM data ----------------------------------------------#

prism_folder <- "C:/Users/erin/Documents/PRISMData/"

# Load and format status-intensity data ---------------------------------------#

# List files with formatted intensity data
intensity_files <- list.files("npn-data/intensity-redbud",
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

# Filter data: plant-phenophase-year combinations -----------------------------#

sites <- si %>%
  group_by(site_id, site_name, state, 
           latitude, longitude, elevation_in_meters) %>%
  summarize(n_plants = n_distinct(individual_id), .groups = "keep") %>%
  rename(lat = latitude,
         lon = longitude,
         elev = elevation_in_meters) %>%
  data.frame()

# Map
leaflet(sites) %>% addTiles() %>%
  addCircleMarkers(lng = ~lon, lat = ~lat, radius = ~n_plants, 
                   fillOpacity = 0.6)
# Are sites with no state listed in the US?
mn_range <- list(
  mn = ~mean(.x, na.rm = TRUE),
  min = ~min(.x, na.rm = TRUE), 
  max = ~max(.x, na.rm = TRUE)
)
si %>%
  filter(is.na(state)) %>%
  summarize(across(latitude:longitude, mn_range))
# Seems like it...

# Limit to US + Ontario and longitude > 100 degW
qcsite <- sites$site_id[sites$state == "QC" & !is.na(sites$state)]
sites <- sites %>%
  filter(lon > (-100)) %>%
  filter(site_id != qcsite)
si <- si %>%
  filter(site_id %in% unique(sites$site_id))

# Some sites are missing elevation - will grab elevations for all sites using
# the elevatr package
elev_fill <- sites %>%
  distinct(site_id, lat, lon) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
elev_fill <- get_elev_point(locations = elev_fill, src = "epqs")
elev_fill <- data.frame(elev_fill) %>%
  mutate(elev_new = round(elevation)) %>%
  select(site_id, elev_new)
elevations <- sites %>%
  left_join(elev_fill, by = "site_id")

# Woah! There are a whole bunch of sites in NN database where the elevation
# listed appears to be in feet rather than meters.
elevations <- elevations %>%
  mutate(convert = round(elev/3.28084)) %>%
  mutate(prob = ifelse(abs(elev_new - convert) < abs(elev_new - elev), 1, 0))

# elevations %>%
#   select(site_id, lat, lon, elev, elev_new, prob) %>%
#   rename(elev_original = elev,
#          elev_USGS = elev_new, 
#          problem = prob) %>%
#   write.csv(file = "C:/Users/erin/Desktop/elevation-issues.csv", 
#             row.names = FALSE)

ggplot(elevations, aes(x = elev, y = elev_new)) +
  geom_point(aes(color = factor(prob))) +
  scale_color_manual(values = c("blue", "red")) +
  geom_abline(slope = 1, intercept = 0, color = "blue") +
  geom_abline(slope = 1/3.28084, intercept = 0, color = "red")

count(elevations, prob)
# >11% of sites have elevation reported in feet

# Also, one reported elevation just seems to be wrong (listed as 257 m, should
# be 634 m):
filter(elevations, elev < 400 & elev_new > 600)

# Given all this, we'll use elevations from the R package for all
sites <- sites %>%
  left_join(select(elevations, site_id, elev_new), by = "site_id") %>% 
  select(-elev) %>%
  rename(elev = elev_new)
rm(elevations)

# Combine observations from both flower phenophases ---------------------------#

# Extract all phenophases except "Flowers of flower buds", "Open flowers"
flowers <- filter(si, grepl("flower", phenophase_description)) %>%
  mutate(php = ifelse(phenophase_description == "Open flowers", 
                      "open", "flowers"))

# Taking out site name and other site columns. Can add them back in from sites
# dataframe later.
flowers <- flowers %>%
  select(-c(site_name, latitude, longitude, elevation_in_meters, state, 
            species_id, genus, species, common_name, 
            phenophase_id, phenophase_description, intensity_category_id, 
            class_id, class_name, intensity_name, intensity_type)) %>%
  rename(status = phenophase_status,
         id = individual_id,
         observer = observedby_person_id)

# Put observations of all phenophases by same observer on same day in one row
flowersw <- flowers %>%
  rename(intensity = intensity_midpoint) %>%
  pivot_wider(id_cols = c(observer, site_id, id, observation_date),
              names_from = php,
              values_from = c(status, intensity)) %>%
  data.frame()

# First resolve any status inconsistencies
count(flowersw, status_open, status_flowers)

# Look at instances where open = 1 and flowers = 0 or NA.
# If there's an intensity value for open flowers, then we'll assume that 
# flower status should be 1. Otherwise we'll delete the observation
  # filter(flowersw, status_open == 1) %>%
  #   count(status_flowers, status_open, !is.na(intensity_open))
flowersw <- flowersw %>%
  mutate(status_flowers = case_when(
    status_open == 1 & 
      (is.na(status_flowers) | status_flowers == 0) & 
      !is.na(intensity_open) ~ 1,
    .default = status_flowers
  )) %>%
  filter(!(!is.na(status_open) & status_open == 1 &
             (is.na(status_flowers) | status_flowers == 0)))

# Look at new totals
count(flowersw, status_flowers, status_open, intensity_open)

# Remove observations with status_flowers = NA & status_open = 0 (no info)
flowersw <- flowersw %>%
  filter(!(is.na(status_flowers) & !is.na(status_open) & status_open == 0))
# Now there are no more status_flowers = NA

# If status_flowers == 0, change all status_open to 0 and all intensity values
# to 0
flowersw <- flowersw %>%
  mutate(status_open = ifelse(status_flowers == 0, 0, status_open))

# Convert intensity values to 0 when status = 0 
flowersw <- flowersw %>%
  mutate(intensity_flowers = ifelse(status_flowers == 0, 0, intensity_flowers)) %>%
  mutate(intensity_open = ifelse(!is.na(status_open) & status_open == 0,
                                 0, intensity_open))

# Remove observations when both intensity values are NA
flowersw <- flowersw %>%
  filter(!(is.na(intensity_flowers) & is.na(intensity_open)))

count(flowersw, status_flowers, status_open, intensity_flowers, intensity_open)
# Now, all  status_flowers = 0 or 1 (no NAs)
# There are 237 observations with status_flowers = 1 and status_open = NA
# All intensity values = 0 when status = 0
# NA intensity values only present when status = 1 and the other intensity value is present

# Are there any observations of the same plant on the same day?
flowersw %>% distinct(id, observation_date) %>% nrow() == nrow(flowersw)
# Yes
flowersw %>%
  group_by(id, observation_date) %>%
  summarize(n_obs = n(), 
            n_status_open = sum(!is.na(status_open)),
            n_int_flowers = sum(!is.na(intensity_flowers)),
            n_int_open = sum(!is.na(intensity_open)),
            .groups = "keep") %>%
  data.frame() %>%
  filter(n_obs > 1) 
# In almost every case, (n = 21) one observation provided an intensity values
# for flowers and the other provided an intensity value for open flowers.
# Will combine them.
max_na <- function(x) {
  ifelse(sum(is.na(x)) == length(x), NA, max(x, na.rm = TRUE))
}
flowersw <- flowersw %>%
  group_by(site_id, id, observation_date) %>%
  summarize(status_flowers = max_na(status_flowers),
            status_open = max_na(status_open),
            intensity_flowers = max_na(intensity_flowers),
            intensity_open = max_na(intensity_open),
            .groups = "keep") %>%
  data.frame()

# Add year and doy columns
flowersw <- flowersw %>%
  mutate(yr = year(observation_date),
         doy = yday(observation_date))

# Plot intensity values -------------------------------------------------------#

# # How much data each year?
# flowersw %>%
#   group_by(yr) %>%
#   summarize(n_plants = n_distinct(id),
#             n_sites = n_distinct(site_id),
#             n_obs = n(),
#             n_int_flowers = sum(!is.na(intensity_flowers)),
#             n_int_open = sum(!is.na(intensity_open)),
#             n_int_both = sum(!is.na(intensity_flowers) & !is.na(intensity_open))) %>%
#   data.frame()
# # Way more data in 2022-2025 than previous years
# 
# # Plot mean daily intensity values for 2022 (was going to use median, but 
# # there are too many 0s)
# flowersw %>%
#   filter(yr == 2022) %>%
#   group_by(doy) %>%
#   summarize(flowers = mean(intensity_flowers, na.rm = TRUE),
#             open = mean(intensity_open, na.rm = TRUE)) %>%
#   mutate(log_flowers = ifelse(flowers == 0, log(flowers + 0.1), log(flowers))) %>%
#   select(-flowers) %>%
#   pivot_longer(cols = c(log_flowers, open),
#                names_to = "type",
#                values_to = "intensity") %>%
#   filter(doy %in% 25:200) %>%
#   ggplot(aes(x = doy, y = intensity)) +
#   geom_line() +
#   facet_grid(type ~ ., scales = "free_y")
# # Messy. However, even here it looks like we can have moderate %open values
# # late in the flowering season when the number of flowers decreases 
# # Try weekly?
# flowersw %>%
#   filter(yr == 2022) %>%
#   mutate(wk = as.numeric(week(observation_date))) %>%
#   group_by(wk) %>%
#   summarize(flowers = mean(intensity_flowers, na.rm = TRUE),
#             open = mean(intensity_open, na.rm = TRUE)) %>%
#   mutate(log_flowers = ifelse(flowers == 0, log(flowers + 0.1), log(flowers))) %>%
#   select(-flowers) %>%
#   pivot_longer(cols = c(log_flowers, open),
#                names_to = "type",
#                values_to = "intensity") %>%
#   filter(wk %in% 9:25) %>%
#   ggplot(aes(x = wk, y = intensity)) +
#   geom_line() +
#   facet_grid(type ~ ., scales = "free_y")
# # Looking at weekly values averaged across individuals, this seems like less of 
# # a problem...
# 
# # What about looking at data for an individual?
# count(filter(flowersw, yr == 2022), id) %>% arrange(desc(n)) %>% head(10)
# flowersw %>%
#   filter(yr == 2022 & id == 287224) %>%
#   mutate(flowers_log = ifelse(intensity_flowers == 0, 
#                               log(intensity_flowers + 0.1), 
#                               log(intensity_flowers))) %>%
#   select(-intensity_flowers) %>%
#   pivot_longer(cols = c(flowers_log, intensity_open),
#                names_to = "type",
#                values_to = "intensity") %>%
#   # filter(doy %in% 25:200) %>%
#   ggplot(aes(x = doy, y = intensity)) +
#   geom_line() +
#   facet_grid(type ~ ., scales = "free_y")
# # %open stays high as number of flowers drops off

# Calculate estimated number of open flowers ----------------------------------#

# Remove observations where either intensity value is NA, then calculate
# (rounding values up to nearest integer)
flowersw <- flowersw %>%
  filter(!is.na(intensity_flowers) & !is.na(intensity_open)) %>%
  mutate(nopen = ceiling(intensity_flowers * intensity_open / 100))

# A bunch of filtering (WILL NEED TO REVIEW):
  # Remove observations past DOY 180
  # Remove plant-year combinations with <2 non-zero counts
  # Remove plant-year combinations with no observations before DOY 100
  # Remove plant-year combinations that did not start and end with 0 counts
of <- flowersw %>%
  filter(doy <= 180) %>%
  group_by(id, yr) %>%
  mutate(nonzero_counts = sum(nopen > 0),
         min_doy = min(doy),
         start0 = ifelse(nopen[1] == 0, 1, 0),
         end0 = ifelse(last(nopen) == 0, 1, 0)) %>%
  ungroup() %>% 
  filter(nonzero_counts > 1) %>% 
  filter(min_doy <= 100) %>%
  filter(start0 == 1 & end0 == 1) %>%
  data.frame()

# Merge with site information, and limit data to 2016-2025
of <- of %>%
  # Get rid of some columns we no longer need
  select(-c(nonzero_counts, min_doy, start0, end0)) %>%
  left_join(select(sites, site_id, state, lat, lon, elev), 
            by = "site_id") %>%
  filter(yr > 2015)

# Extract coordinates for remaining sites
peaksites <- distinct(of, site_id, lat, lon)

# Basic summaries/visualizations of open flower abundance data ----------------#

# # What do the values look like?
# count(of, nopen)
# # How much data do we have each year?
# of %>%
#   group_by(yr) %>%
#   summarize(nplants = n_distinct(id),
#             nobs = n()) %>%
#   mutate(nobsperplant = nobs / nplants) %>%
#   data.frame()
# 
# # Look at data from 2016 (first year with >= 10 plants)
# ggplot(filter(of, yr == 2016), aes(x = doy, y = nopen)) +
#   geom_line() +
#   geom_point(color = "blue") +
#   facet_wrap(~id)
# # There are massive differences in max counts among plants within and across
# # years.
# ggplot(filter(of, yr == 2016), aes(x = doy, y = log(nopen + 0.1))) +
#   geom_line() +
#   geom_point(color = "blue") +
#   facet_wrap(~id)
# 
# # Look at variation in max counts among years for a given plant
# of %>%
#   group_by(id) %>%
#   summarize(nyrs = n_distinct(yr)) %>%
#   arrange(desc(nyrs)) %>%
#   head(20)
# 
# plant <- 42681 # 25365
# of %>%
#   filter(id == plant) %>%
#   group_by(yr) %>%
#   summarize(yrmax = max(nopen)) %>%
#   data.frame()
# of %>%
#   filter(id == plant) %>%
#   ggplot(aes(x = doy, y = nopen)) +
#   geom_line() +
#   geom_point(color = "steelblue3") +
#   facet_wrap(~yr)
# # Occasionally, max values increase over time, potentially indicating the plant
# # was small/young in first years (eg, 46192). But most often, max values vary
# # with no trend.

# Create weather covariates for each site -------------------------------------#
# All data from PRISM (temps in degC and precip in mm)
# This section only needs to be run once:

update_weather <- FALSE
weather_missing <- ifelse(length(list.files("weather-data/redbud")) == 0, 
                          TRUE, FALSE)

if(update_weather | weather_missing) {
  
  # First, put site locations to a SpatVector, and convert to NAD83 to match
  # PRISM data
  peaksitesv <- vect(peaksites, geom = c("lon", "lat"), crs = "epsg:4326")
  # Convert to NAD83 to match PRISM data
  peaksitesv <- terra::project(peaksitesv, "epsg:4269")
  
  # Specify years
  yrs <- 2016:2025
  
  for (yr in yrs) {
  
    # Precipitation for winter (Dec-Feb) and 6-mo (Nov-Apr) ---------------------#
    # Get daily precipitation in Nov-Dec of year prior
    rppt <- readRDS(paste0(prism_folder, "ppt-", yr - 1, ".rds"))
    rppt_nov <- rppt[[str_sub(names(rppt), 9, 10) == "11"]]
    ppt_nov <- terra::extract(rppt_nov, peaksitesv, ID = FALSE)
    rppt_dec <- rppt[[str_sub(names(rppt), 9, 10) == "12"]]
    ppt_dec <- terra::extract(rppt_dec, peaksitesv, ID = FALSE)
    
    # Get daily precipitation in Jan-Apr of same year
    rppt <- readRDS(paste0(prism_folder, "ppt-", yr, ".rds"))
    rppt_jan <- rppt[[str_sub(names(rppt), 9, 10) == "01"]]
    ppt_jan <- terra::extract(rppt_jan, peaksitesv, ID = FALSE)
    rppt_feb <- rppt[[str_sub(names(rppt), 9, 10) == "02"]]
    ppt_feb <- terra::extract(rppt_feb, peaksitesv, ID = FALSE)
    rppt_mar <- rppt[[str_sub(names(rppt), 9, 10) == "03"]]
    ppt_mar <- terra::extract(rppt_mar, peaksitesv, ID = FALSE)
    rppt_apr <- rppt[[str_sub(names(rppt), 9, 10) == "04"]]
    ppt_apr <- terra::extract(rppt_apr, peaksitesv, ID = FALSE)
    
    # Get monthly totals
    ppt <- data.frame(site_id = peaksites$site_id,
                      season = yr,
                      nov = apply(ppt_nov, 1, sum),
                      dec = apply(ppt_dec, 1, sum),
                      jan = apply(ppt_jan, 1, sum),
                      feb = apply(ppt_feb, 1, sum),
                      mar = apply(ppt_mar, 1, sum),
                      apr = apply(ppt_apr, 1, sum))
    rm(rppt, rppt_jan, rppt_feb, rppt_mar, rppt_apr, rppt_nov, rppt_dec,
       ppt_jan, ppt_feb, ppt_mar, ppt_apr, ppt_nov, ppt_dec)
    
    # Winter minimum temperatures (Dec-Feb) -------------------------------------#
    # Get daily min temperatures in Dec of year prior
    rtmin <- readRDS(paste0(prism_folder, "tmin-", yr - 1, ".rds"))
    rtmin_dec <- rtmin[[str_sub(names(rtmin), 10, 11) == "12"]]
    tmin_dec <- terra::extract(rtmin_dec, peaksitesv, ID = FALSE)
  
    # Get daily min temperatures in Jan-Feb of same year
    rtmin <- readRDS(paste0(prism_folder, "tmin-", yr, ".rds"))
    rtmin_janfeb <- rtmin[[str_sub(names(rtmin), 10, 11) %in% c("01", "02")]]
    tmin_janfeb <- terra::extract(rtmin_janfeb, peaksitesv, ID = FALSE)
    
    tmin_win <- cbind(site_id = peaksites$site_id, tmin_dec, tmin_janfeb)
    tmin_win <- data.frame(site_id = tmin_win$site_id,
                           season = yr,
                           tmin_win = apply(tmin_win[, -1], 1, mean))
    rm(rtmin_dec, rtmin_janfeb, tmin_dec, tmin_janfeb)
  
    # GDD values (base 0) from Jan 1 or Mar 1 -----------------------------------#
    # Aleady have raster with tmin values; load raster with tmax values
    rtmax <- readRDS(paste0(prism_folder, "tmax-", yr, ".rds"))
    # Stack them, and use the tapp() function to create a SpatRaster with the 
    # mean temp for each day
    rtmean <- c(rtmin, rtmax)
    ndays <- nlyr(rtmax)
    rtmean <- tapp(rtmean, index = c(1:ndays, 1:ndays), fun = "mean")
    names(rtmean) <- str_replace_all(names(rtmax), "tmax", "tmean")
    
    # Extract mean temperatures through end of June
    rtmean <- rtmean[[str_sub(names(rtmean), 11, 12) %in%
                        c("01", "02", "03", "04", "05", "06")]]
    
    # Calculate daily GDD values, starting on Jan 1
    gdd_Jan1 <- terra::extract(rtmean, peaksitesv, ID = FALSE)
    gdd_Jan1[gdd_Jan1 < 0] <- 0
    # Only need through DOY 180
    gdd_Jan1_180 <- gdd_Jan1[,1:180]
    # Will use cumsum(), but note that it works on columns of a dataframe but not 
    # rows so need to do some transposing
    gdd_Jan1_t <- t(gdd_Jan1_180) %>% data.frame()
    agdd_Jan1 <- cumsum(gdd_Jan1_t) %>% t() %>% data.frame()
    names(agdd_Jan1) <- paste0("d", 1:ncol(agdd_Jan1))
    agdd_Jan1 <- cbind(site_id = peaksites$site_id,
                       season = yr,
                       agdd_Jan1)
  
    # Calculate daily GDD values, starting on Mar 1
    gdd_Mar1 <- gdd_Jan1[!str_sub(names(gdd_Jan1), 11, 12) %in% c("01", "02")]
    # To make sure all years have the same number of columns
    gdd_Mar1_t <- t(gdd_Mar1) %>% data.frame()
    agdd_Mar1 <- cumsum(gdd_Mar1_t) %>% t() %>% data.frame()
    names(agdd_Mar1) <- str_replace_all(names(agdd_Mar1), paste0("tmean_", yr), "agdd_")
    agdd_Mar1 <- cbind(site_id = peaksites$site_id,
                       season = yr,
                       agdd_Mar1)
  
    rm(rtmin, rtmax, rtmean, gdd_Jan1, gdd_Jan1_t, gdd_Jan1_180, 
       gdd_Mar1, gdd_Mar1_t)
    
    if (yr == min(yrs)) {
      agdd_jan1_all <- agdd_Jan1
      agdd_mar1_all <- agdd_Mar1
      tmin_win_all <- tmin_win
      ppt_all <- ppt
    } else {
      agdd_jan1_all <- rbind(agdd_jan1_all, agdd_Jan1)
      agdd_mar1_all <- rbind(agdd_mar1_all, agdd_Mar1)
      tmin_win_all <- rbind(tmin_win_all, tmin_win)
      ppt_all <- rbind(ppt_all, ppt)
    }
  }
  
  # Save to file
  write.csv(agdd_jan1_all,
            "weather-data/redbud/agdd-jan1.csv",
            row.names = FALSE)
  write.csv(agdd_mar1_all,
            "weather-data/redbud/agdd-mar1.csv",
            row.names = FALSE)
  write.csv(tmin_win_all,
            "weather-data/redbud/tmin-win.csv",
            row.names = FALSE)
  write.csv(ppt_all,
            "weather-data/redbud/ppt-monthly.csv",
            row.names = FALSE)
}

# Load weather data:
ppt <- read.csv("weather-data/redbud/ppt-monthly.csv")
tmin_win <- read.csv("weather-data/redbud/tmin-win.csv")
agdd_jan1 <- read.csv("weather-data/redbud/agdd-jan1.csv")
agdd_mar1 <- read.csv("weather-data/redbud/agdd-mar1.csv")

# Modeling variation in max annual counts -------------------------------------# 

# Merge with site information, and limit data to 2016-2025
ofmax <- of %>%
  group_by(site_id, id, yr) %>%
  summarize(nobs = n(),
            maxcount = max(nopen),
            .groups = "keep") %>%
  data.frame() %>%
  left_join(select(sites, site_id, state, lat, lon, elev), 
            by = "site_id") %>%
  filter(yr > 2015) %>%
  mutate(fyr = factor(yr))

# Doesn't look like there are obvious differences among years across plants.
# 2020 and 2022 had lower max counts, but they're not hugely different.
ofmax %>%
  group_by(yr) %>%
  summarize(nplants = n_distinct(id),
            nobs = n(),
            of_mean = round(mean(maxcount), 1),
            of_sd = round(sd(maxcount), 1),
            of_med = median(maxcount),
            of_min = min(maxcount),
            of_max = max(maxcount)) %>%
  data.frame()
# Could look for evidence of masting in much smaller area, but it's tough
# as sample size become prohibitively small in early years:
ofmax %>%
  filter(lat > 37 & lat < 40 & lon > (-80) & lon < (-75)) %>%
  group_by(yr) %>%
  summarize(nplants = n_distinct(id),
            nobs = n(),
            of_mean = round(mean(maxcount), 1),
            of_sd = round(sd(maxcount), 1),
            of_med = median(maxcount),
            of_min = min(maxcount),
            of_max = max(maxcount)) %>%
  data.frame()

# Start simple, using lmer with log(maxcount)?
library(lme4)
m1 <- lmer(log(maxcount) ~ lat + lon + elev + (1|fyr) + (1|site_id), 
           data = ofmax)
summary(m1)
# Max counts increased with latitude; possibly a negative effect of longitude
# No effect of elevation. 
# Residual variance (individual) greater than site variance
# Random year effect tiny.

# Would like to get some kind of weather variables (relating to 
# forcing/chilling/precip?)

# Could also bin counts and then use ordinal regression models....


# Evaluating seasonal variation in counts by plant-year -----------------------# 

# See project notes, but I don't think it makes much sense to fit smooths/GAMs
# to the number of open flowers for each plant-year, and there aren't any 
# natural ways to group plants like there were with gardens/accessions in the 
# Rauschkolb preprint.

# Do we need to filter out plant-years where there are long gaps between
# observations? This could really bias estimates of peak date (eg, 150802 in
# 2018). Really, we just need to find those series where the gap is adjacent to
# max values.

of <- of %>%
  arrange(id, observation_date) %>%
  mutate(obsdate = ymd(observation_date)) %>%
  mutate(interval_before = NA) %>%
  select(-observation_date)
for (i in 2:nrow(of)) {
  of$interval_before[i] <- ifelse(
    of$id[i] == of$id[i - 1] & of$yr[i] == of$yr[i - 1],
    as.numeric(of$obsdate[i] - of$obsdate[i - 1]),
    NA
  )
}
nrows <- nrow(of)
of$interval_after <- c(of$interval_before[2:nrows], NA)

of_intermediate <- of %>%
  filter(yr > 2015) %>%
  group_by(yr, id) %>%
  mutate(maxvalue = max(nopen)) %>%
  ungroup() %>%
  filter(nopen == maxvalue) %>%
  group_by(yr, id, maxvalue) %>%
  summarize(first_max = min(doy),
            last_max = max(doy),
            .groups = "keep") %>%
  data.frame()

of_plantyr <- of %>%
  filter(yr > 2015) %>%
  left_join(select(sites, site_id, state, lat, lon, elev), by = "site_id") %>%
  left_join(of_intermediate, by = c("yr", "id")) %>%
  group_by(yr, id, site_id, lat, lon, elev, maxvalue, first_max, last_max) %>%
  summarize(nobs = n(),
            nflower = sum(status_flowers),
            nopen = sum(status_open),
            first = min(doy),
            last = max(doy), 
            interval_mn = mean(interval_before, na.rm = TRUE),
            interval_mx = max(interval_before, na.rm = TRUE),
            interval_beforemax = interval_before[doy == first_max],
            interval_aftermax = interval_after[doy == last_max],
            .groups = "keep") %>%
  data.frame() %>%
  mutate(filter30 = ifelse(interval_beforemax > 30 | interval_aftermax > 30, 1, 0),
         filter21 = ifelse(interval_beforemax > 21 | interval_aftermax > 21, 1, 0))
# There are a total of 525 plant-years in the filtered dataset for 2016-2025

# Look at filtering results
of_plantyr %>%
  group_by(yr) %>%
  summarize(nplants = n(),
            filter30 = n() - sum(filter30),
            filter21 = n() - sum(filter21))

# For now, will use 30-day interval filter
# Calculate peak date as mean date betewen first and last dates with max count
of_plantyr30 <- of_plantyr %>%
  filter(filter30 == 0) %>%
  select(-c(interval_mn, interval_mx, filter30, filter21)) %>%
  mutate(max_span = last_max - first_max,
         peak = floor((first_max + last_max)/2),
         fyr = factor(yr),
         yr0 = yr - min(yr),
         lat_z = (lat - mean(lat)) / sd(lat),
         lon_z = (lon - mean(lon)) / sd(lon),
         elev_z = (elev - mean(elev)) / sd(elev))
 
# Really basic model:
# (note: can't include site id, id, and year REs in the same model. Singular
# fit warnings because there aren't any DF left for residual error)
mpeak1 <- lmer(peak ~ lat_z + lon_z + elev_z + yr0 + (1|fyr) + (1|site_id), 
           data = of_plantyr30)
summary(mpeak1)
# As expected, strong positive association of peak date with latitude
# Marginal, positive elevation effect
# No effect of longitude
# No evidence of a trend from 2016-2025
# Pretty big random year, site, and residual/individual effects 

# Looked at Tariq's paper on timing of flowering onset in eastern redbud.
# Chilling requirement, starting around mid-Oct
# Forcing requirement with ~0 deg base, starting after exceeding "photoperiod" 
# threshold of ~ 11 hours (which is roughly around 1 March for this region)

# Need to get weather data for sites
peaksites <- distinct(of_plantyr30, site_id, lat, lon)


