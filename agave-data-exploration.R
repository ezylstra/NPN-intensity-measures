# Exploring intensity data for Palmer's agave
# ER Zylstra

library(tidyverse)
library(leaflet)
library(elevatr)
library(sf)
library(ordbetareg)   # Loads brms 
library(tidybayes)    # Manipulate Stan objects in a tidy way
# library(broom)        # Convert model objects to data frames
# library(broom.mixed)  # Convert brms model objects to data frames
# library(emmeans)      # Calculate marginal effects in even fancier ways

# Load and format status-intensity data ---------------------------------------#

# List files with formatted intensity data
intensity_files <- list.files("npn-data",
                              pattern = "intensity-agave",
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
  mutate(n_obs = n()) %>%
  ungroup() %>%
  data.frame() %>%
  arrange(individual_id, observation_date, phenophase_id, 
          desc(phenophase_status), desc(intensity_midpoint)) %>%
  distinct(individual_id, phenophase_id, observation_date, .keep_all = TRUE) %>%
  dplyr::select(-n_obs)

# Doublecheck that there's only one observation of each plant-phenophase per day
if (nrow(si) != nrow(distinct(si, individual_id, 
                              phenophase_id, observation_date))) {
  warning("There is more than one observation of some plant phenophases at ",
          sitecap, " in a day")
}

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

# Make intensity midpoints = 0 if status = 0 and add intensity category labels
si <- si %>%
  mutate(intensity_midpoint = ifelse(phenophase_status == 0, 
                                     0, intensity_midpoint)) %>%
  dplyr::select(-c(intensity_name, intensity_type, intensity_label)) %>%
  left_join(spil, by = c("common_name", "phenophase_description"))
# Now the only NAs left in the intensity_midpoint column occur when the status 
# is yes but no intensity value was provided (relative few of these)
# Check:
# count(si, phenophase_description, phenophase_status, intensity_midpoint)

# Filter data: plant-phenophase-year combinations -----------------------------#

sites <- si %>%
  group_by(site_id, site_name, latitude, longitude) %>%
  summarize(n_plants = n_distinct(individual_id), .groups = "keep") %>%
  data.frame()

# Map
leaflet(sites) %>% addTiles() %>%
  addCircleMarkers(lng = ~longitude, lat = ~latitude, radius = ~n_plants, 
                   fillOpacity = 0.6) %>%
  # Sites in Fort Huachuca in Red
  addCircleMarkers(data = filter(sites, 
                                 longitude > -110.4 & longitude < (-110.2),
                                 latitude > 31.48 & latitude < 31.6),
                   lng = ~longitude,
                   lat = ~latitude,
                   fillOpacity = 1, 
                   color = "red",
                   radius = ~n_plants)

# Create indicator for Fort Huachuca sites
sites <- sites %>%
  mutate(fort = ifelse(longitude > -110.4 & longitude < (-110.2) &
                       latitude > 31.48 & latitude < 31.6, 1, 0))

# Will remove plants from CA, TX, and the Phoenix area
sites <- sites %>%
  filter(latitude < 33) %>%
  filter(longitude > (-114) & longitude < (-109))
si <- si %>%
  filter(site_id %in% unique(sites$site_id)) %>%
  left_join(select(sites, site_id, fort), by = "site_id")

# Just going to focus on open flowers for now, so remove everything else
si <- filter(si, phenophase_description == "Open flowers")

# Previous work on Flowers for Bats showed that plants don't flower before 
# day 75, so removing observations before then
si <- filter(si, day_of_year >= 75)

# Calculate the interval between observations of the plant, phenophase
si <- si %>%
  arrange(common_name, individual_id, yr, day_of_year)

si$interval_raw <- c(NA, si$day_of_year[2:nrow(si)] - si$day_of_year[1:(nrow(si) - 1)])
same_ind <- 1 * (si$individual_id[2:nrow(si)] == si$individual_id[1:(nrow(si) - 1)])
same_yr <- 1 * (si$yr[2:nrow(si)] == si$yr[1:(nrow(si) - 1)])
si$same_ind <- c(NA, same_ind)
si$same_yr <- c(NA, same_yr)
si <- si %>%
  mutate(interval = case_when(
    same_ind == 0 ~ NA,
    same_yr == 0 ~ NA,
    is.na(interval_raw) ~ NA,
    .default = interval_raw
  )) %>%
  dplyr::select(-c(interval_raw, same_ind, same_yr))

# Want to:
# remove plant-year combos with no observations in phase
# remove plant-year combos with no observations with intensity values
# remove plant-year combos with < 5 observations between days 100 and 300
# remove plant-year combos when max interval between days 120 and 220 is >21 days
#### May want to revisit this interval length filter

# Summarize amount and quality of information for each plant, year
pl_yr <- si %>%
  arrange(individual_id, observation_date) %>%
  group_by(site_id, individual_id, yr) %>%
  summarize(nobs = n(),
            nobs_100_300 = sum(day_of_year %in% 100:300),
            nobs_120_220 = sum(day_of_year %in% 120:220),
            first_obs = min(day_of_year),
            last_obs = max(day_of_year),
            mean_int = round(mean(interval, na.rm = TRUE), 2),
            max_int = ifelse(nobs == 1, NA, max(interval, na.rm = TRUE)),
            n_inphase = sum(phenophase_status),
            n_intvalue = sum(!is.na(intensity_value)),
            prop_inphase = round(n_inphase / nobs, 2),
            prop_intvalue = round(n_intvalue / n_inphase, 2),
            .groups = "keep") %>%
  data.frame()

# Evaluating interval lengths between days 120 and 220 (need to use days 140 to
# 240 since interval is the length of time since previous observation)
pl_yr2 <- si %>%
  arrange(individual_id, observation_date) %>%
  filter(day_of_year %in% 140:240) %>%
  group_by(site_id, individual_id, yr) %>%
  summarize(nobs = n(),
            max_int_after120 = ifelse(nobs == 1 & anyNA(interval), NA, 
                                      max(interval, na.rm = TRUE)),
            .groups = "keep") %>%
  data.frame()

pl_yr <- pl_yr %>%
  left_join(select(pl_yr2, -nobs), by = c("site_id", "individual_id", "yr"))

# Try 2 different filters...
pl_yr <- pl_yr %>%
  mutate(remove1 = case_when(
    n_inphase == 0 ~ 1,
    n_intvalue == 0 ~ 1,
    nobs_100_300 < 5 ~ 1,
    max_int_after120 > 21 ~ 1,
    .default = 0
  )) %>%
  mutate(remove2 = case_when(
    n_inphase == 0 ~ 1,
    n_intvalue == 0 ~ 1,
    nobs_100_300 < 5 ~ 1,
    nobs_120_220 < 5 ~ 1,
    # max_int_after120 > 21 ~ 1,
    .default = 0
  ))

si <- si %>%
  left_join(select(pl_yr, individual_id, yr, remove1, remove2),
            by = c("individual_id", "yr")) 

# What's left if I use each filter?
si %>%
  filter(remove1 == 0) %>%
  group_by(yr) %>%
  summarize(n_obs = n(),
            n_sites = n_distinct(site_id),
            n_plants = n_distinct(individual_id),
            obs_per_plant = n_obs / n_plants) %>%
  data.frame()
si %>%
  filter(remove2 == 0) %>%
  group_by(yr) %>%
  summarize(n_obs = n(),
            n_sites = n_distinct(site_id),
            n_plants = n_distinct(individual_id),
            obs_per_plant = n_obs / n_plants) %>%
  data.frame()

# Remove data from 2017 (just one plant; all other years 17-55 plants)
si <- filter(si, yr != 2017)

# Extract information about sites ---------------------------------------------#

sites <- si %>%
  filter(remove1 == 0 | remove2 == 0) %>%
  distinct(latitude, longitude, site_id)
# write.table(sites, "weather-data/agave-sites.csv", sep = ",",
#             row.names = FALSE,
#             col.names = FALSE)

# Plot raw intensity data -----------------------------------------------------#
# Using 2nd filter for now

# Create filename for png
png_name <- "output/IntensityData-OpenFlowers-agave.png" 

# Create ggplot object and save to file (if file doesn't already exist)
if (!file.exists(png_name)) {

  # Figures with lines and points, with size proportional to no. of plants 
  agg <- si %>%
    filter(remove2 == 0) %>%
    group_by(yr, day_of_year, intensity_midpoint) %>%
    summarize(n_indiv = n(), .groups = "keep") %>%
    data.frame()
  
  iplot <- ggplot(agg, aes(x = day_of_year, y = intensity_midpoint)) +
    facet_grid(yr ~ .) +
    geom_line(data = filter(si, remove2 == 0), 
              show.legend = FALSE, alpha = 0.5,
              aes(color = factor(individual_id))) +
    geom_point(data = agg, aes(size = n_indiv), shape = 16, alpha = 0.4) +
    labs(title = "Palmer's century plant, open flowers",
         y = "Open flowers (%)", x = "Day of year") +
    theme_bw()
  iplot
  
  ggsave(png_name,
         plot = iplot,
         width = 6.5,
         height = 9,
         units = "in")
}

# Look at individual curves for a given year
# si %>%
#   filter(remove2 == 0) %>%
#   filter(yr == 2023) %>%
#   ggplot(aes(x = day_of_year, y = intensity_midpoint)) +
#   geom_line() +
#   facet_wrap(~ individual_id) +
#   theme_bw()

# Not sure how I want to deal with the data after peak intensity value, since I 
# really just want to model the increase.....
  # Make all intensity values = NA when status goes back down to zero?

# Prepping data for ordbetareg models -----------------------------------------#

# First, remove plant years that would be excluded with both filters
sif <- si %>%
  filter(remove1 == 0 | remove2 == 0)

# For each plant and year, will find last instance with status = 1 and make all 
# intensity values after that date = NA
sif <- sif %>%
  group_by(individual_id, yr) %>%
  mutate(last_yes = max(day_of_year[phenophase_status == 1])) %>%
  ungroup() %>%
  data.frame()
sif <- sif %>%
  mutate(intensity = ifelse(day_of_year > last_yes, NA, intensity_midpoint))

# Check:
# sif %>%
#   filter(remove2 == 0) %>%
#   filter(yr == 2023) %>%
#   ggplot(aes(x = day_of_year, y = intensity)) +
#   geom_line() +
#   facet_wrap(~ individual_id) +
#   theme_bw()

# Convert intensities to proportions and remove NAs
sif <- sif %>%
  mutate(prop = intensity / 95) %>%
  filter(!is.na(prop))

# Fill in missing elevation data ----------------------------------------------#

elev_fill <- filter(sif, is.na(elevation_in_meters)) %>%
  distinct(site_id, latitude, longitude) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) 
elev_fill <- get_elev_point(locations = elev_fill, src = "epqs")
elev_fill <- data.frame(elev_fill) %>%
  mutate(elev = round(elevation)) %>%
  select(site_id, elev)
sif <- sif %>%
  left_join(elev_fill, by = "site_id") %>%
  mutate(elev = ifelse(!is.na(elev), elev, elevation_in_meters)) %>%
  select(-elevation_in_meters)

# Formatting weather data -----------------------------------------------------#

# Weather data was obtained on this website (much faster than using API or prism
# package for a limited number of sites at high spatial, temporal resolution):
# https://prism.oregonstate.edu/explorer/bulk.php

prism_folder <- "weather-data/agaves-prism/"

# Daily data
prism_files <- list.files(prism_folder,
                          full.names = TRUE,
                          pattern = "stable")
for (i in 1:length(prism_files)) {
  prism1 <- read.csv(prism_files[i],
                     header = FALSE,
                     skip = 11,
                     col.names = c("site_id", "lon", "lat", "elev",
                                   "date", "ppt", "tmin", "tmean", "tmax"))
  if (i == 1) {
    weather <- prism1
  } else {
    weather <- rbind(weather, prism1)
  }
}
rm(prism1)

# Try a couple simple variables? Generally looks like plants start flowering
# in June.
# Daily variables to explain seasonal curve
  # GDD from start of year with 0 degC base
  # GDD from start of year with 50 degF (10 degC) base
# Yearly variables to explain inter-annual differences
  # Cumulative precipitation for previous 3, 6, 9, or 12 months
  # Mean spring temperatures
  # Winter min temperatures

# Calculate daily GDD values (gdd = tmean - base)
weather <- weather %>%
  mutate(gdd0 = ifelse(tmean < 0, 0, tmean - 0)) %>%
  mutate(gdd10 = ifelse(tmean < 10, 0, tmean - 10)) %>%
  mutate(yr = year(date),
         doy = yday(date)) %>%
  arrange(site_id, date)

# Calculate AGDD values
weather <- weather %>%
  group_by(site_id, yr) %>%
  mutate(agdd0 = cumsum(gdd0),
         agdd10 = cumsum(gdd10)) %>%
  ungroup() %>%
  data.frame()

# Calculate seasonal variables
weathermonths <- weather %>%
  mutate(month = month(date)) %>%
  # Create seasonyr to match up with flowering year
  mutate(seasonyr = ifelse(month %in% 6:12, yr + 1, yr)) %>%
  mutate(season = case_when(
    month %in% 3:5 ~ "sp",
    month %in% 6:8 ~ "su",
    month %in% 9:11 ~ "fa",
    .default = "wi"
  )) %>%
  group_by(site_id, seasonyr) %>%
  summarize(tmin_wi = mean(tmin[season == "wi"]),
            ppt_3m = sum(ppt[season == "sp"]),
            ppt_6m = sum(ppt[season %in% c("wi", "sp")]),
            ppt_9m = sum(ppt[season != "su"]),
            ppt_12m = sum(ppt),
            temp_sp = mean(tmean[season == "sp"]),
            .groups = "keep") %>%
  filter(seasonyr %in% unique(sif$yr)) %>%
  data.frame()

# Attach weather data to phenology data
sif <- sif %>%
  left_join(select(weather, site_id, date, agdd0, agdd10), 
            by = c("observation_date" = "date",
                   "site_id" = "site_id")) %>%
  left_join(weathermonths, by = c("yr" = "seasonyr", "site_id" = "site_id"))

# Run model for one year ------------------------------------------------------#

# Use one of the filters, create factor variables, and standardize DOY
sif2_2023 <- sif %>%
  filter(remove2 == 0) %>%
  filter(yr == 2023) %>%
  mutate(doyz = (day_of_year - mean(day_of_year)) / sd(day_of_year)) %>%
  mutate(id = factor(individual_id),
         # fyr = factor(yr),
         site = factor(site_id))

#### Basic DOY model
m_2023 <- ordbetareg(prop ~ doyz + (1|id),
                     data = sif2_2023,
                     control = list(adapt_delta = 0.9),
                     cores = 4, chains = 4)
summary(m_2023)

# For a good explanation of different predictions types, see:
# https://www.andrewheiss.com/blog/2021/11/10/ame-bayes-re-guide/

# Expected proportion open by date, ignoring individual effects (ie, grand mean)
doy_gm <- m_2023 %>%
  epred_draws(newdata = data.frame(doyz = seq(min(sif2_2023$doyz),
                                              max(sif2_2023$doyz),
                                              length = 100)),
              re_formula = NA)
plot_doy_gm <- ggplot(doy_gm, 
                      aes(x = doyz, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Day of year", y = "Predicted proportion of flowers open (%)",
       fill = "Credible interval", title = "Ignoring individual effects") +
  theme_bw() +
  theme(legend.position = "bottom")
plot_doy_gm

# Conditional effects for a new, typical individual (new individual is based on
# variation among existing individuals)
doy_typicalplant <- m_2023 %>%
  epred_draws(newdata = data.frame(doyz = seq(min(sif2_2023$doyz),
                                              max(sif2_2023$doyz),
                                              length = 100),
                                   individual_id = "new plant"),
              re_formula = NULL,
              allow_new_levels = TRUE,
              sample_new_levels = "uncertainty")
plot_doy_typicalplant <- ggplot(doy_typicalplant, 
                                aes(x = doyz, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Day of year", y = "Predicted proportion of flowers open (%)",
       fill = "Credible interval", title = "New, typical plant") +
  theme_bw() +
  theme(legend.position = "bottom")
plot_doy_typicalplant

# Conditional effects for a brand new individual (new individual is based on
# random draws from model)
doy_newplant <- m_2023 %>%
  epred_draws(newdata = data.frame(doyz = seq(min(sif2_2023$doyz),
                                              max(sif2_2023$doyz),
                                              length = 100),
                                   individual_id = "new plant"),
              re_formula = NULL,
              allow_new_levels = TRUE,
              sample_new_levels = "gaussian")
plot_doy_newplant <- ggplot(doy_newplant, 
                            aes(x = doyz, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Day of year", y = "Predicted proportion of flowers open (%)",
       fill = "Credible interval", title = "Brand new plant") +
  theme_bw() +
  theme(legend.position = "bottom")
plot_doy_newplant

#### Run model with quadratic term for DOY
start_time <- Sys.time()
m_2023q <- ordbetareg(prop ~ doyz + I(doyz^2) + (1|id),
                      data = sif2_2023,
                      control = list(adapt_delta = 0.9),
                      cores = 4, chains = 4)
end_time <- Sys.time()
end_time - start_time
summary(m_2023q)

# Expected proportion open by date, ignoring individual effects (ie, grand mean)
doyq_gm <- m_2023q %>%
  epred_draws(newdata = data.frame(doyz = seq(min(sif2_2023$doyz),
                                              max(sif2_2023$doyz),
                                              length = 100)),
              re_formula = NA)
plot_doyq_gm <- ggplot(doyq_gm, 
                       aes(x = doyz, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Day of year", y = "Predicted proportion of flowers open (%)",
       fill = "Credible interval", 
       title = "Quadatic model, ignoring individual effects") +
  theme_bw() +
  theme(legend.position = "bottom")
plot_doyq_gm

#### Run model with AGDD, base 0
sif2_2023 <- sif2_2023 %>%
  mutate(agdd0z = (agdd0 - mean(agdd0)) / sd(agdd0))

start_time <- Sys.time()
m_2023gdd0 <- ordbetareg(prop ~ agdd0z + I(agdd0z^2) + (1|id),
                         data = sif2_2023,
                         control = list(adapt_delta = 0.9),
                         cores = 4, chains = 4)
end_time <- Sys.time()
end_time - start_time
summary(m_2023gdd0)

#### Run model with AGDD, base 10
sif2_2023 <- sif2_2023 %>%
  mutate(agdd10z = (agdd10 - mean(agdd10)) / sd(agdd10))

start_time <- Sys.time()
m_2023gdd10 <- ordbetareg(prop ~ agdd10z + I(agdd10z^2) + (1|id),
                          data = sif2_2023,
                          control = list(adapt_delta = 0.9),
                          cores = 4, chains = 4)
end_time <- Sys.time()
end_time - start_time
summary(m_2023gdd10)

loo_1 <- loo(m_2023)
loo_q <- loo(m_2023q)
loo_gdd0 <- loo(m_2023gdd0)
loo_gdd10 <- loo(m_2023gdd10)
loo_1
loo_q
loo_gdd0
loo_gdd10
loo_compare(loo_1, loo_q, loo_gdd0, loo_gdd10)
# LOOCV suggests the quadratic DOY model is the best fit to data, followed 
# closely by GDD:base0 and then GDD:base10. 
# All models took less than 5 minutes to run.

#### Model with elevation and quadratic DOY effect
sif2_2023 <- sif2_2023 %>%
  mutate(elevz = (elev - mean(elev)) / sd(elev))

start_time <- Sys.time()
m_2023qe <- ordbetareg(prop ~ doyz + I(doyz^2) + elevz + (1|id),
                       data = sif2_2023,
                       control = list(adapt_delta = 0.9),
                       cores = 4, chains = 4)
end_time <- Sys.time()
end_time - start_time
summary(m_2023qe)

loo_qe <- loo(m_2023qe)
loo_compare(loo_qe, loo_q)
# Elevation doesn't seem to contribute anything. Credible interval widely 
# overlaps zeroo and LOOCV suggests the model is no better than the model 
# without elevation.

# Run models for multiple years -----------------------------------------------#

# Use one of the filters, create factor variables, and standardize variables
sif2 <- sif %>%
  filter(remove2 == 0) %>%
  mutate(doyz = (day_of_year - mean(day_of_year)) / sd(day_of_year),
         doyz2 = doyz * doyz,
         agdd0z = (agdd0 - mean(agdd0)) / sd(agdd0),
         agdd0z2 = agdd0z * agdd0z,
         agdd10z = (agdd10 - mean(agdd10)) / sd(agdd10),
         tmin_wiz = (tmin_wi - mean(tmin_wi)) / sd(tmin_wi),
         ppt_3mz = (ppt_3m - mean(ppt_3m)) / sd(ppt_3m),
         ppt_6mz = (ppt_6m - mean(ppt_6m)) / sd(ppt_6m),
         ppt_9mz = (ppt_9m - mean(ppt_9m)) / sd(ppt_9m),
         ppt_12mz = (ppt_12m - mean(ppt_12m)) / sd(ppt_12m),
         temp_spz = (temp_sp - mean(temp_sp)) / sd(temp_sp)) %>%
  mutate(id = factor(individual_id),
         fyr = factor(yr),
         site = factor(site_id))

# DOY model with spring (3-month) precipitation
start_doyp3 <- Sys.time()
m_doyp3 <- ordbetareg(prop ~ doyz + I(doyz^2) + ppt_3mz + (1|id),
                      data = sif2,
                      control = list(adapt_delta = 0.99),
                      iter = 4000, cores = 4, chains = 4)
end_doyp3 <- Sys.time()
save(m_doyp3, file = "output/agave-multiyr-models/doyp3.RData")
end_doyp3 - start_doyp3

# DOY model with 6-month precipitation
start_doyp6 <- Sys.time()
m_doyp6 <- ordbetareg(prop ~ doyz + I(doyz^2) + ppt_6mz + (1|id),
                      data = sif2,
                      control = list(adapt_delta = 0.99),
                      iter = 4000, cores = 4, chains = 4)
end_doyp6 <- Sys.time()
save(m_doyp6, file = "output/agave-multiyr-models/doyp6.RData")
end_doyp6 - start_doyp6

# DOY model with 9-month precipitation
start_doyp9 <- Sys.time()
m_doyp9 <- ordbetareg(prop ~ doyz + I(doyz^2) + ppt_9mz + (1|id),
                      data = sif2,
                      control = list(adapt_delta = 0.99),
                      iter = 4000, cores = 4, chains = 4)
end_doyp9 <- Sys.time()
save(m_doyp9, file = "output/agave-multiyr-models/doyp9.RData")
end_doyp9 - start_doyp9

# DOY model with 12-month precipitation
start_doyp12 <- Sys.time()
m_doyp12 <- ordbetareg(prop ~ doyz + I(doyz^2) + ppt_12mz + (1|id),
                       data = sif2,
                       control = list(adapt_delta = 0.99),
                       iter = 4000, cores = 4, chains = 4)
end_doyp12 <- Sys.time()
save(m_doyp12, file = "output/agave-multiyr-models/doyp12.RData")
end_doyp12 - start_doyp12

# DOY model with spring temperatures
start_doyspt <- Sys.time()
m_doyspt <- ordbetareg(prop ~ doyz + I(doyz^2) + temp_spz + (1|id),
                       data = sif2,
                       control = list(adapt_delta = 0.99),
                       iter = 4000, cores = 4, chains = 4)
end_doyspt <- Sys.time()
save(m_doyspt, file = "output/agave-multiyr-models/doyspt.RData")
end_doyspt - start_doyspt

# DOY model with winter temperatures
start_doywit <- Sys.time()
m_doywit <- ordbetareg(prop ~ doyz + I(doyz^2) + tmin_wiz + (1|id),
                       data = sif2,
                       control = list(adapt_delta = 0.99),
                       iter = 4000, cores = 4, chains = 4)
end_doywit <- Sys.time()
save(m_doywit, file = "output/agave-multiyr-models/doywit.RData")
end_doywit - start_doywit

# DOY model with random year effects
start_doyfyr <- Sys.time()
m_doyfyr <- ordbetareg(prop ~ doyz + doyz2 + (1 + doyz + doyz2|fyr) + (1|id),
                       data = sif2,
                       control = list(adapt_delta = 0.99),
                       iter = 4000, cores = 4, chains = 4)
end_doyfyr <- Sys.time()
save(m_doyfyr, file = "output/agave-multiyr-models/doyfyr.RData")
end_doyfyr - start_doyfyr

# GDD model, assuming effects of gdd are the same each year
start_gdd <- Sys.time()
m_gdd <- ordbetareg(prop ~ agdd0z + I(agdd0z^2) + (1|id),
                    data = sif2,
                    control = list(adapt_delta = 0.99),
                    iter = 4000, cores = 4, chains = 4)
end_gdd <- Sys.time()
save(m_gdd, file = "output/agave-multiyr-models/gdd.RData")
end_gdd - start_gdd

# GDD model, random yearly intercepts
start_gddREi <- Sys.time()
m_gddREi <- ordbetareg(prop ~ agdd0z + I(agdd0z^2) + (1|fyr) + (1|id),
                       data = sif2,
                       control = list(adapt_delta = 0.99),
                       iter = 4000, cores = 4, chains = 4)
end_gddREi <- Sys.time()
save(m_gddREi, file = "output/agave-multiyr-models/gddREi.RData")
end_gddREi - start_gddREi

# GDD model, random yearly intercepts and slopes
start_gddREis <- Sys.time()
m_gddREis <- ordbetareg(prop ~ agdd0z + agdd0z2 + (1 + agdd0z + agdd0z2|fyr) + (1|id),
                        data = sif2,
                        control = list(adapt_delta = 0.99),
                        iter = 4000, cores = 4, chains = 4)
end_gddREis <- Sys.time()
save(m_gddREis, file = "output/agave-multiyr-models/gddREis.RData")
end_gddREis - start_gddREis

# Compare multi-year models
loo_doyp3 <- loo(m_doyp3)
loo_doyp6 <- loo(m_doyp6)
loo_doyp9 <- loo(m_doyp9)
loo_doyp12 <- loo(m_doyp12)
loo_doyspt <- loo(m_doyspt)
loo_doywit <- loo(m_doywit)
loo_doyfyr <- loo(m_doyfyr)
loo_gdd <- loo(m_gdd)
loo_gddREi <- loo(m_gddREi)
loo_gddREis <- loo(m_gddREis)
loo_compare(loo_doyp3, loo_doyp6, loo_doyp9, loo_doyp12, 
            loo_doyspt, loo_doywit, loo_doyfyr,
            loo_gdd, loo_gddREi, loo_gddREis)
# DOY with random slopes and intercepts by year is best, followed by GDD model 
# with random slopes and intercepts by year.


# Visualize results from highest-ranked model

# Expected proportion open by year, ignoring individual REs
doyz <- seq(min(sif2$doyz), max(sif2$doyz), length = 100)

newdata <- expand_grid(doyz = doyz,
                       fyr = levels(sif2$fyr)) %>%
  mutate(doyz2 = doyz * doyz) %>%
  mutate(doy = round(doyz * sd(sif2$day_of_year) + mean(sif2$day_of_year))) %>%
  data.frame()

doy_yr <- m_doyfyr %>%
  epred_draws(newdata = newdata,
              re_formula = ~ (1 + doyz + doyz2|fyr)) 
plot_doy_yr <- ggplot(doy_yr, 
                      aes(x = doy, y = .epred)) +
  stat_lineribbon(aes(color = fyr, fill = after_scale(alpha(color, 0.2))), 
                  .width = 0.80) +
  labs(x = "Day of year", y = "Predicted proportion of flowers open (%)", 
       color = "Year",
       fill = "Year", 
       title = "By year, correlated intercepts/slopes") +
  theme_bw() +
  theme(legend.position = "bottom")
plot_doy_yr
# ggsave("output/agave-multiyr-models/doy-fyr.png",
#        plot_doy_yr,
#        width = 6.5,
#        height = 4,
#        units = "in")

plot_doy_yr2 <- ggplot(filter(doy_yr, doy >= 175), 
                       aes(x = doy, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  # Add a horizontal line for easier comparison
  geom_hline(yintercept = 0.25, linetype = "dotted") +
  facet_wrap(~ fyr, ncol = 2) +
  labs(x = "Day of year", y = "Predicted proportion of flowers open (%)", 
       fill = "Credible interval", 
       title = "By year, correlated intercepts/slopes") +
  theme_bw() +
  theme(legend.position = "bottom")
plot_doy_yr2
# ggsave("output/agave-multiyr-models/doy-fyr-faceted.png",
#        plot_doy_yr2,
#        width = 6.5,
#        height = 6,
#        units = "in")
