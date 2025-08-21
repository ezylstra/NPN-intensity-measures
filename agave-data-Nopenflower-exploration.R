# Exploring intensity data for Palmer's agave
# ER Zylstra

library(tidyverse)
library(leaflet)
library(elevatr)
library(sf)
library(ordbetareg)   # Loads brms 
library(tidybayes)    # Manipulate Stan objects in a tidy way

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
if (sum(si$amended_status) == 0) {
  si$amended_status <- NULL
}

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

# Combine observations from both flower phenophases ---------------------------#

# Remove all phenophases except "Flowers of flower buds", "Open flowers"
si <- filter(si, grepl("flower", phenophase_description)) %>%
  mutate(php = ifelse(phenophase_description == "Open flowers", 
                      "open", "flowers"))

# Taking out site name and other site columns. Can add them back in from sites
# dataframe later.
sib <- si %>%
  select(-c(site_name, latitude, longitude, elevation_in_meters, state, 
            species_id, genus, species, common_name, species_functional_type,
            phenophase_id, phenophase_description, intensity_category_id, 
            class_id, class_name, intensity_name, intensity_type)) %>%
  rename(status = phenophase_status,
         id = individual_id,
         observer = observedby_person_id)

# Put observations of all phenophases by same observer on same day in one row
sibw <- sib %>%
  rename(intensity = intensity_midpoint) %>%
  pivot_wider(id_cols = c(observer, site_id, id, observation_date),
              names_from = php,
              values_from = c(status, intensity)) %>%
  data.frame()

# First resolve any status inconsistencies
count(sibw, status_flowers, status_open)
# look at instances where open = 1 and flowers = 0 or NA
filter(sibw, status_open == 1 & (is.na(status_flowers) | status_flowers == 0))
# If there's an intensity value for open flowers, then we'll assume that 
# flower status should be 1. Otherwise we'll delete the observation
sibw <- sibw %>%
  mutate(status_flowers = case_when(
    status_open == 1 & 
      (is.na(status_flowers) | status_flowers == 0) & 
      !is.na(intensity_open) ~ 1,
    .default = status_flowers
  )) %>%
  filter(!(!is.na(status_flowers) & !is.na(status_open) & 
           status_flowers == 0 & status_open == 1))

# Check for any intensity value inconsistencies 
# (right now, intensity = NA whenever status = 0)
count(sibw, status_flowers, intensity_flowers)
count(sibw, status_open, intensity_open)

# Convert intensity values to 0 when status = 0 and then remove any rows that
# have either intensity_value = NA
sibw <- sibw %>%
  mutate(intensity_open = ifelse(status_open == 0, 0, intensity_open),
         intensity_flowers = ifelse(status_flowers == 0, 0, intensity_flowers)) %>%
  filter(!is.na(intensity_open) & !is.na(intensity_flowers))
# Check:
# count(sibw, status_flowers, status_open, 
#       intensity_flowers > 0, intensity_open > 0)

# Are there any observations of the same plant on the same day?
sibw %>% distinct(id, observation_date) %>% nrow() == nrow(sibw)
# Nope

# Calculate the interval between observations of the plant
sibw <- sibw %>%
  mutate(yr = year(observation_date),
         doy = yday(observation_date)) %>%
  arrange(id, yr, doy)

sibw$interval_raw <- c(NA, sibw$doy[2:nrow(sibw)] - sibw$doy[1:(nrow(sibw) - 1)])
same_ind <- 1 * (sibw$id[2:nrow(sibw)] == sibw$id[1:(nrow(sibw) - 1)])
same_yr <- 1 * (sibw$yr[2:nrow(sibw)] == sibw$yr[1:(nrow(sibw) - 1)])
sibw$same_ind <- c(NA, same_ind)
sibw$same_yr <- c(NA, same_yr)
sibw <- sibw %>%
  mutate(interval = case_when(
    same_ind == 0 ~ NA,
    same_yr == 0 ~ NA,
    is.na(interval_raw) ~ NA,
    .default = interval_raw
  )) %>%
  dplyr::select(-c(interval_raw, same_ind, same_yr))

# What's the distribution of observation intervals?
hist(sibw$interval, breaks = 50)
head(count(sibw, interval), 25)
# Big drop off after 21 days

# Filter the data -------------------------------------------------------------#

# First, remove observations before doy 75 (based on Flowers for Bats analysis)
sibw <- filter(sibw, doy >= 75)

# Then, want to:
# remove plant-year combos with no observations in open flower phase
# remove plant-year combos with < 3 observations between days 100 and 300
# remove plant-year combos with a max interval of > 21 days? 30 days?

# Summarize amount and quality of information for each plant, year
pl_yr <- sibw %>%
  arrange(id, observation_date) %>%
  group_by(site_id, id, yr) %>%
  summarize(nobs = n(),
            nobs_100_300 = sum(doy %in% 100:300),
            first_obs = min(doy),
            last_obs = max(doy),
            mean_int = round(mean(interval, na.rm = TRUE), 2),
            max_int = ifelse(nobs == 1, NA, max(interval, na.rm = TRUE)),
            n_inphase = sum(status_open),
            prop_inphase = round(n_inphase / nobs, 2),
            .groups = "keep") %>%
  data.frame()
pl_yr <- pl_yr %>%
  mutate(inphase = ifelse(n_inphase == 0, 0, 1),
         minobs3 = ifelse(nobs_100_300 < 3, 0, 1),
         maxint21 = ifelse(max_int > 21, 0, 1),
         maxint30 = ifelse(max_int > 30, 0, 1))

count(pl_yr, inphase, minobs3, maxint21)
# How much data would we have if we exlcuded those with max interval > 21 days?
pl_yr %>%
  filter(inphase + minobs3 + maxint21 == 3) %>%
  group_by(yr) %>%
  summarize(nplants = n_distinct(id))
# How much data if we relaxed the max interval length to 30 days?
pl_yr %>%
  filter(inphase + minobs3 + maxint30 == 3) %>%
  group_by(yr) %>%
  summarize(nplants = n_distinct(id))
# Adding 31 plants in 2018, 10 plants in 2019, and 17 plants in 2021

# Remove data from 2017 (just one plant; all other years a lot more)
sibw <- filter(sibw, yr != 2017)

# Attach filtering columns to wide dataframe
sibw <- sibw %>%
  left_join(select(pl_yr, id, yr, inphase, minobs3, maxint21, maxint30), 
            by = c("id", "yr"))

# Calculate number of open flowers --------------------------------------------#

# First pick filter: 21 or 30 day max interval length?
# Picking 30-day for now #######
df <- sibw %>%
  filter(inphase + minobs3 + maxint30 == 3)

# Then, calculate number of open flowers
df <- df %>%
  mutate(nopen = intensity_flowers * intensity_open / 100) %>%
  mutate(nopen = ifelse(nopen > 0 & nopen <= 1, 1, nopen)) %>%
  mutate(nopen = round(nopen))

count(df, nopen)
# 1 "count" = 700 and 3 = 1850; all others < 476

# Unique values:
# 0, 1-10, 18, 31, 42, 48, 70, 185, 310, 420, 475, 700, 1850

bins7 <- c(0, 1, 5, 10, 50, 100, 500, 2000)
table(cut(df$nopen[df$nopen > 0], bins7))

bins4 <- c(0, 1, 10, 100, 2000)
table(cut(df$nopen[df$nopen > 0], bins4))

# What do max counts look like?
dfmax <- df %>%
  group_by(id, yr) %>%
  summarize(maxopen = max(nopen), .groups = "keep") %>%
  data.frame()

count(dfmax, maxopen)

# Add bins (n = 4) for max annual counts
dfmax <- dfmax %>%
  mutate(max_bin = case_when(
    maxopen == 1 ~ "one",
    maxopen %in% 2:10 ~ "few",
    maxopen %in% 11:100 ~ "some",
    maxopen > 100 ~ "many"
  )) %>%
  mutate(max_bin = factor(max_bin,
                          levels = c("one", "few", "some", "many"),
                          ordered = TRUE))

# Plot counts in each category by year
ggplot(dfmax, aes(x = max_bin)) +
  geom_bar() + 
  facet_grid(yr ~ .)
# Why so many max counts of 1 in 2023 and 2024?

# Do we see the same pattern in flower counts?
df %>%
  filter(intensity_flowers > 0) %>%
  group_by(id, yr) %>%
  summarize(maxflower = max(intensity_flowers)) %>%
  ggplot(aes(x = factor(maxflower))) +
  geom_bar() +
  facet_grid(yr ~ .)
# Yes...

# After filtering out plant-years when the plants never had open flowers, then
# the vast majority of plants should only be included in the dataset in one year.
# The exception would be if they had more than one shoot from the same base
# (Also wondering if some people used the same ID# for more than one plant in 
# close proximity)
# Track years each plant was observed
plants <- df %>%
  group_by(site_id, id) %>%
  summarize(y2018 = ifelse(2018 %in% yr, 1, 0),
            y2019 = ifelse(2019 %in% yr, 1, 0),
            y2020 = ifelse(2020 %in% yr, 1, 0),
            y2021 = ifelse(2021 %in% yr, 1, 0),
            y2022 = ifelse(2022 %in% yr, 1, 0),
            y2023 = ifelse(2023 %in% yr, 1, 0),
            y2024 = ifelse(2024 %in% yr, 1, 0),
            nyrs = n_distinct(yr)) %>%
  data.frame()
count(plants, nyrs) 
# 13 of 206 plants had open flower counts in multiple years

# Were different areas monitored each year?
df %>%
  distinct(site_id, id, yr) %>%
  left_join(select(sites, site_id, latitude, longitude), by = "site_id") %>%
  group_by(yr) %>%
  summarize(nplants = n_distinct(id),
            nsites = n_distinct(site_id),
            minlat = min(latitude),
            maxlat = max(latitude),
            minlon = min(longitude),
            maxlon = max(longitude)) %>%
  data.frame()
# Fewer sites monitored in 2022-2024 (8-11 versus 17-61 in other years)

# Plot locations each year:
df %>%
  group_by(site_id, yr) %>%
  summarize(nplants = n_distinct(id)) %>%
  left_join(select(sites, site_id, latitude, longitude), by = "site_id") %>%
  ggplot(aes(x = longitude, y = latitude)) +
  geom_jitter(aes(size = nplants)) +
  facet_wrap(yr ~ .) +
  theme_bw()
# There were areas in Chiricahua NM and Fort Bowie NHS where lots of plots were
# observed in 2019, 2023, and 2024. None of these areas were surveyed in other 
# years.

# Was the sampling freq in those years also different?
pl_yr %>%
  group_by(yr) %>%
  summarize(nobs_mn = round(mean(nobs)),
            nobs_min = min(nobs),
            nobs_max = max(nobs),
            firstobs_mn = mean(first_obs),
            lastobs_mn = mean(last_obs))
# Not really.

# Look at seasonal trajectories 
# First log values...
df <- df %>%
  mutate(nopen_log = ifelse(nopen == 0, NA, log(nopen)),
         nopen_log0 = ifelse(nopen == 0, log(0.5), log(nopen)))

# Too difficult to pick out any patterns when overlaying individuals in each yr
ggplot(df, aes(x = doy, y = nopen_log0)) +
  geom_line(aes(color = factor(id))) +
  facet_wrap(~ yr) +
  theme_bw() +
  theme(legend.position = "none")
# Though you can see that 2023-2024 had many individuals with very few flowers

# Trajectories for one year at a time
ggplot(filter(df, yr == 2024), aes(x = doy, y = nopen_log0)) +
  geom_line() +
  facet_wrap(~ id) +
  theme_bw()
# 2018: pretty good. Many individuals have clear up then down.
# 2019: similar. If we wanted to pick out a "peak" date, then may want to
  # think about filters for individuals that had just a single non-0 value,
  # particularly if that value is first or last observation date
# 2020: pretty good
# 2021: pretty good
# 2022: many fewer individuals, but trajectories ok
# 2023: trajectories ok, but many plants with just a single open flower
# 2024: trajectories ok, but many plants with just a single open flower

# Look at 2024 plants
filter(dfmax, yr == 2024) %>%
  left_join(distinct(sibw, id, site_id), by = "id") %>%
  left_join(sites, by = "site_id") %>%
  arrange(site_id, maxopen)
# 16 of 31 at CHIR/FOBO and all but one of them have a max open flower count = 1

# Look at 2023 plants
filter(dfmax, yr == 2023) %>%
  left_join(distinct(sibw, id, site_id), by = "id") %>%
  left_join(sites, by = "site_id") %>%
  arrange(site_id, maxopen)
# 25 of 33 plants at CHIR/FOBO/CORO and all but 2 of them have a max open flower count = 1

# Look at 2019 plants
filter(dfmax, yr == 2019) %>%
  left_join(distinct(sibw, id, site_id), by = "id") %>%
  left_join(sites, by = "site_id") %>%
  arrange(site_id, maxopen)
# 11 of 38 plants at CHIR/FOBO/CORO but max counts are varied (good!)
# 12 of 38 plants at Karchner Caverns and they all have max counts of 1

# How are we ending up with so many max open flower counts = 1?
  # What were the original intensity categories again?
  count(filter(sib, !is.na(intensity_midpoint)), 
        php, intensity_midpoint, intensity_value)
  
  # What %open values were provided when the flower count was reported to be 
  # "Less than 3" (midpoint = 1). There can only be 1 or 2 flowers, so the 
  # only possible correct options if any flowers are open are:
    # 50-74% (midpoint = 62%) if 1/2 flowers are open
    # 95% or more (midpoint = 95%) if 1/1 or 2/2 are open
  count(filter(df, intensity_flowers == 1 & df$intensity_open > 0), 
        intensity_flowers, intensity_open, nopen) %>%
    mutate(prop = round(n/sum(df$intensity_flowers == 1 & df$intensity_open > 0), 2))
  # Most common intensity values are: 
    # 37% (25-49%)
    # 14% (5-24%)
    # other impossible values like <5%, 75-94% also reported
  
######## 
# Just realized that a lot of observers are reporting the number of flowering 
# stalks for the flowers intensity value (as they're instructed to do in the 
# Flowers for Bats campaign and elsewhere). Notably, observers are asked to 
# estimate the proportions of flowers on each stalk that are open. So the two
# intensity values are taken are different scales (no. flowers = stalk;
# %open = individual flowers).
########

# CAN'T USE THIS DATASET TO UNDERSTAND VARIATION IN THE NUMBER OF OPEN FLOWERS
