# Exploring intensity data for eastern redbud
# ER Zylstra

library(tidyverse)
library(leaflet)
library(elevatr)
library(sf)
# library(ordbetareg)   # Loads brms 
# library(tidybayes)    # Manipulate Stan objects in a tidy way

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
  group_by(site_id, site_name, state, latitude, longitude) %>%
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
  filter(longitude > (-100)) %>%
  filter(site_id != qcsite)
si <- si %>%
  filter(site_id %in% unique(sites$site_id))

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

# How much data each year?
flowersw %>%
  group_by(yr) %>%
  summarize(n_plants = n_distinct(id),
            n_sites = n_distinct(site_id),
            n_obs = n(),
            n_int_flowers = sum(!is.na(intensity_flowers)),
            n_int_open = sum(!is.na(intensity_open)),
            n_int_both = sum(!is.na(intensity_flowers) & !is.na(intensity_open))) %>%
  data.frame()
# Way more data in 2022-2025 than previous years

# Plot mean daily intensity values for 2022 (was going to use median, but 
# there are too many 0s)
flowersw %>%
  filter(yr == 2022) %>%
  group_by(doy) %>%
  summarize(flowers = mean(intensity_flowers, na.rm = TRUE),
            open = mean(intensity_open, na.rm = TRUE)) %>%
  mutate(log_flowers = ifelse(flowers == 0, log(flowers + 0.1), log(flowers))) %>%
  select(-flowers) %>%
  pivot_longer(cols = c(log_flowers, open),
               names_to = "type",
               values_to = "intensity") %>%
  filter(doy %in% 25:200) %>%
  ggplot(aes(x = doy, y = intensity)) +
  geom_line() +
  facet_grid(type ~ ., scales = "free_y")
# Messy. However, even here it looks like we can have moderate %open values
# late in the flowering season when the number of flowers decreases 
# Try weekly?
flowersw %>%
  filter(yr == 2022) %>%
  mutate(wk = as.numeric(week(observation_date))) %>%
  group_by(wk) %>%
  summarize(flowers = mean(intensity_flowers, na.rm = TRUE),
            open = mean(intensity_open, na.rm = TRUE)) %>%
  mutate(log_flowers = ifelse(flowers == 0, log(flowers + 0.1), log(flowers))) %>%
  select(-flowers) %>%
  pivot_longer(cols = c(log_flowers, open),
               names_to = "type",
               values_to = "intensity") %>%
  filter(wk %in% 9:25) %>%
  ggplot(aes(x = wk, y = intensity)) +
  geom_line() +
  facet_grid(type ~ ., scales = "free_y")
# Looking at weekly values averaged across individuals, this seems like less of 
# a problem...



# What about looking at data for an individual?
count(filter(flowersw, yr == 2022), id) %>% arrange(desc(n)) %>% head(10)
flowersw %>%
  filter(yr == 2022 & id == 287224) %>%
  mutate(flowers_log = ifelse(intensity_flowers == 0, 
                              log(intensity_flowers + 0.1), 
                              log(intensity_flowers))) %>%
  select(-intensity_flowers) %>%
  pivot_longer(cols = c(flowers_log, intensity_open),
               names_to = "type",
               values_to = "intensity") %>%
  # filter(doy %in% 25:200) %>%
  ggplot(aes(x = doy, y = intensity)) +
  geom_line() +
  facet_grid(type ~ ., scales = "free_y")
# %open stays high as number of flowers drops off
