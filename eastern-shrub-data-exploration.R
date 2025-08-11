# Exploring intensity data for shrubs in (south?)eastern US
# ER Zylstra

library(tidyverse)
library(leaflet)
library(ordbetareg)   # Loads brms 
library(tidybayes)    # Manipulate Stan objects in a tidy way
# library(broom)        # Convert model objects to data frames
# library(broom.mixed)  # Convert brms model objects to data frames
# library(emmeans)      # Calculate marginal effects in even fancier ways

# Load and format status-intensity data ---------------------------------------#

# List files with formatted intensity data
intensity_files <- list.files("npn-data",
                              pattern = "intensity-shrub",
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

# Calculate the interval between observations of the plant, phenophase
si <- si %>%
  arrange(common_name, individual_id, phenophase_description, yr, day_of_year)

si$interval_raw <- c(NA, si$day_of_year[2:nrow(si)] - si$day_of_year[1:(nrow(si) - 1)])
same_ind <- 1 * (si$individual_id[2:nrow(si)] == si$individual_id[1:(nrow(si) - 1)])
same_php <- 1 * (si$phenophase_id[2:nrow(si)] == si$phenophase_id[1:(nrow(si) - 1)])
same_yr <- 1 * (si$yr[2:nrow(si)] == si$yr[1:(nrow(si) - 1)])
si$same_ind <- c(NA, same_ind)
si$same_php <- c(NA, same_php)
si$same_yr <- c(NA, same_yr)
si <- si %>%
  mutate(interval = case_when(
    same_ind == 0 ~ NA,
    same_php == 0 ~ NA,
    same_yr == 0 ~ NA,
    is.na(interval_raw) ~ NA,
    .default = interval_raw
  )) %>%
  dplyr::select(-c(interval_raw, same_ind, same_php, same_yr))

# si %>% count(interval) %>% mutate(prop = n / 63829) %>% round(3)

# Check that there's only one intensity category for each species-phenophase?
spil <- si %>%
  filter(!is.na(intensity_category_id)) %>%
  distinct(common_name, phenophase_description, intensity_name,
           intensity_type, intensity_label)
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
# is yes but no intensity value was provided

# Filter data: plant-phenophase-year combinations -----------------------------#

si %>%
  group_by(state) %>%
  summarize(nplants = n_distinct(individual_id),
            beautyberry = n_distinct(individual_id[common_name == "American beautyberry"]),
            buttonbush = n_distinct(individual_id[common_name == "common buttonbush"])) %>%
  data.frame() %>%
  arrange(desc(nplants))

sites <- si %>%
  group_by(site_id, site_name, latitude, longitude) %>%
  summarize(n_plants = n_distinct(individual_id), 
            ambe = ifelse("American beautyberry" %in% common_name, 1, 0),
            cobu = ifelse("common buttonbush" %in% common_name, 1, 0),
            .groups = "keep") %>%
  data.frame()

# Map
leaflet(sites) %>% addTiles() %>%
  addCircleMarkers(lng = ~longitude, lat = ~latitude, radius = ~n_plants, 
                   fillOpacity = 0.6)

leaflet(sites) %>% addTiles() %>%
  addCircles(lng = ~longitude, lat = ~latitude, 
             data = filter(sites, ambe == 1), 
             group = "Beautyberry",
             color = "red", fillOpacity = 1, weight = 10) %>%
  addCircles(lng = ~longitude, lat = ~latitude, 
             data = filter(sites, cobu == 1), 
             group = "Buttonbush",
             color = "blue", fillOpacity = 1, weight = 10) %>%
  addLayersControl(overlayGroups = c("Beautyberry", "Buttonbush"),
                   options = layersControlOptions(collapse = FALSE))

# For now, just remove observations in the west (CA, OR), but may want to later 
# constrain beautyberry locations to the SE (FL, LA, TX, NC, MS, AR, AL, GA, OK, 
# SC, TN, maybe VA?)
si <- si %>%
  filter(longitude > -100) %>%
  mutate(se = ifelse(state %in% c("FL", "LA", "TX", "NC", "MS", "AR", "AL",
                                  "GA", "OK", "SC", "TN", "VA"), 1, 0))

# Just going to focus on open flowers for now, so eliminating everything else
si <- filter(si, phenophase_description == "Open flowers")

# First glance at dates plants were in open flower phenophase
si %>%
  filter(phenophase_status == 1) %>%
  ggplot(aes(x = day_of_year)) +
  geom_histogram() +
  facet_grid(common_name ~.)

# Want to:
# remove plant-year combos with no observations in phase
# remove plant-year combos with no observations with intensity values
# remove plant-year combos with < 5 observations
# remove plant-year combos when max interval > 21 days

# Summarize amount and quality of information for each plant, year
pl_yr <- si %>%
  group_by(site_id, common_name, individual_id, yr) %>%
  summarize(nobs = n(),
            first_obs = min(day_of_year),
            last_obs = max(day_of_year),
            mean_int = round(mean(interval, na.rm = TRUE), 2),
            max_int = ifelse(nobs ==1, NA, max(interval, na.rm = TRUE)),
            n_inphase = sum(phenophase_status),
            n_intvalue = sum(!is.na(intensity_value)),
            prop_inphase = round(n_inphase / nobs, 2),
            prop_intvalue = round(n_intvalue / n_inphase, 2),
            .groups = "keep") %>%
  data.frame()

pl_yr <- pl_yr %>%
  mutate(remove = case_when(
    n_inphase == 0 ~ 1,
    n_intvalue == 0 ~ 1,
    nobs < 5 ~ 1,
    # max_int > 21 ~ 1,
    max_int > 30 ~ 1,
    .default = 0
  ))
sif <- si %>%
  left_join(select(pl_yr, individual_id, yr, remove),
            by = c("individual_id", "yr")) %>%
  filter(remove == 0) %>%
  select(-remove) %>%
  mutate(indyr = paste0(individual_id, "_", yr))

sif %>%
  group_by(common_name) %>%
  summarize(nobs = n(),
            nseries = n_distinct(indyr),
            nplants = n_distinct(individual_id),
            nsites = n_distinct(site_id),
            nyrs = n_distinct(yr)) %>%
  data.frame()

sif %>%
  filter(common_name == "American beautyberry") %>%
  ggplot(aes(x = day_of_year, y = intensity_midpoint)) +
  geom_line() +
  facet_wrap(~ indyr)
sif %>%
  filter(common_name == "common buttonbush") %>%
  ggplot(aes(x = day_of_year, y = intensity_midpoint)) +
  geom_line() +
  facet_wrap(~ indyr)
