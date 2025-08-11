# Exploring intensity data for Palmer's agave
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
count(si, phenophase_description, phenophase_status, intensity_midpoint)

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
                   fillOpacity = 0.6)

# Will remove plants from CA, TX, and the Phoenix area
sites <- sites %>%
  filter(latitude < 33) %>%
  filter(longitude > (-114) & longitude < (-109))
si <- si %>%
  filter(site_id %in% unique(sites$site_id))

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
# remove plant-year combos with < 5 observations between days 100 and 200
# remove plant-year combos when max interval between days 100 and 200 is >21 days
#### May want to revisit this interval length filter

# Summarize amount and quality of information for each plant, year
pl_yr <- si %>%
  arrange(individual_id, observation_date) %>%
  group_by(site_id, individual_id, yr) %>%
  summarize(nobs = n(),
            nobs_100_300 = sum(day_of_year %in% 100:300),
            nobs_120_210 = sum(day_of_year %in% 120:210),
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

pl_yr2 <- si %>%
  arrange(individual_id, observation_date) %>%
  filter(day_of_year %in% 120:210) %>%
  group_by(site_id, individual_id, yr) %>%
  summarize(nobs = n(),
            max_int_pre200 = ifelse(nobs == 1 & anyNA(interval), NA, 
                                    max(interval, na.rm = TRUE)),
            .groups = "keep") %>%
  data.frame()

pl_yr <- pl_yr %>%
  left_join(select(pl_yr2, -nobs), by = c("site_id", "individual_id", "yr"))

pl_yr <- pl_yr %>%
  mutate(remove = case_when(
    n_inphase == 0 ~ 1,
    n_intvalue == 0 ~ 1,
    nobs_100_300 < 5 ~ 1,
    max_int_pre200 > 21 ~ 1,
    .default = 0
  ))
sif <- si %>%
  left_join(select(pl_yr, individual_id, yr, remove),
            by = c("individual_id", "yr")) %>%
  filter(remove == 0) %>%
  select(-remove)


# What's left?
sif %>%
  group_by(yr) %>%
  summarize(n_obs = n(),
            n_sites = n_distinct(site_id),
            n_plants = n_distinct(individual_id),
            obs_per_plant = n_obs / n_plants) %>%
  data.frame()

# Remove data from 2017 (just one plant; all other years 17-55 plants)
sif <- filter(sif, yr > 2017)

# Extract information about sites ---------------------------------------------#

sites <- sif %>%
  distinct(latitude, longitude, site_id)
# write.table(sites, "weather-data/agave-sites.csv", sep = ",",
#             row.names = FALSE,
#             col.names = FALSE)

# Plot raw intensity data -----------------------------------------------------#

# Create filename for png
png_name <- "output/IntensityData-OpenFlowers-agave.png" 

# Create ggplot object and save to file (if file doesn't already exist)
if (!file.exists(png_name)) {

  # Figures with lines and points, with size proportional to no. of plants 
  agg <- sif %>%
    group_by(yr, day_of_year, intensity_midpoint) %>%
    summarize(n_indiv = n(), .groups = "keep") %>%
    data.frame()
  
  iplot <- ggplot(agg, aes(x = day_of_year, y = intensity_midpoint)) +
    facet_grid(yr ~ .) +
    geom_line(data = sif, show.legend = FALSE, alpha = 0.5,
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

sif %>%
  filter(yr == 2023) %>%
  ggplot(aes(x = day_of_year, y = intensity_midpoint)) +
  geom_line() +
  facet_wrap(~ individual_id) +
  theme_bw()

# These plots are pretty interesting -- gradual increase in % open for 
# individual plants.

# Not sure how I want to deal with the data after peak intensity value, since I 
# really just want to model the increase.....
  # Make all intensity values = NA when status goes back down to zero?



