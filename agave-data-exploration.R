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

si %>%
  filter(remove2 == 0) %>%
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

# Run model for one year ------------------------------------------------------#

# May want to fill in missing elevation values

# Use one of the filters, create factor variables, and standardize DOY
sif2_2023 <- sif %>%
  filter(remove2 == 0) %>%
  filter(yr == 2023) %>%
  mutate(doyz = (day_of_year - mean(day_of_year)) / sd(day_of_year)) %>%
  mutate(id = factor(individual_id),
         # fyr = factor(yr),
         site = factor(site_id))

m_2023 <- ordbetareg(prop ~ doyz + (1|id),
                     data = sif2_2023,
                     control = list(adapt_delta = 0.9),
                     cores = 4, chains = 4)
summary(m_2023)

# For a good explanation of the different predictions types (grand mean,
# marginal effects), see:
# https://www.andrewheiss.com/blog/2021/11/10/ame-bayes-re-guide/

# Expected proportion open by date, ignoring individual effects (ie, grand mean)
doy_gm <- m_2023 %>%
  epred_draws(newdata = data.frame(doyz = seq(min(sif2_2023$doyz),
                                              max(sif2_2023$doyz),
                                              length = 100)),
              re_formula = NA)
# This creates a large tibble, with nrows = unique(doy) * no. post samples
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

# Run models for multiple years -----------------------------------------------#

# Use one of the filters, create factor variables, and standardize DOY
sif2 <- sif %>%
  filter(remove2 == 0) %>%
  mutate(doyz = (day_of_year - mean(day_of_year)) / sd(day_of_year)) %>%
  mutate(id = factor(individual_id),
         fyr = factor(yr),
         site = factor(site_id))

# Modeling year as a random effect
m_allyrs <- ordbetareg(prop ~ doyz + (1 + doyz|fyr) + (1|id),
                       data = sif2,
                       # control = list(adapt_delta = 0.9),
                       cores = 4, chains = 4)
summary(m_allyrs)
# Less than 10 min to run, but problems with divergent transitions and ESS
coef(m_allyrs)$fyr

# Expected proportion open by year, ignoring individual REs
doy_yr <- m_allyrs %>%
  epred_draws(newdata = expand_grid(doyz = seq(min(sif2$doyz),
                                               max(sif2$doyz),
                                               length = 100),
                                    fyr = levels(sif2$fyr)),
              re_formula = ~ (1 + doyz | fyr))
plot_doy_yr <- ggplot(doy_yr, 
                      aes(x = doyz, y = .epred)) +
  stat_lineribbon(aes(color = fyr, fill = after_scale(alpha(color, 0.2))), 
                  .width = 0.80) +
  labs(x = "Day of year", y = "Predicted proportion of flowers open (%)", 
       color = "Year",
       fill = "Year", 
       title = "By year, correlated intercepts/slopes") +
  theme_bw() +
  theme(legend.position = "bottom")
plot_doy_yr

plot_doy_yr2 <- ggplot(filter(doy_yr, doyz > (-1)), 
                       aes(x = doyz, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  facet_wrap(~ fyr, ncol = 2) +
  labs(x = "Day of year", y = "Predicted proportion of flowers open (%)", 
       fill = "Credible interval", 
       title = "By year, correlated intercepts/slopes") +
  theme_bw() +
  theme(legend.position = "bottom")
plot_doy_yr2


