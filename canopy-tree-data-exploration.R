# Exploring intensity data for canopy trees in southern Appalachian Trail area
# ER Zylstra

library(tidyverse)
library(leaflet)
library(ordbetareg)   # Loads brms 
library(tidybayes)    # Manipulate Stan objects in a tidy way
library(broom)        # Convert model objects to data frames
library(broom.mixed)  # Convert brms model objects to data frames
library(emmeans)      # Calculate marginal effects in even fancier ways

# Load and format status-intensity data ---------------------------------------#

# List files with formatted intensity data
intensity_files <- list.files("npn-data",
                              pattern = "intensity-canopy",
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

# si %>% count(interval) %>% mutate(prop = n / 280000) %>% round(3)

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

# Most observations in TN, NC, and VA, with others in WV, SC. One site/plant 
# listed as ND, which I'll remove just in case. Keeping NAs, since lat/lon seem
# ok and it's much easier to understand how that happens
si <- filter(si, is.na(state) | state != "ND")

# Only a few red mapl etrees east of -78 longitude, so we'll filter them out
si <- filter(si, longitude <= (-78))

# Just going to focus on canopy fullness for now, so eliminating everything else
si <- filter(si, phenophase_description == "Leaves")

# Going to focus on canopy closure in spring, so limiting observations to 
# the first half of the year (day 182 is June 30/July 1)
si <- filter(si, day_of_year <= 182)

# Want to:
  # remove plant-year combos with no observations in phase
  # remove plant-year combos with no observations with intensity values
  # remove plant-year combos with < 5 observations
  # remove plant-year combos when max interval > 21 days
  #### May want to revisit this interval length filter

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
    max_int > 21 ~ 1,
    # max_int > 30 ~ 1,
    .default = 0
  ))
si <- si %>%
  left_join(select(pl_yr, individual_id, yr, remove),
            by = c("individual_id", "yr")) %>%
  filter(remove == 0) %>%
  select(-remove)

# Looks like there might be instances where an individual tree had high canopy
# values with a 0 value in-between. This isn't really possible, at least in such
# a short amount of time, so we'll try to filter these data out since they're
# likely due to observation error.

# Look for patterns with consecutive intensity values >80, then <50, then
# back up
sif <- si %>%
  mutate(intensity_cat = ifelse(intensity_midpoint > 80, 2,
                                ifelse(intensity_midpoint < 50, 0, 1))) %>%
  # Remove observations in phase with no intensity value
  filter(!is.na(intensity_midpoint)) %>%
  arrange(individual_id, observation_date)

# Where do high-low-high patterns occur?
pattern <- c(2, 0, 2)
starts <- which(zoo::rollapply(sif$intensity_cat, 
                               length(pattern), 
                               function(x) all(x == pattern)))

# Identify anomalously low observations so we can filter them out
sif$anomalous <- 0
sif$anomalous[starts + 1] <- 1

# Filter them out
sif <- sif %>%
  filter(anomalous == 0) %>%
  select(-anomalous)

# What's left?
spp_yr <- sif %>% 
  group_by(common_name, yr) %>%
  summarize(n_trees = n_distinct(individual_id),
            .groups = "keep") %>%
  data.frame()
spp_yr
count(spp_yr, n_trees)
sum(spp_yr$n_trees >= 5) / nrow(spp_yr)
# 48/53 species-years (90%) have >= 5 trees. 
# American basswood has < 10 trees every year, only 5 years with >= 5 trees

# Extract information about sites ---------------------------------------------#

sites <- sif %>%
  distinct(latitude, longitude, site_id)
# write.table(sites, "weather-data/canopy-plant-sites.csv", sep = ",",
#             row.names = FALSE,
#             col.names = FALSE)

sites2 <- sif %>%
  group_by(site_id, site_name, latitude, longitude) %>%
  summarize(n_spp = n_distinct(common_name),
            n_plants = n_distinct(individual_id),
            basswood = ifelse("American basswood" %in% common_name, 1, 0),
            oak = ifelse("northern red oak" %in% common_name, 1, 0),
            beech = ifelse("American beech" %in% common_name, 1, 0),
            sugarm = ifelse("sugar maple" %in% common_name, 1, 0),
            redm = ifelse("red maple" %in% common_name, 1, 0),
            stripedm = ifelse("striped maple" %in% common_name, 1, 0),
            .groups = "keep") %>%
  data.frame()

# Map
leaflet(sites2) %>% addTiles() %>%
  addCircleMarkers(lng = ~longitude, lat = ~latitude, radius = ~n_spp, 
                   fillOpacity = 0.6)
# Map by species
leaflet(sites2) %>% addTiles() %>%
  addCircles(lng = ~longitude, lat = ~latitude, 
                   data = filter(sites2, redm == 1), 
                   group = "Red maple",
                   color = "red", fillOpacity = 1, weight = 10) %>%
  addCircles(lng = ~longitude, lat = ~latitude, 
                   data = filter(sites2, sugarm == 1), 
                   group = "Sugar maple",
                   color = "black", fillOpacity = 1, weight = 10) %>%
  addCircles(lng = ~longitude, lat = ~latitude, 
                   data = filter(sites2, oak == 1), 
                   group = "Northern red oak",
                   color = "purple", fillOpacity = 1, weight = 10) %>%
  addCircles(lng = ~longitude, lat = ~latitude, 
                   data = filter(sites2, beech == 1), 
                   group = "American beech",
                   color = "blue", fillOpacity = 1, weight = 10) %>%
  addCircles(lng = ~longitude, lat = ~latitude, 
                   data = filter(sites2, basswood == 1), 
                   group = "American basswood",
                   color = "gray", fillOpacity = 1, weight = 10) %>%
  addCircles(lng = ~longitude, lat = ~latitude, 
                   data = filter(sites2, stripedm == 1), 
                   group = "Striped maple",
                   color = "orange", fillOpacity = 1, weight = 10) %>%
  addLayersControl(overlayGroups = c("American basswood",
                                     "American beech",
                                     "Northern red oak",
                                     "Red maple",
                                     "Striped maple",
                                     "Sugar maple"),
                   options = layersControlOptions(collapse = FALSE))

sites2 %>% arrange(desc(n_plants)) %>% head(20)
# Two NEON sites have way more plants (oak, red maple, striped maple) than 
# any other sites (Great Smoky Mtns NP; Mountain Lake Biol Station):
# GRSM_068.phenology.phe - primary
# MLBS_077.phenology.phe - primary
sites2 %>% filter(grepl(".phenology.", site_name))
# Also one NEON phenocam site:
# GRSM_067.phenology.phe - phenocam

# Plot raw intensity data -----------------------------------------------------#

# Identify the min/max number of plants per species and year to keep point size 
# consistent among plots
nplants_size <- sif %>%
  group_by(intensity_label, common_name, yr) %>%
  summarize(n_plants = n_distinct(individual_id), .groups = "keep") %>%
  data.frame()

# Loop through species
spps <- unique(sif$common_name)
  
for (spp in spps) {
    
  # Create filename for png
  spp_nospace <- str_replace_all(spp, " ", "_")
  png_name <- paste0("output/canopy-intensities-by-spp-yr/IntensityData-CanopyFullness-", 
                     spp_nospace, ".png")
    
  # Create ggplot object and save to file (if file doesn't already exist)
  if (!file.exists(png_name)) {
    
    si_int_spp <- filter(sif, common_name == spp)

    # Figures with lines and points, with size proportional to no. of plants 
    agg <- si_int_spp %>%
      group_by(yr, day_of_year, intensity_midpoint) %>%
      summarize(n_indiv = n(), .groups = "keep") %>%
      data.frame()
    
    iplot <- ggplot(agg, aes(x = day_of_year, y = intensity_midpoint)) +
      facet_grid(yr ~ .) +
      geom_line(data = si_int_spp, show.legend = FALSE, alpha = 0.5,
                aes(color = factor(individual_id))) +
      geom_point(data = agg, aes(size = n_indiv), shape = 16, alpha = 0.4) +
      # scale_size_continuous(limits = c(1, max(nplants_size$n_plants))) +
      labs(title = paste0(str_to_sentence(spp), ",  Canopy fullness"),
           y = "Leaf canopy fullness (%)", x = "Day of year") +
      theme_bw()
    iplot
    
    ggsave(png_name,
           plot = iplot,
           width = 6.5,
           height = 9,
           units = "in")
  }
}

# American beech is an outlier, since looks like some of trees retained their
# leaves through the winter, especially in later years.

# Look at:
# Weird American beech in 2017
# Weird red oak in 2016, 2017, 2019
# Occasional other odd trajectories (maybe make a rule about excluding any 
# plants that have decreasing or drop offs in intensity values? They could be 
# real if the plant is sick or if there's some weather event, but it could
# also be observer error)

# Look at some anomalies:
sif %>%
  filter(common_name == "northern red oak",
         yr == 2016) %>%
  ggplot(aes(x = day_of_year, y = intensity_midpoint)) +
    geom_line() +
    facet_wrap(~ individual_id) +
    labs(x = "Day of year", y = "Canopy fullness (%)",
         title = "Red oak, 2016")
# Several declining trajectories (29929, 45483, 49999, 50001, 103295)

sif %>%
  filter(common_name == "northern red oak",
         yr == 2017) %>%
  ggplot(aes(x = day_of_year, y = intensity_midpoint)) +
  geom_line() +
  facet_wrap(~ individual_id) +
  labs(x = "Day of year", y = "Canopy fullness (%)",
       title = "Red oak, 2017")
# Bunch that bounce between very low and very high (23062, 23067, 29929, 45483)

sif %>%
  filter(common_name == "northern red oak",
         yr == 2019) %>%
  ggplot(aes(x = day_of_year, y = intensity_midpoint)) +
  geom_line() +
  facet_wrap(~ individual_id) +
  labs(x = "Day of year", y = "Canopy fullness (%)",
       title = "Red oak, 2019")
# Couple that bounce between very low and very high (23062, 166287)
# One that only has observations over short period with highest value

# Create model for one year and make (inverse) predictions --------------------#

# Red maple in 2017
rema17 <- sif %>%
  filter(common_name == "red maple",
         yr == 2017)
# Convert 0-95 percentages to 0-1 proportions
rema17 <- rema17 %>%
  mutate(prop = intensity_midpoint / 95)
# For prediction later, want to know what proportion = 50% intensity value
prop50 <- 50/95

m_rema17 <- ordbetareg(prop ~ day_of_year + (1|individual_id),
                      data = rema17, 
                      # true_bounds = c(0, 95),
                      # control = list(adapt_delta = 0.9),
                      cores = 4, chains = 4,
                      backend = "cmdstanr")
summary(m_rema17)
plot(m_rema17) # posterior distributions and trace plots

# For a good explanation of the different predictions types (grand mean,
# marginal effects), see:
# https://www.andrewheiss.com/blog/2021/11/10/ame-bayes-re-guide/

# Expected canopy fullness by date, ignoring individual effects (ie, grand mean)
doy_gm <- m_rema17 %>%
  epred_draws(newdata = data.frame(day_of_year = seq(0, 180, by = 5)),
              re_formula = NA)
# This creates a large tibble, with nrows = unique(doy) * no. post samples
plot_doy_gm <- ggplot(filter(doy_gm, day_of_year %in% 25:180), 
                      aes(x = day_of_year, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Day of year", y = "Predicted canopy fullness (%)",
       fill = "Credible interval", title = "Ignoring individual effects") +
  theme_bw() +
  theme(legend.position = "bottom")
plot_doy_gm

# Conditional effects for a new, typical individual (new individual is based on
# variation among existing individuals)
doy_typicalplant <- m_rema17 %>%
  epred_draws(newdata = data.frame(day_of_year = seq(0, 180, by = 5),
                                   individual_id = "new plant"),
              re_formula = NULL,
              allow_new_levels = TRUE,
              sample_new_levels = "uncertainty")
plot_doy_typicalplant <- ggplot(filter(doy_typicalplant, day_of_year %in% 25:180), 
                                aes(x = day_of_year, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Day of year", y = "Predicted canopy fullness (%)",
       fill = "Credible interval", title = "New, typical plant") +
  theme_bw() +
  theme(legend.position = "bottom")
plot_doy_typicalplant

# Conditional effects for a brand new individual (new individual is based on
# random draws from model)
doy_newplant <- m_rema17 %>%
  epred_draws(newdata = data.frame(day_of_year = seq(0, 180, by = 5),
                                   individual_id = "new plant"),
              re_formula = NULL,
              allow_new_levels = TRUE,
              sample_new_levels = "gaussian")
plot_doy_newplant <- ggplot(filter(doy_newplant, day_of_year %in% 25:180), 
                                aes(x = day_of_year, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Day of year", y = "Predicted canopy fullness (%)",
       fill = "Credible interval", title = "Brand new plant") +
  theme_bw() +
  theme(legend.position = "bottom")
plot_doy_newplant


# Inverse prediction: what day of year will predicted canopy fullness be 50%?
# See: https://discourse.mc-stan.org/t/inverse-predictions-from-posterior/35350

# Without taking random effects into account:
d <- as_draws_rvars(m_rema17)
# d 
# Note that this will hang for a while to summarize random intercepts
# Don't need to wait though - can just grab parameter names

intercept_p <- d$`b_Intercept`
slope_p <- d$`b_day_of_year`

x_50_p = (qlogis(prop50) - intercept_p) / slope_p
str(posterior::draws_of(x_50_p))
posterior::summarize_draws(x_50_p)
quantile(posterior::draws_of(x_50_p)[, 1], probs = c(0.025, 0.5, 0.975))

# But what if we do want to take random effects into account?
# Need to add gaussian noise by sampling from normal(0, sd(indidividual_id))
samples <- coef(m_rema17, summary = FALSE) 
# This is a list of 1 array (name = individual_id). 
# Dims of array: rows = iter (4000); cols = individual (83); slice = param (2 [int, doy])

oint <- rowMeans(samples$individual_id[,,1])
# Summarizing the mean intercept across individuals for each iteration
# gets you very close to summary of Intercept in summary(m_rema17)

# So to take variation across individuals into account, use matrix of intercept
# values instead of mean across individuals
intercepts <- as.vector(samples$individual_id[,,1])
slopes <- as.vector(samples$individual_id[,,2])
x_50_pRE = (qlogis(prop50) - intercepts) / slopes
summary(x_50_pRE)
quantile(x_50_pRE, probs = c(0.025, 0.50, 0.975))

# Look at variation among plants in/out of NEON sites -------------------------#
# Still working with red maple 2017 data

# First, identify NEON sites:
rema17 <- rema17 %>%
  mutate(neon = ifelse(grepl(".phenology.", site_name), 1, 0))
rema17 %>%
  distinct(individual_id, site_name, neon) %>%
  count(neon, site_name)
# NEON site has 30 plants
# GRSM-Tremont-Lumber Ridge has 14 plants
# All other sites have <= 3 plants

rema17 <- rema17 %>%
  mutate(group3 = case_when(
    grepl(".phenology.", site_name) ~ "neon",
    grepl("Lumber Ridge", site_name) ~ "lumber",
    .default = "other"
  ))

# Summarize how much of the data come from NEON sites
rema17 %>%
  group_by(group3) %>%
  summarize(n_sites = n_distinct(site_name),
            n_plants = n_distinct(individual_id),
            n_obs = n()) %>%
  data.frame()

# Where is the NEON site relative to the others?
rema17sites <- rema17 %>%
  distinct(site_name, latitude, longitude, group3)
leaflet(filter(rema17sites, group3 == "other")) %>% 
  addTiles() %>%
  addCircleMarkers(lng = ~longitude, 
                   lat = ~latitude, 
                   fillColor = "yellow", 
                   fillOpacity = 1,
                   radius = 5,
                   color = NA) %>%
  addCircleMarkers(data = filter(rema17sites, group3 == "neon"),
                   lng = ~longitude, 
                   lat = ~latitude, 
                   fillColor = "red", 
                   fillOpacity = 1,
                   radius = 5,
                   color = NA) %>%
  addCircleMarkers(data = filter(rema17sites, group3 == "lumber"),
                   lng = ~longitude, 
                   lat = ~latitude, 
                   fillColor = "purple", 
                   fillOpacity = 1,
                   radius = 5,
                   color = NA)
# NEON and Lumber Ridge in GSMNP, with other sites. Remaining sites to the east

# Plot raw data
ggplot(data = rema17, aes(x = day_of_year, y = intensity_midpoint)) +
  geom_line(aes(color = factor(individual_id))) +
  facet_grid(group3 ~ .) +
  theme_bw() +
  theme(legend.position = "none")
# With just one year of data, it doesn't necessarily look like more sites 
# means more variation. Decent amount of variation in timing among plants 
# within NEON and Lumber Ridge sites. 
ggplot(data = filter(rema17, neon == 1),
       aes(x = day_of_year, y = intensity_midpoint)) +
  geom_line() +
  facet_wrap(~factor(individual_id))
# Is there something weird going on with the data around day 115? Drop in 
# canopy fullness for many trees. Did some event happen?
rema17 %>%
  group_by(observation_date) %>%
  summarize(nobs = n(),
            fullness = round(mean(intensity_midpoint))) %>%
  data.frame()
# April 24-25

# Multiple-year model ---------------------------------------------------------#
# Starting with red maples
rema <- sif %>%
  filter(common_name == "red maple")

# How are the data distributed across years, sites, individuals?
rema %>%
  group_by(yr) %>%
  summarize(nsites = n_distinct(site_name),
            nplants = n_distinct(individual_id),
            nobs = n()) %>%
  data.frame()
# Way less data in 2020

rema %>%
  group_by(site_name) %>%
  summarize(nyrs = n_distinct(yr),
            first_yr = min(yr),
            last_yr = max(yr),
            yr2020 = ifelse(2020 %in% yr, 1, 0),
            nplants = n_distinct(individual_id),
            nobs = n()) %>%
  arrange(desc(nplants), desc(nyrs), desc(nobs)) %>%
  data.frame()
# Two NEON sites have 30+ plants; Lumber Ridge has 14 plants
# All other sites have 10 or fewer plants

# Convert 0-95 percentages to 0-1 proportions
rema <- rema %>%
  mutate(prop = intensity_midpoint / 95)
# For prediction later, want to know what proportion = 50% intensity value
prop50 <- 50/95

# Create year factor (after removing 2020)
rema <- rema %>%
  filter(yr != 2020) %>%
  mutate(fyr = factor(yr))

# Model with year as fixed effect:
m_rema <- ordbetareg(prop ~ day_of_year * fyr + (1|individual_id),
                     data = rema, 
                     # control = list(adapt_delta = 0.9),
                     cores = 4, chains = 4,
                     backend = "cmdstanr")
# Took about an hour!
summary(m_rema)

# Expected canopy fullness by date, ignoring individual effects (ie, grand mean)
doy_yr_gm <- m_rema %>%
  epred_draws(newdata = expand.grid(day_of_year = seq(0, 180, by = 5),
                                    fyr = unique(rema$fyr),
                                    KEEP.OUT.ATTRS = FALSE),
              re_formula = NA)
# This creates a large tibble, with nrows = unique(doy) * unique(yr) * no. post samples
plot_doy_gm <- ggplot(filter(doy_yr_gm, day_of_year %in% 100:125, 
                             fyr %in% c("2016", "2017", "2018", "2019")), 
                      aes(x = day_of_year, y = .epred)) +
  stat_lineribbon(aes(color = fyr, fill = fyr), alpha = 1, .width = 0) +
  # facet_wrap(~ fyr) +
  labs(x = "Day of year", y = "Predicted canopy fullness (%)", color = "Year",
       fill = "Year",
       title = "Ignoring individual effects") +
  theme_bw() +
  theme(legend.position = "bottom")
plot_doy_gm
# It didn't look like things were matching up very well with the raw data plots.
# For example, this shows 2023 values increased earlier than in other years, 
# but looking closely at the raw plots, I'm wondering whether that's a function
# of a weird trajectory that was at 95% and help steady from very early in the
# year. So maybe the estimates are okay, but demostrate the need for more careful
# cleaning/filtering?

# Saving model object (actually entire workspace) for further exploration adn
# comparison:
# save.image("output/rema-multiyear-model.RData")

# Multiple species models -----------------------------------------------------#

# Find sites in GRSM
sites2 <- sites2 %>%
  mutate(grsm = ifelse(grepl("GRSM", site_name), 1, 0))
# Check that this designation is correct:
leaflet(sites2) %>% addTiles() %>%
  addCircles(lng = ~longitude, lat = ~latitude, 
             data = filter(sites2, grsm == 0), 
             group = "Other",
             color = "red", fillOpacity = 1, weight = 10) %>%
  addCircles(lng = ~longitude, lat = ~latitude, 
             data = filter(sites2, grsm == 1), 
             group = "GRSM",
             color = "black", fillOpacity = 1, weight = 10) %>%
  addLayersControl(overlayGroups = c("GRSM", "Other"),
                   options = layersControlOptions(collapse = FALSE))

grsm <- sif %>%
  filter(grepl("GRSM", site_name))
grsm %>% 
  group_by(yr) %>%
  summarize(nsites = n_distinct(site_name),
            nspp = n_distinct(common_name),
            nplants = n_distinct(individual_id),
            nobs = n()) %>%
  data.frame()
# 2020 has much less data (and only 4 species)
# 2021 seems like decent starting year?
grsm %>% 
  group_by(common_name) %>%
  summarize(nsites = n_distinct(site_name),
            nplants = n_distinct(individual_id),
            nyrs = n_distinct(yr),
            nobs = n()) %>%
  data.frame()
# Red maple most common, but most species represented decently well.

# Try one year to start:
grsm21 <- grsm %>%
  filter(yr == 2021) %>%
  mutate(spp = factor(common_name),
         id = factor(individual_id),
         prop = intensity_midpoint / 95) 
grsm21 %>%
  group_by(common_name) %>%
  summarize(nplants = n_distinct(id),
            nsites = n_distinct(site_name))
# Only 5 sugar maple and 6 basswood trees, so definitely should treat species
# as random effect

m_grsm21 <- ordbetareg(prop ~ day_of_year + (1|spp/id) ,
                       data = grsm21, 
                       cores = 4, chains = 4)
# Took 20 min to run. ESS for individual RE = 391 and warnings about 
# 8 divergent transitions, so probably should increase adapt_delta and maybe
# number of samples
summary(m_grsm21)

# Expected canopy fullness by date, ignoring random effects
doy_gm <- m_grsm21 %>%
  epred_draws(newdata = data.frame(day_of_year = seq(0, 180, by = 5)),
              re_formula = NA)
# This creates a large tibble, with nrows = unique(doy) * no. post samples
plot_doy_gm <- ggplot(filter(doy_gm, day_of_year %in% 25:180), 
                      aes(x = day_of_year, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Day of year", y = "Predicted canopy fullness (%)",
       fill = "Credible interval", title = "Ignoring individual effects") +
  theme_bw() +
  theme(legend.position = "bottom")
plot_doy_gm

# Expected canopy fullness for each species, ignoring individual REs
doy_spp <- m_grsm21 %>%
  epred_draws(newdata = expand_grid(day_of_year = seq(0, 180, by = 5),
                                    spp = levels(grsm21$spp)),
              re_formula = ~ (1 | spp))
plot_doy_spp <- ggplot(filter(doy_spp, day_of_year %in% 25:180), 
                              aes(x = day_of_year, y = .epred)) +
  # stat_lineribbon(aes(color = spp, fill = spp), alpha = 0.6, .width = 0.95) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  facet_wrap(~ spp) +
  labs(x = "Day of year", y = "Predicted canopy fullness (%)",
       fill = "Credible interval", title = "By species") +
  theme_bw() +
  theme(legend.position = "bottom")
plot_doy_spp
# Uncertainty varies as expected (lots of data, little uncertainty)
# But clear that we need to have random slopes for each species!

# Model with random slopes by species, random intercepts for species and tree
m_grsm21_rs <- ordbetareg(prop ~ day_of_year + (0 + day_of_year|spp) + (1|spp/id) ,
                          data = grsm21, 
                          iter = 1000,
                          cores = 4, chains = 4,
                          control = list(adapt_delta = 0.9))
# Took 35 min to run, but with higher delta and 1000 iterations. 
# ESS low for several parameters and warnings about 4 divergent transitions
summary(m_grsm21_rs)

# Species-level estimates from uncorrelated random slopes model
coef(m_grsm21_rs, summary = TRUE) 
# Is this specified correctly? It seems weird that day_of_year effects are 
# identical for each individual, but things look ok at the species level, so
# maybe it makes sense. Note that this model assumes no correlation between
# species-level intercepts and slopes, while the one below allows for correlations

# Trying again...
m_grsm21_rs2 <- ordbetareg(prop ~ day_of_year + (1 + day_of_year|spp) + (1|id),
                           data = grsm21, 
                           iter = 2000,
                           cores = 4, chains = 4,
                           control = list(adapt_delta = 0.9))
# Took 26 min to run. 
# ESS low for several parameters and warnings about 43 divergent transitions
summary(m_grsm21_rs2)
coef(m_grsm21_rs2)

# Look at predictions for both models.
# Expected canopy fullness for each species, ignoring individual REs
doy_rs_spp <- m_grsm21_rs %>%
  epred_draws(newdata = expand_grid(day_of_year = seq(0, 180, by = 5),
                                    spp = levels(grsm21$spp)),
              re_formula = ~ (0 + day_of_year | spp) + (1 | spp))
plot_doy_rs_spp <- ggplot(filter(doy_rs_spp, day_of_year %in% 25:180), 
                          aes(x = day_of_year, y = .epred)) +
  stat_lineribbon(aes(color = spp, fill = after_scale(alpha(color, 0.2))), 
                  .width = 0.80) +
  # stat_lineribbon() +
  # scale_fill_brewer(palette = "Blues") +
  # facet_wrap(~ spp) +
  labs(x = "Day of year", y = "Predicted canopy fullness (%)", color = "Species",
       fill = "Species", title = "By species") +
  theme_bw() +
  theme(legend.position = "bottom")
plot_doy_rs_spp

# Expected canopy fullness for each species, ignoring individual REs
doy_rs2_spp <- m_grsm21_rs2 %>%
  epred_draws(newdata = expand_grid(day_of_year = seq(0, 180, by = 5),
                                    spp = levels(grsm21$spp)),
              re_formula = ~ (1 + day_of_year | spp))
plot_doy_rs2_spp <- ggplot(filter(doy_rs2_spp, day_of_year %in% 25:180), 
                          aes(x = day_of_year, y = .epred)) +
  stat_lineribbon(aes(color = spp, fill = after_scale(alpha(color, 0.2))), 
                  .width = 0.80) +
  # stat_lineribbon() +
  # scale_fill_brewer(palette = "Blues") +
  # facet_wrap(~ spp) +
  labs(x = "Day of year", y = "Predicted canopy fullness (%)", color = "Species",
       fill = "Species", title = "By species, correlated intercepts/slopes") +
  theme_bw() +
  theme(legend.position = "bottom")
plot_doy_rs2_spp
# This actually looks quite reasonable, sugar maple is before the other species,
# while beech and northern red oak are the latest (and slowest)

# Tried a model with random slopes and intercepts for species and individuals
# but that took a lot longer and had a lot more issues with ESS and divergent
# transitions. For now at least, I think it's fine to assume slope doesn't vary 
# among individual trees.

# Re-running a model for 2021 to improve estimates
m_grsm21_rs3 <- ordbetareg(prop ~ day_of_year + (1 + day_of_year|spp) + (1|id),
                           data = grsm21, 
                           iter = 4000,
                           cores = 4, chains = 4,
                           control = list(adapt_delta = 0.99))

# Expected canopy fullness for each species, ignoring individual REs
doy_rs3_spp <- m_grsm21_rs3 %>%
  epred_draws(newdata = expand_grid(day_of_year = seq(0, 180, by = 5),
                                    spp = levels(grsm21$spp)),
              re_formula = ~ (1 + day_of_year | spp))
plot_doy_rs3_spp <- ggplot(filter(doy_rs3_spp, day_of_year %in% 25:180), 
                           aes(x = day_of_year, y = .epred)) +
  stat_lineribbon(aes(color = spp, fill = after_scale(alpha(color, 0.2))), 
                  .width = 0.80) +
  # stat_lineribbon() +
  # scale_fill_brewer(palette = "Blues") +
  # facet_wrap(~ spp) +
  labs(x = "Day of year", y = "Predicted canopy fullness (%)", color = "Species",
       fill = "Species", title = "By species, correlated intercepts/slopes") +
  theme_bw() +
  theme(legend.position = "bottom")
plot_doy_rs3_spp

# Multiple species and years? -------------------------------------------------#
# Not sure this is worth running given that we'd have to keep model complexity
# in check and thus, restrict how year might affect trajectories. 

grsm3yr <- grsm %>%
  filter(yr %in% 2021:2023) %>%
  mutate(fyr = factor(yr),
         spp = factor(common_name),
         id = factor(individual_id),
         prop = intensity_midpoint / 95) 

m_grsm3yr_rs <- ordbetareg(prop ~ day_of_year + fyr + (1 + day_of_year|spp) + (1|id),
                           data = grsm3yr, 
                           iter = 4000,
                           cores = 4, chains = 4,
                           control = list(adapt_delta = 0.99))

# save.image("output/grsm-multispp-models.RData")
summary(m_grsm3yr_rs)

# Expected canopy fullness for each species, ignoring individual REs
doy_3yr_spp <- m_grsm3yr_rs %>%
  epred_draws(newdata = expand_grid(day_of_year = seq(0, 180, by = 5),
                                    fyr = levels(grsm3yr$fyr),
                                    spp = levels(grsm3yr$spp)),
              re_formula = ~ (1 + day_of_year | spp))
plot_doy_3yr_spp <- ggplot(filter(doy_3yr_spp, day_of_year %in% 75:175), 
                           aes(x = day_of_year, y = .epred)) +
  stat_lineribbon(aes(color = spp, fill = after_scale(alpha(color, 0.2))), 
                  .width = 0.80) +
  facet_grid(fyr ~ .) +
  geom_hline(yintercept = 0.50, linetype = "dashed") +
  labs(x = "Day of year", y = "Predicted canopy fullness (%)", color = "Species",
       fill = "Species", title = "By species, correlated intercepts/slopes") +
  theme_bw() +
  theme(legend.position = "bottom")
plot_doy_3yr_spp
# Results ok -- clear that canopy filled slightly later in 2022 and much 
# earlier in 2023. However, I think the model structure may not be flexible 
# enough since differences among species in 2021 are not the same as they were 
# in the 2021 model. My guess is that they're constrained to be similar in all 
# years and so they may not be great for a single year. Could allow random 
# species-level effects for year, but a better choice might be using weather 
# variables instead of year.

# Formatting weather data -----------------------------------------------------#

# Stuff commented out below was only run once, to create monthly summaries
prism_folder <- "weather-data/canopy-plant-prism/"

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

# Calculate daily GDD values, with base temp of 0 deg C (gdd = tmean - base)
base <- 0
weather <- weather %>%
  mutate(gdd = ifelse(tmean < base, 0, tmean - base)) %>%
  mutate(yr = year(date),
         doy = yday(date)) %>%
  filter(doy <= 200) %>%
  arrange(site_id, date)

# Calculated AGDD values
weather <- weather %>%
  group_by(site_id, yr) %>%
  mutate(agdd = cumsum(gdd)) %>%
  ungroup() %>%
  data.frame()

# GDD models for red maple ----------------------------------------------------#

# Attach AGDD values to red maple data (2016-2019, 2021-2024; all sites)
agdds <- weather %>%
  select(site_id, date, agdd) %>%
  rename(observation_date = date)
  
rema <- rema %>%
  left_join(agdds, by = c("site_id", "observation_date"))

# Model with GDD instead of doy and year
m_rema_gdd <- ordbetareg(prop ~ agdd + (1|individual_id),
                         data = rema, 
                         control = list(adapt_delta = 0.99),
                         iter = 4000, cores = 4, chains = 4)
# Took over 7 hours!
summary(m_rema_gdd)
# save(m_rema_gdd, file = "output/rema-gdd-model.RData")

# Expected canopy fullness for each species, ignoring individual REs
gdd_rema <- m_rema_gdd %>%
  epred_draws(newdata = expand_grid(agdd = seq(0, 2800, by = 100)),
              re_formula = NA)
plot_gdd_rema <- ggplot(filter(gdd_rema, agdd %in% 0:1800), 
                        aes(x = agdd, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  geom_hline(yintercept = 0.50, linetype = "dashed") +
  labs(x = "AGDD (deg C)", y = "Predicted canopy fullness (%)",
       fill = "Credible interval",
       title = "Red maple, ignoring individual effects") +
  theme_bw() +
  theme(legend.position = "bottom")
plot_gdd_rema

# Conditional effects for a new, typical individual (new individual is based on
# variation among existing individuals)
gdd_typicalplant <- m_rema_gdd %>%
  epred_draws(newdata = expand_grid(agdd = seq(0, 2800, by = 50),
                                    individual_id = "new plant"),
              re_formula = NULL,
              allow_new_levels = TRUE,
              sample_new_levels = "uncertainty")
plot_gdd_typicalplant <- ggplot(filter(gdd_typicalplant, agdd %in% 0:1800), 
                                aes(x = agdd, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "AGDD (deg C)", y = "Predicted canopy fullness (%)",
       fill = "Credible interval", 
       title = "New, typical red maple") +
  theme_bw() +
  theme(legend.position = "bottom")
plot_gdd_typicalplant

# GDD model for multiple species in GRSM --------------------------------------#

# Attach AGDD values to GRSM data (2016-2024)
agdds <- weather %>%
  select(site_id, date, agdd) %>%
  rename(observation_date = date)

grsm <- grsm %>%
  left_join(agdds, by = c("site_id", "observation_date"))

# Standardized AGDD values and create proportion, factor variables
grsm <- grsm %>%
  # Removing 2020 data since it's more sparse (in case I want to compare with yr model)
  filter(yr != 2020) %>%
  mutate(agddz = (agdd - mean(agdd)) / sd(agdd),
         spp = factor(common_name),
         id = factor(individual_id),
         prop = intensity_midpoint / 95)

# Multi-species model with GDD instead of doy and year
m_grsm_gdd <- ordbetareg(prop ~ agddz + (1 + agddz|spp) + (1|id),
                         data = grsm, 
                         control = list(adapt_delta = 0.99),
                         iter = 4000, cores = 4, chains = 4,
                         backend = "cmdstanr")
#
# save(m_grsm_gdd, file = "output/grsm-gdd-model.RData")
#

# Expected canopy fullness for each species, ignoring individual REs
gdd_spp <- m_grsm_gdd %>%
  epred_draws(newdata = expand_grid(agddz = seq(-1.654, 3.1133, length = 100),
                                    spp = levels(grsm$spp)),
              re_formula = ~ (1 + agddz | spp))
plot_gdd_spp <- ggplot(filter(gdd_spp, agddz < 2), 
                       aes(x = agddz, y = .epred)) +
  stat_lineribbon(aes(color = spp, fill = after_scale(alpha(color, 0.2))), 
                  .width = 0.80) +
  # stat_lineribbon() +
  # scale_fill_brewer(palette = "Blues") +
  # facet_wrap(~ spp) +
  labs(x = "Day of year", y = "Predicted canopy fullness (%)", color = "Species",
       fill = "Species", title = "By species, correlated intercepts/slopes") +
  theme_bw() +
  theme(legend.position = "bottom")
plot_gdd_spp
# Don't have other models to compare with since this is the first time I used
# all years for GRSM sites, but the plot seems interesting. Sugar and striped
# maples early (very similar) then red maple, then basswood and red oak reach
# 50% (though basswood starts later and ramps up faster), then beech, which
# fills slowly.

# Multi-species model with doy and year
grsm <- grsm %>%
  mutate(fyr = factor(yr),
         doyz = (day_of_year - mean(day_of_year)) / sd(day_of_year)) 

start_time <- Sys.time()
m_grsm_doy <- ordbetareg(prop ~ doyz + fyr + (1 + doyz|spp) + (1|id),
                         data = grsm,
                         control = list(adapt_delta = 0.99),
                         iter = 4000, cores = 4, chains = 4,
                         backend = "cmdstanr")
end_time <- Sys.time()
end_time - start_time
save(m_grsm_doy, file = "output/grsm-doy-model.RData")
#
summary(m_grsm_doy)
#