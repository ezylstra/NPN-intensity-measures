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

m_rema1 <- ordbetareg(prop ~ day_of_year + (1|individual_id),
                      data = rema17, 
                      # true_bounds = c(0, 95),
                      # control = list(adapt_delta = 0.9),
                      cores = 4, chains = 4,
                      backend = "cmdstanr")
summary(m_rema1)
plot(m_rema1) # posterior distributions and trace plots

# For a good explanation of the different predictions types (grand mean,
# marginal effects), see:
# https://www.andrewheiss.com/blog/2021/11/10/ame-bayes-re-guide/

# Expected canopy fullness by date, ignoring individual effects (ie, grand mean)
doy_gm <- m_rema1 %>%
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
doy_typicalplant <- m_rema1 %>%
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
doy_newplant <- m_rema1 %>%
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
d <- as_draws_rvars(m_rema1)
d # Note that this will hang for a while to summarize random intercepts
  # Don't need to wait though - can just grab parameter names

intercept_p <- d$`b_Intercept`
slope_p <- d$`b_day_of_year`

x_50_p = (qlogis(prop50) - intercept_p) / slope_p
str(posterior::draws_of(x_50_p))
posterior::summarize_draws(x_50_p)
quantile(posterior::draws_of(x_50_p)[, 1], probs = c(0.025, 0.975))

# But what if we do want to take random effects into account?
# Need to add gaussian noise by sampling from normal(0, sd(indidividual_id))
samples <- coef(m_rema1, summary = FALSE) 
# This is a list of 1 array (name = individual_id). 
# Dims of array: rows = iter (4000); cols = individual (83); slice = param (2 [int, doy])

oint <- rowMeans(samples$individual_id[,,1])
# Summarizing the mean intercept across individuals for each iteration
# gets you very close to summary of Intercept in summary(m_rema1)

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
ggplot(data = filter(rema17, neon == 1),
       aes(x = day_of_year, y = intensity_midpoint)) +
  geom_line() +
  facet_wrap(~factor(individual_id))
# With just one year of data, it doesn't necessarily look like more sites 
# means more variation. Decent amount of variation in timing among plants 
# within NEON and Lumber Ridge sites. 

