# Floral abundance, saguaros
# ER Zylstra
# 29 October 2025

library(rnpn)
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(ordinal)
library(terra)
library(tidyterra)

# Load shapefile with US state boundaries -------------------------------------#

states <- vect("states/cb_2017_us_state_500k.shp")

# Download data, basic formatting (if not done already) -----------------------#

saguaro_data_file <- "npn-data/saguaro-flowers.csv"

if(!file.exists(saguaro_data_file)) {
  
  # Get species ID for saguaro
  saguaro_id <- npn_species() %>%
    filter(common_name == "saguaro") %>%
    pull(species_id)
  
  # Download most recent 10 years of data
  dl <- npn_download_status_data(
    request_source = "erinz",
    years = 2016:2025,
    species_ids = saguaro_id,
    climate_data = FALSE) %>%
    data.frame()
  
  # Extract data for flowers/flower buds
  df <- dl %>%
    filter(phenophase_description == "Flowers or flower buds")
  
  # There are some observations where the state field is missing. Will get state
  # for each site using shapefiles from the census bureau.
  state_fill <- df %>%
    distinct(site_id, latitude, longitude, state)
  state_fillv <- vect(state_fill, 
                      geom = c("longitude", "latitude"), 
                      crs = "epsg:4326")
  state_new <- terra::extract(states, state_fillv)
  state_fill <- cbind(state_fill, state_new = state_new$STUSPS)
  # Check how these compare to original assignments:
  count(state_fill, state, state_new)
  
  # Attach new state assignments to the saguaro data and restrict sites to Arizona
  df <- df %>%
    left_join(select(state_fill, site_id, state_new), by = "site_id") %>%
    select(-state) %>%
    rename(state = state_new) %>%
    filter(!is.na(state), state == "AZ")
  
  # Remove unnecessary columns
  df <- df %>%
    select(-c(update_datetime, genus, species, kingdom, phenophase_id,
              intensity_category_id, abundance_value)) 
  
  write.csv(df, saguaro_data_file, row.names = FALSE)
  rm(df, dl, state_fill, state_fillv, state_new, saguaro_id)
}

# Load saguaro flowering data and simplify ------------------------------------#

df <- read.csv(saguaro_data_file)

# Remove unnecessary columns, rename others, and create year column
df <- df %>%
  select(-c(observation_id, species_id)) %>%
  rename(site = site_id, 
         lat = latitude,
         lon = longitude, 
         elev = elevation_in_meters, 
         id = individual_id,
         php = phenophase_description, 
         obsdate = observation_date,
         doy = day_of_year,
         status = phenophase_status) %>%
  mutate(yr = year(obsdate))

# Remove observations where status was unknown (-1)
# sum(df$status == -1)/nrow(df) * 100 # 0.17% of observations (n = 44)
df <- filter(df, status != -1)

# Create numeric column with ~midpoints for each intensity category (for easier 
# sorting and coding)
df <- df %>%
  mutate(intensity_midpoint = case_when(
    intensity_value == "Less than 3" ~ 1,
    intensity_value == "3 to 10" ~ 5,
    intensity_value == "11 to 100" ~ 50,
    intensity_value == "101 to 1,000" ~ 500,
    intensity_value == "More than 1,000" ~ 1001,
    .default = NA
  ))
  
# Occasionally there are two records for a plant in one day with 
# different intensity or status values. Keeping record that was in phase with 
# highest intensity value (or record that doesn't have NAs)
df <- df %>%
  group_by(id, obsdate) %>%
  mutate(n_obs = n()) %>%
  ungroup() %>%
  data.frame() %>%
  arrange(id, obsdate, desc(status), desc(intensity_midpoint)) %>%
  distinct(id, obsdate, .keep_all = TRUE) %>%
  dplyr::select(-n_obs)

# Check that there aren't any observations where the plant was out of phase 
# but an intensity value was reported.
# count(df, status, intensity_midpoint)

# Make intensity midpoints = 0 if status = 0
df <- df %>%
  mutate(intensity_midpoint = ifelse(status == 0, 0, intensity_midpoint))
# Now the only NAs left in the intensity_midpoint column occur when the status 
# is yes but no intensity value was provided 

# Filter data -----------------------------------------------------------------#

# Want to exclude any saguaros that weren't observed multiple times during the 
# typical flowering period (since maximum flower counts could be biased low).

# Specifically, will identify the days-of-year that bracket the period when 80% 
# of the flowering observations occurred (flowering period), and the 
# days-of-year that bracket the period when 50% of the flowering occurred (peak
# flowering period). 

# Then:
# Exclude plant-year combos when <4 observations were made during the flowering
# period and <2 observations were made during the peak flowering period. 

# Identify when saguaros are typically flowering
flowering_periods <- df %>%
  filter(status == 1) %>%
  group_by(common_name) %>%
  summarize(nyrs = n_distinct(yr),
            nplants = n_distinct(id),
            nplantyrs = n_distinct(paste(yr, id)),
            mean = round(mean(doy)),
            q0.10 = quantile(doy, probs = 0.10),
            q0.25 = quantile(doy, probs = 0.25),
            q0.50 = quantile(doy, probs = 0.50),
            q0.75 = quantile(doy, probs = 0.75),
            q0.90 = quantile(doy, probs = 0.90),
            .groups = "keep") %>%
  mutate(across(mean:q0.90, round)) %>%
  mutate(length_50 = q0.75 - q0.25,
         length_80 = q0.90 - q0.10) %>%
  data.frame()
flowering_periods

q0.10 <- flowering_periods$q0.10[1]
q0.25 <- flowering_periods$q0.25[1]
q0.75 <- flowering_periods$q0.75[1]
q0.90 <- flowering_periods$q0.90[1]

# Summarize amount and quality of information for each plant and year
pl_yr <- df %>%
  mutate(in50 = ifelse(doy >= q0.25 & doy <= q0.75, 1, 0)) %>%
  mutate(in80 = ifelse(doy >= q0.10 & doy <= q0.90, 1, 0)) %>%
  group_by(id, site, yr) %>%
  summarize(nobs = n(),
            first_obs = min(doy),
            last_obs = max(doy),
            nobs50 = sum(in50),
            nobs80 = sum(in80),
            n_inphase = sum(status),
            n_intvalue = sum(!is.na(intensity_value)),
            prop_inphase = round(n_inphase / nobs, 2),
            prop_intvalue = round(n_intvalue / n_inphase, 2),
            max_intensity = ifelse(sum(is.na(intensity_midpoint)) == n(),
                                   NA, max(intensity_midpoint, na.rm = TRUE)),
            .groups = "keep") %>%
  data.frame()

# Identify which plant-yr combinations to filter out
pl_yr <- pl_yr %>%
  mutate(remove = case_when(
    is.na(max_intensity) ~ 1,
    nobs80 < 4 ~ 1,
    nobs50 < 2 ~ 1,
    .default = 0
  ))
# What proportion of plant-years will be included/excluded and what do the 
# distribution of max counts look like in included/excluded plant-years?
pl_yr %>% 
  group_by(remove) %>%
  summarize(n = n(),
            nobs_mn = median(nobs),
            nobs50_md = median(nobs50),
            nobs80_md = median(nobs80),
            max0_n = sum(max_intensity == 0 & !is.na(max_intensity)),
            max10_n = sum(max_intensity %in% c(1, 5) & !is.na(max_intensity)),
            max50_n = sum(max_intensity == 50 & !is.na(max_intensity)),
            max500_n = sum(max_intensity >= 500 & !is.na(max_intensity))) %>%
  data.frame()

# Filter out those plant-years from dataset
df_filter <- df %>%
  left_join(select(pl_yr, id, yr, remove), by = c("id", "yr")) %>%
  filter(remove == 0) %>%
  select(-remove)

# Aggregate data for each plant-year ------------------------------------------#

flowers <- df_filter %>%
  group_by(id, site, lat, lon, elev, yr) %>%
  summarize(n_obs = n(),
            n_inphase = sum(status),
            max_count = max(intensity_midpoint, na.rm = TRUE),
            .groups = "keep") %>%
  data.frame()
head(flowers)

count(flowers, yr) # 11-53 saguaros/yr
count(flowers, max_count)

# Create four count categories
  # midpoint = 0 (none)
  # midpoints = 1 or 5 (1-10; few)
  # midpoint = 50 (11-100; some)
  # midpoints = 500 or 1001 (> 100; many)
flowers <- flowers %>%
  mutate(abund4 = case_when(
    max_count == 0 ~ "none",
    max_count %in% 1:5 ~ "few",
    max_count == 50 ~ "some",
    max_count >= 500 ~ "many",
  )) %>%
  mutate(abund4 = factor(abund4, 
                         levels = c("none", "few", "some", "many"),
                         ordered = TRUE))

# Create three count categories
  # midpoint = 0 or 1 or 5 (0-10; few or none)
  # midpoint = 50 (11-100; some)
  # midpoints = 500 or 1001 (> 100; many)
flowers <- flowers %>%
  mutate(abund3 = case_when(
    max_count %in% 0:5 ~ "few or none",
    max_count == 50 ~ "some",
    max_count >= 500 ~ "many",
  )) %>%
  mutate(abund3 = factor(abund3, 
                         levels = c("few or none", "some", "many"),
                         ordered = TRUE))

# Convert year variable as a factor
flowers$fyr <- factor(flowers$yr)

# Extract information about sites ---------------------------------------------#

sites <- flowers %>%
  group_by(site, lat, lon, elev) %>%
  summarize(n_plants = n_distinct(id),
            n_yrs = n_distinct(yr),
            n_plantyrs = n_distinct(paste0(id, "_", yr)),
            .groups = "keep") %>%
  data.frame()

# Write sites information to file to get PRISM weather data
# Downloading PRISM data through web interface, so we could get data at 800m 
# resolution without having to download daily rasters for CONUS, which takes a
# long time (https://prism.oregonstate.edu/explorer/bulk.php)
# write.table(select(sites, lat, lon, site), "weather-data/saguaro-sites.csv", 
#             sep = ",", row.names = FALSE, col.names = FALSE)

# Process and attach weather data ---------------------------------------------#

# Load daily data
prism_files <- list.files("weather-data/saguaros",
                          full.names = TRUE,
                          pattern = "stable|provisional")
for (i in 1:length(prism_files)) {
  prism1 <- read.csv(prism_files[i],
                     header = FALSE,
                     skip = 11,
                     col.names = c("site", "lon", "lat", "elev",
                                   "date", "ppt", "tmin", "tmax"))
  if (i == 1) {
    weather <- prism1
  } else {
    weather <- rbind(weather, prism1)
  }
}
rm(prism1)

# Summarize by month
weather <- weather %>%
  select(-c(lon, lat, elev)) %>%
  mutate(tmean = (tmin + tmax) / 2,
         date = ymd(date),
         yr = year(date),
         month = month(date)) %>%
  group_by(site, yr, month) %>%
  summarize(tmin_mn = mean(tmin),
            tmax_mn = mean(tmax),
            tmean_mn = mean(tmean),
            ppt = sum(ppt),
            .groups = "keep") %>%
  filter(site %in% sites$site) %>%
  data.frame()

# Load 30-year (1991-2020) monthly normals:
normsfile <- list.files("weather-data/saguaros/", 
                        full.names = TRUE,
                        pattern = "30yr")
norms <- read.csv(normsfile,
                  header = FALSE,
                  skip = 11,
                  col.names = c("site", "lon", "lat", "elev",
                                "mon", "ppt30", "tmin30", "tmax30")) %>%
  filter(mon != "Annual") %>%
  mutate(month = as.integer(factor(mon, levels = month.name))) %>%
  select(-c(lon, lat, elev, mon)) %>%
  filter(site %in% sites$site)

# Weather variables:
  # Min temps in winter (DJF) [mean of daily minimum temperatures]
  # Precip in fall (SON)
  # Precip, 6 mo (fall + winter; Sep-Feb)
  # Precip, 9 mo (summer + fall + winter; Jun-Feb)
  # Annual mean temperature
  # Annual precipitation
weathervars <- weather %>%
  rename(tmin = tmin_mn,
         tmax = tmax_mn,
         tmean = tmean_mn) %>%
  # Create seasonyr to match up with flowering year
  mutate(seasonyr = ifelse(month %in% 6:12, yr + 1, yr)) %>%
  mutate(season = case_when(
    month %in% 3:5 ~ "sp",
    month %in% 6:8 ~ "su",
    month %in% 9:11 ~ "fa",
    .default = "wi"
  )) %>%
  group_by(site, seasonyr) %>%
  summarize(tmin_wi = mean(tmin[season == "wi"]),
            ppt_fa = sum(ppt[season == "fa"]),
            ppt_6 = sum(ppt[season %in% c("fa", "wi")]),
            ppt_9 = sum(ppt[season != "sp"]),
            ann_temp = mean(tmean),
            ann_ppt = sum(ppt),
            .groups = "keep") %>%
  filter(seasonyr %in% unique(flowers$yr)) %>%
  data.frame()

# Calculate 30-year normals
normvars <- norms %>%
  mutate(tmean30 = (tmin30 + tmax30) / 2) %>%
  group_by(site) %>%
  summarize(tmin_wi_30 = mean(tmin30[month %in% c(12, 1:2)]),
            ppt_fa_30 = sum(ppt30[month %in% 9:11]),
            ppt_6_30 = sum(ppt30[month %in% c(9:12, 1:2)]),
            ppt_9_30 = sum(ppt30[month %in% c(6:12, 1:2)]),
            ann_temp_30 = mean(tmean30),
            ann_ppt_30 = sum(ppt30)) %>%
  data.frame()

# Merge normals with annual weather and calculate anomalies for temperature
# and % of 30-year normals for precipitation. (If used anomalies for precip
# then you see the range of anomalies is smaller for low elevation/drier sites)
weathervars <- weathervars %>%
  left_join(normvars, by = "site") %>%
  mutate(tmin_wi_a = tmin_wi - tmin_wi_30,
         ppt_fa_p = (ppt_fa / ppt_fa_30) * 100,
         ppt_6_p = (ppt_6 / ppt_6_30) * 100,
         ppt_9_p = (ppt_9 / ppt_9_30) * 100)

# Add climate normals to sites dataframe amd look at patterns with elevation
sites <- sites %>%
  left_join(normvars, by = "site")
ggplot(sites) +
  geom_point(aes(x = ann_temp_30, y = ann_ppt_30, color = elev),
             size = 2) +
  scale_color_viridis_c() +
  labs(x = "Mean annual temperature (degC)", 
       y = "Mean annual precipitation (mm)",
       color = "Elevation (m)") +
  theme_bw()

# Will remove single site that is 350 m higher elevation than the rest 
# because climate there is much different (wetter and colder). Site only had
# 2 years of data for one saguaro anyways.
sites <- sites %>%
  filter(elev < 1100)
flowers <- flowers %>%
  filter(elev < 1100)

# Add weather data to flowering dataset (and standardize weather, elev data)
flowers <- flowers %>%
  left_join(weathervars, by = c("site" = "site", "yr" = "seasonyr")) %>%
  mutate(tmin_wi_z = (tmin_wi_a - mean(tmin_wi_a)) / sd(tmin_wi_a),
         ppt_fa_z = (ppt_fa_p - mean(ppt_fa_p)) / sd(ppt_fa_p),
         ppt_6_z = (ppt_6_p - mean(ppt_6_p)) / sd(ppt_6_p),
         ppt_9_z = (ppt_9_p - mean(ppt_9_p)) / sd(ppt_9_p),
         elev_z = (elev - mean(elev)) / sd(elev))

# Plot temp, precip z-scores by elevation
flowers %>%
  mutate(elev_cats = cut(elev, breaks = c(0, 600, 900, 1100),
                         labels = c("<600", "600-900", ">900"))) %>%
  ggplot(aes(x = ppt_9_z, y = tmin_wi_z, color = elev)) +
  geom_point() +
  facet_grid(~elev_cats) +
  scale_color_viridis_c() +
  theme_bw()
# All sites had average - hotter than average winters, but only the mid-elev
# sites had cooler than average winters.
# At the high elevation sites, there were no years when precipitation was 
# much lower than average. At all other sites, precipitation values look good.
# This has implications for prediction plot....

# Elevation/location figures --------------------------------------------------#

# Elevation
ggplot(sites, aes(x = elev)) +
  geom_histogram(bins = 30, fill = "gray") +
  labs(x = "Elevation (m)", y = "No. sites") +
  geom_vline(aes(xintercept = mean(elev)), color = "steelblue4",
             linetype = "dashed") +
  geom_vline(aes(xintercept = 1074), color = "salmon4", linetype = "dashed") +
  annotate(geom = "text", label = "Mean", color = "steelblue4", 
           x = mean(sites$elev) - 5, y = 9, hjust = 1) +
  annotate(geom = "text", label = "Renzi", color = "salmon4", 
           x = 1074 - 5, y = 9, hjust = 1) +
  theme_bw()

# Geographic coordinates and plant-years
ggplot(filter(states, STUSPS == "AZ")) +
  geom_spatvector(fill = "transparent") +
  geom_spatvector(data = vect(sites, geom = c("lon", "lat"), crs = "epsg:4326"),
                  aes(size = n_plantyrs), color = "steelblue4", alpha = 0.5) +
  labs(size = "No. plant-years") +
  theme_bw()
# Would be better to have this overlaid on an elevation gradient...

# Run ordinal models ----------------------------------------------------------#

# Using 3 count categories because then sample sizes are a bit more equitable

# Because we'll ultimately want to use clmm2 models to make predictions (if we 
# don't go a Bayesian route), we're restricted to one random effect. 

# Model with year (not trend) effects and site-level random effects
m_yr_site <- clmm2(abund3 ~ fyr, random = factor(site), 
                   Hess = TRUE, data = flowers)
summary(m_yr_site)
# CIs from the profiled likelihood for the SD for the random effect
confint(m_yr_site)

# Model with year (not trend) effects and individual-level random effects
m_yr_ind <- clmm2(abund3 ~ fyr, random = factor(id), 
                   Hess = TRUE, data = flowers)
summary(m_yr_ind)
# CIs from the profiled likelihood for the SD for the random effect
confint(m_yr_ind)

################################################################################
# Is it worth going Bayesian so we can evaluate different random effects?

library(brms)
m_test <- brm(abund3 ~ fyr + (1|site) + (1|id), data = flowers,
              family = cumulative(link = "logit"))
summary(m_test)
m_test_s <- brm(abund3 ~ fyr + (1|site), data = flowers,
              family = cumulative(link = "logit"))
summary(m_test_s)
m_test_i <- brm(abund3 ~ fyr + (1|id), data = flowers,
                family = cumulative(link = "logit"))
summary(m_test_i)

loo_both <- loo(m_test, cores = 4)
loo_both
loo_site <- loo(m_test_s, cores = 4)
loo_site
loo_id <- loo(m_test_i, cores = 4)
loo_id
loo_compare(loo_both, loo_site, loo_id)
# Model with elpd_diff = 0 and/or lowest looic is best
# And this suggests that a model that includes both random effects is much
# better than a model with just site or just individual. 

m_test_full <- brm(abund3 ~ ppt_9_z * tmin_wi_z * elev_z + (1|site) + (1|id),
                   data = flowers, family = cumulative(link = "logit"))
summary(m_test_full)
# Estimates are pretty similar to those from a clmm2 model with site random effects.
################################################################################

# Try different precip models
m_ppt_fa <- clmm2(abund3 ~ ppt_fa_z, random = factor(site), Hess = TRUE,
                  data = flowers)
m_ppt_6 <- clmm2(abund3 ~ ppt_6_z, random = factor(site), Hess = TRUE,
                 data = flowers)
m_ppt_9 <- clmm2(abund3 ~ ppt_9_z, random = factor(site), Hess = TRUE,
                 data = flowers)
AIC(m_ppt_fa, m_ppt_6, m_ppt_9, m_yr)
# 9-mo ppt model better than other ppt models (though not as good as year)
summary(m_ppt_9)
# Negative effect of precipitation

# Winter minimum temperatures
m_tmin <- clmm2(abund3 ~ tmin_wi_z, random = factor(site), Hess = TRUE,
                data = flowers)
summary(m_tmin) 
# Negative effect of winter temperatures

# Precipitation and winter temperatures
m_ppt9_tmin <- clmm2(abund3 ~ ppt_9_z * tmin_wi_z, random = factor(site), 
                     Hess = TRUE, data = flowers)
summary(m_ppt9_tmin)
# All significant, negative effects (including interaction)

# Precipitation and winter temperatures, allowing effects to vary with elevation
m_full <- clmm2(abund3 ~ ppt_9_z * tmin_wi_z * elev_z, random = factor(site), 
                Hess = TRUE, data = flowers)
summary(m_full)
# Most interaction effects that include elevation aren't significant, but at
# one is close (elev*temp; P = 0.08). Still worth including elevation to 
# highlight how weather effects vary regionally.
anova(m_full, m_ppt9_tmin)
# Chi-squared not signif (P = 0.19), but -logLik lower with full model

# Precipitation and winter temperatures, allowing effects to vary with elevation
# With individual, not site-level, random effects
m_full_ind <- clmm2(abund3 ~ ppt_9_z * tmin_wi_z * elev_z, random = factor(id), 
                    Hess = TRUE, data = flowers)
summary(m_full_ind)
AIC(m_full_ind, m_full)

# Make predictions ------------------------------------------------------------#

# Won't make predictions for colder than average winters since that only 
# occurred at middle elevation sites (not high or low elev). 

# For high elevation sites, we could restrict precip z-scores to observed range 
# since there were no observations in abnormally dry years, but this does make
# things look a little odd. Maybe keep in (and extrapolate for now), but
# 

# See ranges in weather z-scores
summary(flowers$tmin_wi_z)
summary(flowers$ppt_9_z)
highelev_ppt <- range(flowers$ppt_9_z[flowers$elev > 900]) %>% round(1)

# Create dataframe for prediction
newdat <- expand.grid(
  abund3 = unique(flowers$abund3),
  elev_z = c(min(flowers$elev_z), 0, max(flowers$elev_z)),
  tmin_wi_z = c(0, 2),
  ppt_9_z = seq(min(flowers$ppt_9_z), max(flowers$ppt_9_z), length = 100),
  KEEP.OUT.ATTRS = FALSE
)
# Plot predictions (= probability that saguaro will have X flowers
# given 9-month precipitation for cold/average/warm winter min temps at a 
# low/average/high elevation site)
preds <- cbind(newdat, 
               est = predict(m_full, newdata = newdat)) %>%
  arrange(abund3, elev_z, tmin_wi_z, ppt_9_z) %>%
  mutate(ppt_9_p = ppt_9_z * sd(flowers$ppt_9_p) + mean(flowers$ppt_9_p)) %>%
  mutate(loc = case_when(
    elev_z == min(flowers$elev_z) ~ "Low elevation",
    elev_z == 0 ~ "Average elevation",
    elev_z == max(flowers$elev_z) ~ "High elevation",
  )) %>%
  mutate(winter = case_when(
    tmin_wi_z == 0 ~ "Average winter",
    tmin_wi_z == 2 ~ "Warm winter"
  )) %>%
  mutate(loc = factor(loc, 
                      levels = c("High elevation", 
                                 "Average elevation", 
                                 "Low elevation"))) %>%
  mutate(winter = factor(winter, 
                         levels = c("Average winter", "Warm winter"))) %>%
  # filter(!(loc == "High elevation" & 
  #            (ppt_9_z < highelev_ppt[1] | ppt_9_z > highelev_ppt[2]))) %>%
  mutate(abund3 = case_when(
    abund3 == "few or none" ~ "10 or less",
    abund3 == "some" ~ "11 to 100",
    abund3 == "many" ~ "More than 100"
  )) %>%
  mutate(abund3 = factor(abund3, levels = c("10 or less",
                                            "11 to 100",
                                            "More than 100")))

fig_6panel <- ggplot(preds, aes(x = ppt_9_p, y = est)) +
  geom_line(aes(color = abund3), linewidth = 1.3) +
  scale_color_manual(values = c("#d8b365", "#80cdc1", "#018571")) +
  facet_grid(loc ~ winter) +
  labs(x = "Cumulative 9-month precipitation, % of normal", 
       y = "Probability", 
       color = "Flowers") +
  theme_bw() +
  theme(legend.position = "bottom")

# ggsave("output/manuscript/saguaro-predictions-6panel.png",
#        fig_6panel,
#        height = 8, width = 6.5, units = "in", dpi = 600)

fig_4panel <- ggplot(filter(preds, loc != "Average elevation"),
                     aes(x = ppt_9_p, y = est)) +
  geom_line(aes(color = abund3), linewidth = 1.3) +
  scale_color_manual(values = c("#d8b365", "#80cdc1", "#018571")) +
  facet_grid(loc ~ winter) +
  labs(x = "Cumulative 9-month precipitation, % of normal", 
       y = "Probability", 
       color = "Flowers") +
  theme_bw() +
  theme(legend.position = "bottom")

# ggsave("output/manuscript/saguaro-predictions-4panel.png",
#        fig_4panel,
#        height = 5.5, width = 6.5, units = "in", dpi = 600)

# Other model/figure options --------------------------------------------------#

# Could also remove non-signficant terms and make predictions from simpler model
m_simple <- clmm2(abund3 ~ ppt_9_z * tmin_wi_z  + tmin_wi_z * elev_z, 
                  random = factor(site), 
                  Hess = TRUE, data = flowers)
summary(m_simple)

preds_simple <- cbind(newdat, 
                      est = predict(m_simple, newdata = newdat)) %>%
  arrange(abund3, elev_z, tmin_wi_z, ppt_9_z) %>%
  mutate(ppt_9_p = ppt_9_z * sd(flowers$ppt_9_p) + mean(flowers$ppt_9_p)) %>%
  mutate(loc = case_when(
    elev_z == min(flowers$elev_z) ~ "Low elevation",
    elev_z == 0 ~ "Average elevation",
    elev_z == max(flowers$elev_z) ~ "High elevation",
  )) %>%
  mutate(winter = case_when(
    tmin_wi_z == 0 ~ "Average winter",
    tmin_wi_z == 2 ~ "Warm winter"
  )) %>%
  mutate(loc = factor(loc, 
                      levels = c("High elevation", 
                                 "Average elevation", 
                                 "Low elevation"))) %>%
  mutate(winter = factor(winter, 
                         levels = c("Average winter", "Warm winter"))) %>%
  # filter(!(loc == "High elevation" & 
  #            (ppt_9_z < highelev_ppt[1] | ppt_9_z > highelev_ppt[2]))) %>%
  mutate(abund3 = str_to_sentence(abund3)) %>%
  mutate(abund3 = factor(abund3, levels = c("Few or none", "Some", "Many")))

ggplot(preds_simple, aes(x = ppt_9_p, y = est)) +
  geom_line(aes(color = abund3), linewidth = 1.3) +
  scale_color_manual(values = c("#d8b365", "#80cdc1", "#018571")) +
  facet_grid(loc ~ winter) +
  labs(x = "Cumulative 9-month precipitation, % of normal", 
       y = "Probability", 
       color = "Flowers") +
  theme_bw() +
  theme(legend.position = "bottom")

# Use full model, but plot predictions for one set of conditions with CIs
m_noRE <- clm(abund3 ~ ppt_9_z * tmin_wi_z * elev_z, data = flowers)
preds_noRE <- predict(m_noRE, newdata = newdat, interval = TRUE) 
preds_noRE <- cbind(newdat, fit = preds_noRE$fit, 
                    lwr = preds_noRE$lwr, upr = preds_noRE$upr)

preds_noRE %>%
  mutate(loc = case_when(
    elev_z == min(flowers$elev_z) ~ "low elev",
    elev_z == 0 ~ "average elev",
    elev_z == max(flowers$elev_z) ~ "high elev",
  )) %>%
  mutate(winter = case_when(
    tmin_wi_z == 0 ~ "average winter",
    tmin_wi_z == 2 ~ "warm winter"
  )) %>%
  mutate(loc = factor(loc, 
                      levels = c("high elev", "average elev", "low elev"))) %>%
  mutate(winter = factor(winter, 
                         levels = c("average winter", 
                                    "warm winter"))) %>%
  filter(loc == "high elev" & winter == "warm winter") %>%
  mutate(abund3 = str_to_sentence(abund3)) %>%
  mutate(abund3 = factor(abund3, levels = c("Few or none", "Some", "Many"))) %>%
  ggplot(aes(x = ppt_9_z, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = abund3), alpha = 0.3) +
  geom_line(aes(color = abund3), linewidth = 1.3) +
  scale_color_manual(values = c("#d8b365", "#80cdc1", "#018571")) +
  scale_fill_manual(values = c("#d8b365", "#80cdc1", "#018571")) +
  facet_grid(abund3 ~ .) +
  labs(x = "Cumulative 9-month precipitation, % 30-yr normals", 
       y = "Probability", 
       color = "Flowers") +
  theme_bw() +
  theme(legend.position = "bottom")
