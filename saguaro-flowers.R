# Floral abundance, saguaros
# ER Zylstra
# 28 October 2025

library(rnpn)
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(ordinal)
library(terra)

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

# Elevation
ggplot(sites, aes(x = elev)) +
  geom_histogram(bins = 20, width = 0.8)

# Geographic coordinates and plant-years
ggplot(sites, aes(x = lon, y = lat, size = n_plantyrs)) +
  geom_point()

# Process and attach weather data ---------------------------------------------#

# Daily data
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

weather <- weather %>%
  select(-c(lon, lat, elev)) %>%
  mutate(date = ymd(date),
         yr = year(date),
         month = month(date)) %>%
  group_by(site, yr, month) %>%
  summarize(tmin_mn = mean(tmin),
            tmax_mn = mean(tmax),
            ppt = sum(ppt),
            .groups = "keep") %>%
  filter(site %in% sites$site) %>%
  data.frame()

# 30-year (1991-2020) monthly normals:
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
weathervars <- weather %>%
  rename(tmin = tmin_mn,
         tmax = tmax_mn) %>%
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
            .groups = "keep") %>%
  filter(seasonyr %in% unique(flowers$yr)) %>%
  data.frame()

# Calculate 30-year normals
normvars <- norms %>%
  group_by(site) %>%
  summarize(tmin_wi_30 = mean(tmin30[month %in% c(12, 1:2)]),
            ppt_fa_30 = sum(ppt30[month %in% 9:11]),
            ppt_6_30 = sum(ppt30[month %in% c(9:12, 1:2)]),
            ppt_9_30 = sum(ppt30[month %in% c(6:12, 1:2)])) %>%
  data.frame()

# Merge normals with annual weather and calculate anomalies
weathervars <- weathervars %>%
  left_join(normvars, by = "site") %>%
  mutate(tmin_wi_a = tmin_wi - tmin_wi_30,
         ppt_fa_a = ppt_fa - ppt_fa_30,
         ppt_6_a = ppt_6 - ppt_6_30,
         ppt_9_a = ppt_9 - ppt_9_30)

# Add weather data to flowering dataset (and standardize)
flowers <- flowers %>%
  left_join(weathervars, by = c("site" = "site", "yr" = "seasonyr")) %>%
  mutate(tmin_wi_z = (tmin_wi_a - mean(tmin_wi_a)) / sd(tmin_wi_a),
         ppt_fa_z = (ppt_fa_a - mean(ppt_fa_a)) / sd(ppt_fa_a),
         ppt_6_z = (ppt_6_a - mean(ppt_6_a)) / sd(ppt_6_a),
         ppt_9_z = (ppt_9_a - mean(ppt_9_a)) / sd(ppt_9_a))

# Run ordinal models ----------------------------------------------------------#

# Standardize site elevations, so we can add that into models
flowers <- flowers %>%
  mutate(elev_z = (elev - mean(elev)) / sd(elev))


######## THINK IT MIGHT BE BETTER TO USE 3 COUNT CATEGORIES RATHER THAN 4
######## ADJUST MODELS BELOW ACCORDINGLY.
# sample sizes for 4 categories: 43, 20, 172, 125
# sample sizes for 3 categories: 63, 172, 125 (combining none + few)



# Model with year and plant, site random effects
m_yr <- clmm(abund4 ~ fyr + (1|id) + (1|site), Hess = TRUE, data = flowers)
summary(m_yr)

# Try different precip models
m_ppt_fa <- clmm(abund4 ~ ppt_fa_z + (1|id) + (1|site), Hess = TRUE,
                 data = flowers)
m_ppt_6 <- clmm(abund4 ~ ppt_6_z + (1|id) + (1|site), Hess = TRUE,
                data = flowers)
m_ppt_9 <- clmm(abund4 ~ ppt_9_z + (1|id) + (1|site), Hess = TRUE,
                data = flowers)
AIC(m_ppt_fa, m_ppt_6, m_ppt_9, m_yr)
# 9-mo ppt model better than other ppt models (though not as good as year)
summary(m_ppt_9)
# Negative effect of precipitation

# Winter minimum temperatures
m_tmin <- clmm(abund4 ~ tmin_wi_z + (1|id) + (1|site), Hess = TRUE,
               data = flowers)
summary(m_tmin) 
# Negative effect of winter temperatures

# Precipitation and winter temperatures
m_ppt9_tmin <- clmm(abund4 ~ ppt_9_z * tmin_wi_z + (1|id) + (1|site), 
                    Hess = TRUE, data = flowers)
summary(m_ppt9_tmin)

# Precipitation and winter temperatures, allowing effects to vary with elevation
m_3way <- clmm(abund4 ~ ppt_9_z * tmin_wi_z * elev_z + (1|id) + (1|site), 
               Hess = TRUE, data = flowers)
summary(m_3way)

# Need to use clmm2 to make predictions, and that function only allows one
# random effect:
m_3way2 <- clmm2(abund4 ~ ppt_9_z * tmin_wi_z * elev_z, random = factor(site), 
                 nAGQ = 10, Hess = TRUE, data = flowers)
summary(m_3way2)
# Higher-level interaction effects aren't significant, but still worth 
# including to highlight how weather effects vary regionally.

# Create dataframe for prediction
newdat <- expand.grid(
  abund4 = unique(flowers$abund4),
  elev_z = c(min(flowers$elev_z), 0, max(flowers$elev_z)),
  tmin_wi_z = c(min(flowers$tmin_wi_z), 0, max(flowers$tmin_wi_z)),
  ppt_9_z = seq(min(flowers$ppt_9_z), max(flowers$ppt_9_z), length = 100),
  KEEP.OUT.ATTRS = FALSE
)
# Plot predictions (= probability that saguaro will have X flowers
# given 9-month precipitation for cold/average/warm winter min temps at a 
# low/average/high elevation site)
preds <- cbind(newdat, 
               est = predict(m_3way2, newdata = newdat)) %>%
  arrange(abund4, elev_z, tmin_wi_z, ppt_9_z) %>%
  mutate(loc = case_when(
    elev_z == min(flowers$elev_z) ~ "low elev",
    elev_z == 0 ~ "average elev",
    elev_z == max(flowers$elev_z) ~ "high elev",
  )) %>%
  mutate(winter = case_when(
    tmin_wi_z == min(flowers$tmin_wi_z) ~ "cold winter",
    tmin_wi_z == 0 ~ "average winter",
    tmin_wi_z == max(flowers$tmin_wi_z) ~ "warm winter"
  )) %>%
  mutate(loc = factor(loc, 
                      levels = c("high elev", "average elev", "low elev"))) %>%
  mutate(winter = factor(winter, 
                         levels = c("cold winter", "average winter", 
                                    "warm winter")))

ggplot(preds, aes(x = ppt_9_z, y = est)) +
  geom_line(aes(color = abund4), linewidth = 1.3) +
  scale_color_brewer(palette = "BrBG") +
  facet_grid(loc ~ winter) +
  labs(x = "Cumulative 9-monthprecipitation, anomaly (standardized)", 
       y = "Probability", 
       color = "Flowers") +
  theme_bw() +
  theme(legend.position = "bottom")


m_noRE <- clm(abund4 ~ ppt_9_z * tmin_wi_z * elev_z, data = flowers)
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
    tmin_wi_z == min(flowers$tmin_wi_z) ~ "cold winter",
    tmin_wi_z == 0 ~ "average winter",
    tmin_wi_z == max(flowers$tmin_wi_z) ~ "warm winter"
  )) %>%
  mutate(loc = factor(loc, 
                      levels = c("high elev", "average elev", "low elev"))) %>%
  mutate(winter = factor(winter, 
                         levels = c("cold winter", "average winter", 
                                    "warm winter"))) %>%
  filter(loc == "high elev" & winter == "warm winter") %>%
  ggplot(aes(x = ppt_9_z, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = abund4), alpha = 0.3) +
  geom_line(aes(color = abund4), linewidth = 1.3) +
  scale_color_brewer(palette = "BrBG") +
  scale_fill_brewer(palette = "BrBG") +
  facet_wrap(~abund4) +
  labs(x = "Cumulative 9-monthprecipitation, anomaly (standardized)", 
       y = "Probability", 
       color = "Flowers") +
  theme_bw() +
  theme(legend.position = "bottom")


# Better fit with just 3 abundance categories?
m3_3way2 <- clmm2(abund3 ~ ppt_9_z * tmin_wi_z * elev_z, random = factor(site), 
                  nAGQ = 10, Hess = TRUE, data = flowers)
summary(m3_3way2)

# Create dataframe for prediction
newdat3 <- expand.grid(
  abund3 = unique(flowers$abund3),
  elev_z = c(min(flowers$elev_z), 0, max(flowers$elev_z)),
  tmin_wi_z = c(min(flowers$tmin_wi_z), 0, max(flowers$tmin_wi_z)),
  ppt_9_z = seq(min(flowers$ppt_9_z), max(flowers$ppt_9_z), length = 100),
  KEEP.OUT.ATTRS = FALSE
)
# Plot predictions (= probability that saguaro will have X flowers
# given 9-month precipitation for cold/average/warm winter min temps at a 
# low/average/high elevation site)
preds3 <- cbind(newdat3, 
               est = predict(m3_3way2, newdata = newdat3)) %>%
  arrange(abund3, elev_z, tmin_wi_z, ppt_9_z) %>%
  mutate(loc = case_when(
    elev_z == min(flowers$elev_z) ~ "low elev",
    elev_z == 0 ~ "average elev",
    elev_z == max(flowers$elev_z) ~ "high elev",
  )) %>%
  mutate(winter = case_when(
    tmin_wi_z == min(flowers$tmin_wi_z) ~ "cold winter",
    tmin_wi_z == 0 ~ "average winter",
    tmin_wi_z == max(flowers$tmin_wi_z) ~ "warm winter"
  )) %>%
  mutate(loc = factor(loc, 
                      levels = c("high elev", "average elev", "low elev"))) %>%
  mutate(winter = factor(winter, 
                         levels = c("cold winter", "average winter", 
                                    "warm winter")))

ggplot(preds3, aes(x = ppt_9_z, y = est)) +
  geom_line(aes(color = abund3), linewidth = 1.3) +
  scale_color_manual(values = c("#d8b365", "#80cdc1", "#018571")) +
  facet_grid(loc ~ winter) +
  labs(x = "Cumulative 9-month precipitation, anomaly (standardized)", 
       y = "Probability", 
       color = "Flowers") +
  theme_bw() +
  theme(legend.position = "bottom")

m3_noRE <- clm(abund3 ~ ppt_9_z * tmin_wi_z * elev_z, data = flowers)
preds3_noRE <- predict(m3_noRE, newdata = newdat3, interval = TRUE) 
preds3_noRE <- cbind(newdat3, fit = preds3_noRE$fit, 
                     lwr = preds3_noRE$lwr, upr = preds3_noRE$upr)

preds3_noRE %>%
  mutate(loc = case_when(
    elev_z == min(flowers$elev_z) ~ "low elev",
    elev_z == 0 ~ "average elev",
    elev_z == max(flowers$elev_z) ~ "high elev",
  )) %>%
  mutate(winter = case_when(
    tmin_wi_z == min(flowers$tmin_wi_z) ~ "cold winter",
    tmin_wi_z == 0 ~ "average winter",
    tmin_wi_z == max(flowers$tmin_wi_z) ~ "warm winter"
  )) %>%
  mutate(loc = factor(loc, 
                      levels = c("high elev", "average elev", "low elev"))) %>%
  mutate(winter = factor(winter, 
                         levels = c("cold winter", "average winter", 
                                    "warm winter"))) %>%
  filter(loc == "high elev" & winter == "warm winter") %>%
  ggplot(aes(x = ppt_9_z, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = abund3), alpha = 0.3) +
  geom_line(aes(color = abund3), linewidth = 1.3) +
  scale_color_manual(values = c("#d8b365", "#80cdc1", "#018571")) +
  scale_fill_manual(values = c("#d8b365", "#80cdc1", "#018571")) +
  facet_wrap(~abund3) +
  labs(x = "Cumulative 9-monthprecipitation, anomaly (standardized)", 
       y = "Probability", 
       color = "Flowers",
       fill = "Flowers") +
  theme_bw() +
  theme(legend.position = "bottom")
