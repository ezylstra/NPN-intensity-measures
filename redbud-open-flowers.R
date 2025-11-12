# Number of open flowers, eastern redbud
# ER Zylstra
# 12 November 2025

library(rnpn)
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(lme4)
library(ggplot2)
library(terra)

# library(sf)

# Load shapefile with US state boundaries -------------------------------------#

states <- vect("states/cb_2017_us_state_500k.shp")

# Download data, basic formatting (if not done already) -----------------------#

redbud_data_file <- "npn-data/redbud-flowers.csv"

if(!file.exists(redbud_data_file)) {
  
  # Get species ID for eastern redbud
  redbud_id <- npn_species() %>%
    filter(common_name == "eastern redbud") %>%
    pull(species_id)

  # Download data since 2016
  dl <- npn_download_status_data(
    request_source = "erinz",
    years = 2016:2025,
    species_ids = redbud_id,
    climate_data = FALSE,
    additional_fields = "observedby_person_id") %>%
    data.frame()
  
  # Extract data for flowers, open flowers phenophases
  df <- dl %>%
    filter(phenophase_description %in% c("Flowers or flower buds", 
                                         "Open flowers"))
  
  # There are some observations where the state field is missing. Will get state
  # for each missing site using shapefiles from the census bureau.
  state_fill <- df %>%
    distinct(site_id, latitude, longitude, state)
  state_fillv <- vect(state_fill, 
                      geom = c("longitude", "latitude"), 
                      crs = "epsg:4326")
  state_new <- terra::extract(states, state_fillv)
  state_fill <- cbind(state_fill, state_new = state_new$STUSPS)

  # Attach new state assignments to the redbud data
  df <- df %>%
    left_join(select(state_fill, site_id, state_new), by = "site_id") %>%
    mutate(state = ifelse(is.na(state), state_new, state)) %>%
    select(-state_new)

  # Limit to US + Ontario and longitude > 100 degW
  df <- df %>%
    filter(!is.na(state) & state != "QC") %>%
    filter(longitude > (-100))

  # Remove unnecessary columns
  df <- df %>%
    select(-c(update_datetime, genus, species, kingdom, phenophase_id,
              intensity_category_id, abundance_value)) 
  
  write.csv(df, redbud_data_file, row.names = FALSE)
  rm(df, dl, redbud_id, state_fill, state_fillv, state_new)
}

# Load redbud flowers data and simplify ---------------------------------------#

df <- read.csv(redbud_data_file)

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
# sum(df$status == -1)/nrow(df) * 100 # 0.37% of observations (n = )
df <- filter(df, status != -1)

# Create numeric column with ~midpoints for each intensity category (for easier 
# sorting and coding)
df <- df %>%
  mutate(intensity_midpoint = case_when(
    intensity_value == "Less than 5%" ~ 2,
    intensity_value == "5-24%" ~ 14,
    intensity_value == "25-49%" ~ 37,
    intensity_value == "50-74%" ~ 62,
    intensity_value == "75-94%" ~ 84,
    intensity_value == "95% or more" ~ 95,
    intensity_value == "Less than 3" ~ 1,
    intensity_value == "3 to 10" ~ 5,
    intensity_value == "11 to 100" ~ 50,
    intensity_value == "101 to 1,000" ~ 500,
    intensity_value == "1,001 to 10,000" ~ 5000,
    intensity_value == "More than 10,000" ~ 10001,
    .default = NA
  ))

# Occasionally there are two records for a plant-phenophase in one day with 
# different intensity or status values. Keeping record that was in phase with 
# highest intensity value (or record that doesn't have NAs)
df <- df %>%
  group_by(id, php, obsdate) %>%
  arrange(id, obsdate, php, desc(status), desc(intensity_midpoint)) %>%
  distinct(id, php, obsdate, .keep_all = TRUE) %>%
  data.frame()

# Are there observations where the plant is out of phase but an intensity value
# is reported? If so, change the status but add a column to note that the data
# were amended.
df <- df %>%
  mutate(amended_status = ifelse(status == 0 & !is.na(intensity_midpoint), 
                                 1, 0)) %>%
  mutate(status = ifelse(status == 0 & !is.na(intensity_midpoint), 1, status))

# Extract information about sites ---------------------------------------------#

sites <- df %>%
  group_by(site, lat, lon, elev, state) %>%
  summarize(n_plants = n_distinct(id),
            n_yrs = n_distinct(yr),
            n_plantyrs = n_distinct(paste0(id, "_", yr)),
            .groups = "keep") %>%
  data.frame()

# Combine observations from both flower phenophases ---------------------------#

# Simplify data
flowers <- df %>%
  mutate(php = ifelse(php == "Open flowers", "open", "flowers")) %>%
  select(-c(lat, lon, elev, state, common_name, amended_status)) %>%
  rename(intensity = intensity_midpoint,
         observer = observedby_person_id)

# Put observations of all phenophases by same observer on same day in one row
flowers <- flowers %>%
  pivot_wider(id_cols = c(observer, site, id, obsdate, doy, yr),
              names_from = php,
              values_from = c(status, intensity)) %>%
  data.frame()

# Look at instances where open = 1 and flowers = 0 or NA.
# If there's an intensity value for open flowers, then we'll assume that 
# flower status should be 1. Otherwise we'll delete the observation
# filter(flowers, status_open == 1) %>%
#   count(status_flowers, status_open, !is.na(intensity_open))
flowers <- flowers %>%
  mutate(status_flowers = case_when(
    status_open == 1 & 
      (is.na(status_flowers) | status_flowers == 0) & !is.na(intensity_open) ~ 1,
    .default = status_flowers
  )) %>%
  filter(!(!is.na(status_open) & status_open == 1 &
             (is.na(status_flowers) | status_flowers == 0)))

# Remove observations with status_flowers = NA & status_open = 0 (no info)
flowers <- flowers %>%
  filter(!(is.na(status_flowers) & !is.na(status_open) & status_open == 0))
# Now there are no more status_flowers = NA

# If status_flowers == 0, change all status_open to 0
flowers <- flowers %>%
  mutate(status_open = ifelse(status_flowers == 0, 0, status_open))

# Convert intensity values to 0 when status = 0 
flowers <- flowers %>%
  mutate(intensity_flowers = ifelse(status_flowers == 0, 0, intensity_flowers)) %>%
  mutate(intensity_open = ifelse(!is.na(status_open) & status_open == 0,
                                 0, intensity_open))

# Remove observations when both intensity values are NA
flowers <- flowers %>%
  filter(!(is.na(intensity_flowers) & is.na(intensity_open)))

# Look at totals
count(flowers, status_flowers, status_open, intensity_flowers, intensity_open)
# Now, all  status_flowers = 0 or 1 (no NAs)
# There are 204 observations with status_flowers = 1 and status_open = NA
# All intensity values = 0 when status = 0
# NA intensity values only present when status = 1 and the other intensity value is present

# Are there any observations of the same plant on the same day?
flowers %>% distinct(id, obsdate) %>% nrow() == nrow(flowers)
# Yes

# flowers %>%
#   group_by(id, obsdate) %>%
#   summarize(n_obs = n(), 
#             n_status_open = sum(!is.na(status_open)),
#             n_int_flowers = sum(!is.na(intensity_flowers)),
#             n_int_open = sum(!is.na(intensity_open)),
#             .groups = "keep") %>%
#   data.frame() %>%
#   filter(n_obs > 1) 

# In almost every case, (n = 21) one observation provided an intensity values
# for flowers and the other provided an intensity value for open flowers.

# Will combine them.
max_na <- function(x) {
  ifelse(sum(is.na(x)) == length(x), NA, max(x, na.rm = TRUE))
}
flowers <- flowers %>%
  group_by(site, id, obsdate, doy, yr) %>%
  summarize(status_flowers = max_na(status_flowers),
            status_open = max_na(status_open),
            intensity_flowers = max_na(intensity_flowers),
            intensity_open = max_na(intensity_open),
            .groups = "keep") %>%
  data.frame()

# Plot intensity values -------------------------------------------------------#

# How much data each year?
flowers %>%
  group_by(yr) %>%
  summarize(n_plants = n_distinct(id),
            n_sites = n_distinct(site),
            n_obs = n(),
            n_int_flowers = sum(!is.na(intensity_flowers)),
            n_int_open = sum(!is.na(intensity_open)),
            n_int_both = sum(!is.na(intensity_flowers) & 
                               !is.na(intensity_open))) %>%
  data.frame()
# Way more data in 2022-2025 than previous years

# Look at the weekly mean number of flowers and open flowers by year (calculated
# open flowers coarsely by multiplying averages across individuals, not within 
# individuals). Can see peak in open flowers does not always coincide with 
# flowers.
wkcounts <- flowers %>%
  mutate(wk = as.numeric(week(obsdate))) %>%
  group_by(yr, wk) %>%
  summarize(flowers = mean(intensity_flowers, na.rm = TRUE),
            open = mean(intensity_open, na.rm = TRUE),
            .groups = "keep") %>%
  mutate(nopen = ceiling(flowers * open / 100)) %>%
  pivot_longer(cols = c(flowers, nopen),
               names_to = "type",
               values_to = "intensity") %>%
  mutate(type = ifelse(type == "nopen", "Open flowers", "Flowers")) %>%
  filter(wk %in% 9:25)
wkcount_summary <- wkcounts %>%
  group_by(yr, type) %>%
  summarize(max_value = max(intensity),
            max_wk = mean(wk[intensity == max_value]),
            .groups = "keep") %>%
  data.frame()

wkcounts_plot <- ggplot(data = wkcounts, aes(x = wk, y = intensity)) +
  geom_line(aes(color = type, linetype = type)) +
  scale_color_manual(values = c("#66c2a5", "#fc8d62")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  geom_segment(data = wkcount_summary, 
               aes(x = max_wk, xend = max_wk, y = 0, yend = max_value,
                   linetype = type),
               show.legend = FALSE,
               color = "gray30") +
  facet_wrap(~yr) +
  labs(x = "Week", y = "Mean count") +
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.85, 0.15),
        legend.title = element_blank(),
        panel.grid = element_blank())
wkcounts_plot
# ggsave("output/manuscript/redbud-weekly-mn-counts-3x4.png",
#        wkcounts_plot,
#        width = 6.5, height = 5, units = "in", dpi = 600)

wkcounts_plot_tall <- ggplot(data = wkcounts, aes(x = wk, y = intensity)) +
  geom_line(aes(color = type, linetype = type)) +
  scale_color_manual(values = c("#66c2a5", "#fc8d62")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  geom_segment(data = wkcount_summary, 
               aes(x = max_wk, xend = max_wk, y = 0, yend = max_value,
                   linetype = type),
               show.legend = FALSE,
               color = "gray30") +
  facet_wrap(~yr, ncol = 2) +
  labs(x = "Week", y = "Mean count") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        legend.margin = margin(c(0, 5, 0.5, 5)),
        legend.title = element_blank(),
        panel.grid = element_blank())
wkcounts_plot_tall
# ggsave("output/manuscript/redbud-weekly-mn-counts-5x2.png",
#        wkcounts_plot_tall,
#        width = 3.5, height = 9, units = "in", dpi = 600)

# Calculate estimated number of open flowers (by individual) ------------------#

# Remove observations where either intensity value is NA, then multiply
# (rounding values up to nearest integer)
flowers <- flowers %>%
  filter(!is.na(intensity_flowers) & !is.na(intensity_open)) %>%
  mutate(nopen = ceiling(intensity_flowers * intensity_open / 100))

# Remove observations past DOY 180
flowers <- filter(flowers, doy <= 180)

# Create rules for filtering plant-years. 
# Exclude plant-years that don't meet these criteria:
  # >=2 non-zero counts
  # At least one observation before DOY 100
  # Start and end observations with 0 counts
of_plantyrs <- flowers %>%
  arrange(id, obsdate) %>%
  group_by(id, yr) %>%
  summarize(nonzero_counts = sum(nopen > 0),
            min_doy = min(doy),
            start0 = ifelse(nopen[1] == 0, 1, 0),
            end0 = ifelse(last(nopen) == 0, 1, 0),
            .groups = "keep") %>%
  data.frame() %>%
  mutate(min_counts = ifelse(nonzero_counts > 1, 1, 0),
         earliest_obs = ifelse(min_doy <= 100, 1, 0),
         startend0 = ifelse(start0 == 1 & end0 == 1, 1, 0)) %>%
  mutate(remove = ifelse(min_counts + earliest_obs + startend0 == 3, 0, 1))

# Assess impact of filtering rules:
  count(of_plantyrs, remove, min_counts)
  # 944 plant-years had 0 or 1 non-zero count
  count(filter(of_plantyrs, min_counts == 1), remove, earliest_obs)
  # Of plant-years with >= 2 counts, 55 had earliest observation after DOY 100
  count(filter(of_plantyrs, min_counts == 1 & earliest_obs == 1), 
        remove, startend0)
  # Of plant-years with >= 2 counts with obs before DOY 100, 190 did not start
  # and end with 0 counts
  
  count(of_plantyrs, remove)
  sum(of_plantyrs$remove == 0) / nrow(of_plantyrs) * 100
  # Left with 528 plant-years (30.8% of 1717 plant years)

# Filter
of <- flowers %>%
  left_join(select(of_plantyrs, id, yr, remove), by = c("id", "yr")) %>%
  filter(remove == 0) %>%
  select(-remove)

# Merge with site information
of <- of %>%
  left_join(select(sites, site, state, lat, lon, elev), by = "site")

# Create dataframe with basic info for each site included in open flower dataset
ofsites <- of %>%
  group_by(site, lat, lon, elev, state) %>%
  summarize(n_plants = n_distinct(id),
            n_yrs = n_distinct(yr),
            n_plantyrs = n_distinct(paste0(id, "_", yr)),
            .groups = "keep") %>%
  data.frame()

# Write sites information to file to get PRISM weather data
# Downloading PRISM data through web interface, so we could get data at 800m 
# resolution without having to download daily rasters for CONUS, which takes a
# long time (https://prism.oregonstate.edu/explorer/bulk.php)
# write.table(select(ofsites, lat, lon, site), "weather-data/redbud-sites.csv",
#             sep = ",", row.names = FALSE, col.names = FALSE)

# Process weather data --------------------------------------------------------#

# Load daily data
prism_files <- list.files("weather-data/redbuds",
                          full.names = TRUE,
                          pattern = "stable|provisional")
for (i in 1:length(prism_files)) {
  prism1 <- read.csv(prism_files[i],
                     header = FALSE,
                     skip = 11,
                     col.names = c("site", "lon", "lat", "elev",
                                   "date", "ppt", "tmin", "tmean", "tmax"))
  if (i == 1) {
    weather <- prism1
  } else {
    weather <- rbind(weather, prism1)
  }
}
rm(prism1)

months.match <- data.frame(month.name = month.name, month = 1:12)

# Load 30-year (1991-2020) daily normals:
normsfile <- list.files("weather-data/redbuds/", 
                        full.names = TRUE,
                        pattern = "30yr")
norms <- read.csv(normsfile,
                  header = FALSE,
                  skip = 11,
                  col.names = c("site", "lon", "lat", "elev", "date",
                                "ppt30", "tmin30", "tmean30", "tmax30")) %>%
  separate_wider_delim(date, "-", names = c("month.name", "day")) %>%
  left_join(months.match, by = "month.name") %>%
  mutate(day = as.numeric(day)) %>%
  select(-month.name) %>%
  # Throw out Feb 29
  filter(!(month == 2 & day == 29))

# Summarize precipitation by month
ppt_month <- weather %>%
  select(site, date, ppt) %>%
  mutate(date = ymd(date),
         yr = year(date),
         month = month(date)) %>%
  group_by(site, yr, month) %>%
  summarize(ppt = sum(ppt),
            .groups = "keep") %>%
  data.frame()
# Add monthly normals
pptnorms_month <- norms %>%
  group_by(site, month) %>%
  summarize(ppt30 = sum(ppt30), .groups = "keep")
ppt_month <- ppt_month %>%
  left_join(pptnorms_month, by = c("site", "month"))

# Calculate winter precipitation (DJF) and % of 30-year norms
winter_ppt <- ppt_month %>%
  filter(month %in% c(12, 1:2)) %>%
  mutate(season = ifelse(month == 12, yr + 1, yr)) %>%
  filter(season %in% 2016:2025) %>%
  group_by(site, season) %>%
  summarize(winter_ppt = sum(ppt),
            winter_ppt30 = sum(ppt30),
            .groups = "keep") %>%
  mutate(winter_ppt_perc = (winter_ppt - winter_ppt30)/winter_ppt30 * 100) %>%
  data.frame()
  
# Summarize winter minimum temperatures (DJF) and calculate anomalies
winter_tmin <- weather %>%
  select(site, date, tmin) %>%
  mutate(date = ymd(date),
         yr = year(date),
         month = month(date)) %>%
  filter(month %in% c(12, 1:2)) %>%
  mutate(season = ifelse(month == 12, yr + 1, yr)) %>%
  filter(season %in% 2016:2025) %>%
  group_by(site, season) %>%
  summarize(winter_tmin = mean(tmin), .groups = "keep") %>%
  data.frame()
# Add normals
winter_tmin_norms <- norms %>%
  filter(month %in% c(12, 1:2)) %>%
  group_by(site) %>%
  summarize(winter_tmin30 = mean(tmin30)) %>%
  data.frame()
winter_tmin <- winter_tmin %>%
  left_join(winter_tmin_norms, by = "site") %>%
  mutate(winter_tmin_anom = winter_tmin - winter_tmin30)

# Calculate AGDD for DOY 1-90 (non-leap years = Jan-Mar) and calculate anomalies
agdd <- weather %>%
  select(site, date, tmean) %>%
  mutate(date = ymd(date),
         season = year(date),
         doy = yday(date)) %>%
  filter(season %in% 2016:2025) %>%
  filter(doy <= 90) %>%
  mutate(gdd = ifelse(tmean < 0, 0, tmean)) %>%
  group_by(site, season) %>%
  summarize(agdd = sum(gdd), .groups = "keep") %>%
  data.frame()
# Add normals
agdd_norms <- norms %>%
  mutate(doy = yday(ymd(paste0("2025-", month, "-", day)))) %>%
  filter(doy <= 90) %>%
  mutate(tmean30 = ifelse(tmean30 < 0, 0, tmean30)) %>%
  group_by(site) %>%
  summarize(agdd30 = sum(tmean30)) %>%
  data.frame()
agdd <- agdd %>%
  left_join(agdd_norms, by = "site") %>%
  mutate(agdd_anom = agdd - agdd30)

# Merge weather data
weather_vars <- winter_ppt %>%
  select(site, season, winter_ppt_perc) %>%
  left_join(select(winter_tmin, site, season, winter_tmin_anom),
            by = c("site", "season")) %>%
  left_join(select(agdd, site, season, agdd_anom), by = c("site" ,"season"))

# Create dataset to evaluate annual/spatial variation in max counts -----------# 

ofmax <- of %>%
  group_by(site, id, state, lat, lon, elev, yr) %>%
  summarize(nobs = n(),
            maxcount = max(nopen),
            .groups = "keep") %>%
  data.frame() %>%
  mutate(fyr = factor(yr)) %>%
  left_join(weather_vars, by = c("site" = "site", "yr" = "season"))

# Systematic differences among years?
ofmax %>%
  group_by(yr) %>%
  summarize(nplants = n_distinct(id),
            of_mean = round(mean(maxcount), 1),
            of_sd = round(sd(maxcount), 1),
            of_med = median(maxcount),
            of_min = min(maxcount),
            of_max = max(maxcount)) %>%
  data.frame()
# 2020 and 2022 had lower max counts, but they're not hugely different.

# Standardize variables
ofmax <- ofmax %>%
  mutate(lat_z = (lat - mean(lat))/sd(lat),
         lon_z = (lon - mean(lon))/sd(lon),
         elev_z = (elev - mean(elev))/sd(elev),
         winter_ppt_z = (winter_ppt_perc - mean(winter_ppt_perc))/sd(winter_ppt_perc),
         winter_tmin_z = (winter_tmin_anom - mean(winter_tmin_anom))/sd(winter_tmin_anom),
         agdd_z = (agdd_anom - mean(agdd_anom))/sd(agdd_anom))

# Correlations among potential covariates?
ofmax %>%
  select(winter_ppt_perc, winter_tmin_anom, agdd_anom, lat, lon, elev) %>%
  cor() %>%
  round(2)
# Nothing worrying here. Highest correlation is 0.55 between winter minimum
# temperatures and AGDD (which isn't too surprising since they both include
# temperature based metrics in Jan-Feb)

# Model annual/spatial variation in max counts --------------------------------# 

# Will want to log maximum counts since they vary over orders of magnitude
ofmax <- ofmax %>%
  mutate(maxcount_log = log(maxcount))

# Will use ML to compare models with different fixed effects, then run model
# using REML for inferences.

# Model that only includes topographic/geographic variables and random effects:
m_noweather <- lmer(maxcount_log ~ lat_z + lon_z + elev_z + (1|fyr) + (1|site),
                    data = ofmax, REML = FALSE)
summary(m_noweather)
# Max counts higher at higher latitudes. 
# No effect of elevation, and not strong evidence for longitude either
# Residual variance (individual) greater than site variance
# Random year effect tiny (which is interesting since there's no annual
# covariates in the model)

# Throw everything into a model:
m_full <- lmer(maxcount_log ~ lat_z + lon_z + elev_z + 
                 winter_ppt_z + winter_tmin_z + agdd_z + (1|fyr) + (1|site),
               data = ofmax, REML = FALSE)
summary(m_full)

# Remove variables that have no explanatory power (|t| < 1; AGDD and elev):
m_winter <- lmer(maxcount_log ~ lat_z + lon_z +
                   winter_ppt_z + winter_tmin_z + (1|fyr) + (1|site),
                 data = ofmax, REML = FALSE)
summary(m_winter)
anova(m_full, m_winter) # Reduced model is good (and has lower AIC)

# Keeping everything else in, even if not signficant at a 0.05 level. Refit with
# REML = TRUE
m_winter <- lmer(maxcount_log ~ lat_z + lon_z +
                   winter_ppt_z + winter_tmin_z + (1|fyr) + (1|site),
                 data = ofmax, REML = TRUE)
summary(m_winter)
confint(m_winter)

# More blooms on trees at higher latitudes, slighter more on trees further west
# More blooms following wetter winters (Dec-Feb) and slightly more after colder
# winters.

# Save table with model estimates to file:
max_fixed <- as.data.frame(summary(m_winter)$coefficients) %>%
  tibble::rownames_to_column("parameter") %>%
  select(-c("t value", "Std. Error")) %>%
  mutate(parameter_type = "fixed", .before = parameter)
max_random <- as.data.frame(VarCorr(m_winter)) %>%
  select(grp, sdcor) %>%
  rename(parameter = grp,
         Estimate = sdcor) %>%
  mutate(parameter_type = "random-SDs", .before = parameter)
max_ci <- as.data.frame(confint(m_winter)) %>%
  rename(lower = "2.5 %",
         upper = "97.5 %")
max_model_ests <- rbind(max_random, max_fixed) %>% cbind(max_ci)
write.csv(max_model_ests,
          file = "output/manuscript/max-model-estimates.csv",
          row.names = FALSE)

# Create dataset to evaluate variation in peak flower timing ------------------# 

# It doesn't make much sense to fit smooths/GAMs to the number of open flowers 
# for each plant-year, and there aren't any natural ways to group plants like 
# there were with gardens/accessions in the Rauschkolb preprint. Probably best
# to use date(s) of highest flower count. Given then, we'll want to filter out 
# plant-years where there were long gaps either before or after the date when
# the highest number of flowers were reported (to limit bias)

# Within each plant-year, calculate number of days since previous observation 
# (interval_before) and number of days until next observation (interval_after)
of <- of %>%
  arrange(id, obsdate) %>%
  mutate(obsdate = ymd(obsdate)) %>%
  mutate(interval_before = NA)
for (i in 2:nrow(of)) {
  of$interval_before[i] <- ifelse(
    of$id[i] == of$id[i - 1] & of$yr[i] == of$yr[i - 1],
    as.numeric(of$obsdate[i] - of$obsdate[i - 1]),
    NA
  )
}
nrows <- nrow(of)
of$interval_after <- c(of$interval_before[2:nrows], NA)

# For each plant-year, identify max flower count and date(s) when they occurred
of_intermediate <- of %>%
  group_by(yr, id) %>%
  mutate(maxvalue = max(nopen)) %>%
  ungroup() %>%
  filter(nopen == maxvalue) %>%
  group_by(yr, id, maxvalue) %>%
  summarize(first_max = min(doy),
            last_max = max(doy),
            .groups = "keep") %>%
  data.frame()

# Summarize information for each plant-year AND identify plant-years when there
# there were more than 21 or 30 days before or after dates when max flower
# counts were observed
of_plantyr <- of %>%
  left_join(of_intermediate, by = c("yr", "id")) %>%
  group_by(yr, id, site, lat, lon, elev, maxvalue, first_max, last_max) %>%
  summarize(nobs = n(),
            firstobs = min(doy),
            lastobs = max(doy), 
            nflower = sum(status_flowers),
            nopen = sum(status_open),
            interval_beforemax = interval_before[doy == first_max],
            interval_aftermax = interval_after[doy == last_max],
            .groups = "keep") %>%
  data.frame() %>%
  mutate(filter30 = ifelse(interval_beforemax > 30 | interval_aftermax > 30, 1, 0),
         filter21 = ifelse(interval_beforemax > 21 | interval_aftermax > 21, 1, 0))

# Look at filtering results by year
of_plantyr %>%
  group_by(yr) %>%
  summarize(nplants = n(),
            filter30 = n() - sum(filter30),
            filter21 = n() - sum(filter21))
# For now, will use 30-day interval filter

# Calculate peak date as mean date between first and last dates with max count
ofpeak <- of_plantyr %>%
  filter(filter30 == 0) %>%
  select(-c(filter30, filter21)) %>%
  mutate(max_span = last_max - first_max,
         peak = floor((first_max + last_max)/2),
         fyr = factor(yr),
         yr0 = yr - min(yr))
# There are a total of 509 plant-year combinations in the filtered dataset

# Attach weather variables
ofpeak <- ofpeak %>%
  left_join(weather_vars, by = c("site" = "site", "yr" = "season"))

# Standardize variables
ofpeak <- ofpeak %>%
  mutate(lat_z = (lat - mean(lat))/sd(lat),
         lon_z = (lon - mean(lon))/sd(lon),
         elev_z = (elev - mean(elev))/sd(elev),
         winter_ppt_z = (winter_ppt_perc - mean(winter_ppt_perc))/sd(winter_ppt_perc),
         winter_tmin_z = (winter_tmin_anom - mean(winter_tmin_anom))/sd(winter_tmin_anom),
         agdd_z = (agdd_anom - mean(agdd_anom))/sd(agdd_anom))

# Correlations among potential covariates?
ofpeak %>%
  select(winter_ppt_perc, winter_tmin_anom, agdd_anom, lat, lon, elev) %>%
  cor() %>%
  round(2)
# Again, nothing worrying here. Highest correlation of 0.55 between AGDD and 
# winter minimum temperatures

# Model variation in peak flower timing ---------------------------------------# 

# Model that only includes topographic/geographic variables and random effects:
mpeak_noweather <- lmer(peak ~ lat_z + lon_z + elev_z + (1|fyr) + (1|site),
                        data = ofpeak, REML = FALSE)
summary(mpeak_noweather)
confint(mpeak_noweather)
# Flowering peak later at higher latitudes, and to a much smaller degree,
# higher elevations. No effect of longitude.
# All random effects (year, site, individual [residual]) seem important

# Evidence of a trend over time?
mpeak_trend <- lmer(peak ~ lat_z + lon_z + elev_z + yr0 + (1|fyr) + (1|site),
                    data = ofpeak, REML = FALSE)
summary(mpeak_trend)
anova(mpeak_trend, mpeak_noweather)
# No strong evidence of a trend (though peak maybe slightly earlier over time?)

# Throw everything (but trend) into a model:
mpeak_full <- lmer(peak ~ lat_z + lon_z + elev_z +
                     winter_ppt_z + winter_tmin_z + agdd_z + (1|fyr) + (1|site),
                   data = ofpeak, REML = FALSE)
summary(mpeak_full)

# Remove variables that have no explanatory power (|t| < 1; lon, winter vars):
mpeak_agdd <- lmer(peak ~ lat_z + elev_z + agdd_z + (1|fyr) + (1|site),
                   data = ofpeak, REML = FALSE)
summary(mpeak_agdd)
anova(mpeak_full, mpeak_agdd) # Reduced model is good (and has lower AIC)

# Keeping everything else in, even if not signficant at a 0.05 level. Refit with
# REML = TRUE
mpeak_agdd <- lmer(peak ~ lat_z + elev_z + agdd_z + (1|fyr) + (1|site),
                   data = ofpeak, REML = TRUE)
summary(mpeak_agdd)
confint(mpeak_agdd)

# Later peak at higher latitudes and elevations. Earlier peak if Jan-Mar is 
# warmer than normal (2.6 days earlier for each 1-SD increase [80 degC] in AGDD)

# Save table with model estimates to file:
peak_fixed <- as.data.frame(summary(mpeak_agdd)$coefficients) %>%
  tibble::rownames_to_column("parameter") %>%
  select(-c("t value", "Std. Error")) %>%
  mutate(parameter_type = "fixed", .before = parameter)
peak_random <- as.data.frame(VarCorr(mpeak_agdd)) %>%
  select(grp, sdcor) %>%
  rename(parameter = grp,
         Estimate = sdcor) %>%
  mutate(parameter_type = "random-SDs", .before = parameter)
peak_ci <- as.data.frame(confint(mpeak_agdd)) %>%
  rename(lower = "2.5 %",
         upper = "97.5 %")
peak_model_ests <- rbind(peak_random, peak_fixed) %>% cbind(peak_ci)
write.csv(peak_model_ests,
          file = "output/manuscript/peak-model-estimates.csv",
          row.names = FALSE)
