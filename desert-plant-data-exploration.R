# Exploring intensity data for 5 desert species that are monitored at McDowell
# ER Zylstra

library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(ordinal)
library(brms)
library(leaflet)
library(geosphere)
library(mgcv)

# Plant phenophase classes
leaf_classes <- 1:5
flower_classes <- 6:9
fruit_classes <- 10:13

# McDowell site IDs
mcdo_sites <- c(24702, 24703, 24705)

# Load and format status-intensity data ---------------------------------------#

# List files with formatted intensity data
intensity_files <- list.files("npn-data",
                              pattern = "intensity-spp",
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

# si %>% count(interval) %>% mutate(prop = n / 212723) %>% round(3)

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

# Vast majority of observations are in AZ, so removing everything else
si <- filter(si, !is.na(state) & state == "AZ")

# Check latitude within Arizona?
# quantile(unique(si$latitude), probs = seq(0.1, 1, by = 0.1))
# Everything is below 34.53, and 90% are below 33.5, so ok as is

# Remove one site with mesquite that's confusing. 
# UTMs near Truth or Consequences, NM but state is listed as AZ. 
# Elevation is > 2000 m, much higher than everywhere else
si <- filter(si, site_id != 40329)

# Making observational filters a bit more liberal than was used for just 
# McDowell data

# Want to:
  # remove any pollen release data (qualitative data)
  # remove any spp-phenophase combos with < 10 plant-years
  # remove plant-php-year combos when fewer than 5 observations were made 
    # between dates associated with the 10th and 90th percentile AND fewer than 
    # 2 observations were made between dates associated with the 25th and 75th 
    # percentiles.

# Identify when species are typically in various phenophases
timings <- si %>%
  filter(phenophase_status == 1) %>%
  group_by(common_name, phenophase_description) %>%
  summarize(nyrs = n_distinct(yr),
            nplants = n_distinct(individual_id),
            nplantyrs = n_distinct(paste(yr, individual_id)),
            mean = round(mean(day_of_year)),
            q0.10 = quantile(day_of_year, probs = 0.10),
            q0.25 = quantile(day_of_year, probs = 0.25),
            q0.50 = quantile(day_of_year, probs = 0.50),
            q0.75 = quantile(day_of_year, probs = 0.75),
            q0.90 = quantile(day_of_year, probs = 0.90),
            .groups = "keep") %>%
  mutate(across(mean:q0.90, round)) %>%
  mutate(length_50 = q0.75 - q0.25,
         length_80 = q0.90 - q0.10) %>%
  data.frame()
timings

# Identify which spp-phenophases have < 10 plant-years
toremove <- timings %>%
  filter(nplantyrs < 10) %>%
  select(common_name, phenophase_description) %>%
  mutate(remove = 1)

# Filter based on phenophases
si <- si %>%
  filter(phenophase_description != "Pollen release (flowers)") %>%
  left_join(toremove, by = c("common_name", "phenophase_description")) %>%
  filter(is.na(remove)) %>%
  select(-remove)

# Summarize amount and quality of information for each plant, phenophase, year
pl_ph_yr <- si %>%
  left_join(select(timings, common_name, phenophase_description, 
                   q0.10, q0.25, q0.75, q0.90), 
            by = c("common_name", "phenophase_description")) %>%
  mutate(in50 = ifelse(day_of_year >= q0.25 & day_of_year <= q0.75, 1, 0)) %>%
  mutate(in80 = ifelse(day_of_year >= q0.10 & day_of_year <= q0.90, 1, 0)) %>%
  group_by(common_name, individual_id, site_id, phenophase_description, yr) %>%
  summarize(nobs = n(),
            first_obs = min(day_of_year),
            last_obs = max(day_of_year),
            nobs50 = sum(in50),
            nobs80 = sum(in80),
            mean_int = round(mean(interval, na.rm = TRUE), 2),
            max_int = ifelse(nobs ==1, NA, max(interval, na.rm = TRUE)),
            n_inphase = sum(phenophase_status),
            n_intvalue = sum(!is.na(intensity_value)),
            prop_inphase = round(n_inphase / nobs, 2),
            prop_intvalue = round(n_intvalue / n_inphase, 2),
            max_intensity = ifelse(sum(is.na(intensity_midpoint)) == n(),
                                   NA, max(intensity_midpoint, na.rm = TRUE)),
            .groups = "keep") %>%
  data.frame()

# Identify which plant-php-yr combinations to filter out
pl_ph_yr <- pl_ph_yr %>%
  mutate(remove = case_when(
    is.na(max_intensity) ~ 1,
    is.na(nobs50) ~ 1,
    is.na(nobs80) ~ 1,
    nobs80 < 5 ~ 1,
    nobs50 < 2 ~ 1,
    .default = 0
  ))

sif <- si %>%
  left_join(select(pl_ph_yr, individual_id, phenophase_description, yr, remove),
            by = c("individual_id", "phenophase_description", "yr")) %>%
  filter(remove == 0) %>%
  select(-remove)

# Extract information about sites ---------------------------------------------#

sites <- sif %>%
  distinct(latitude, longitude, site_id)
# write.table(sites, "weather-data/desert-plant-sites.csv", sep = ",", 
#             row.names = FALSE,
#             col.names = FALSE)

# Load weather data for these sites -------------------------------------------#
# Downloaded PRISM data through web interface, so we could get data at 800m 
# resolution without having to download daily rasters for CONUS, which takes a
# long time (https://prism.oregonstate.edu/explorer/bulk.php)

# Stuff commented out below was only run once, to create monthly summaries

prism_folder <- "weather-data/desert-plant-prism/"

# Daily data
# prism_files <- list.files(prism_folder, 
#                           full.names = TRUE,
#                           pattern = "stable")
# for (i in 1:length(prism_files)) {
#   prism1 <- read.csv(prism_files[i],
#                      header = FALSE,
#                      skip = 11,
#                      col.names = c("site_id", "lon", "lat", "elev",
#                                    "date", "ppt", "tmin", "tmax"))
#   if (i == 1) {
#     weather <- prism1
#   } else {
#     weather <- rbind(weather, prism1)
#   }
# }
# rm(prism1)
# 
# weather <- weather %>%
#   select(-c(lon, lat)) %>%
#   mutate(date = ymd(date),
#          yr = year(date),
#          month = month(date)) %>%
#   group_by(site_id, elev, yr, month) %>%
#   summarize(tmin_mn = mean(tmin),
#             tmax_mn = mean(tmax),
#             ppt = sum(ppt),
#             freezing = sum(tmin < 0),
#             .groups = "keep") %>%
#   data.frame()
#
# Write monthly data to file
# write.csv(weather, "weather-data/desert-plant-20112024.csv", row.names = FALSE)

# Read in monthly weather data
weather <- read.csv("weather-data/desert-plant-20112024.csv")

# 30-year monthly normals:
normsfile <- list.files(prism_folder, 
                        full.names = TRUE,
                        pattern = "30yr")
norms <- read.csv(normsfile,
                  header = FALSE,
                  skip = 11,
                  col.names = c("site_id", "lon", "lat", "elev",
                                "mon", "ppt30", "tmin30", "tmean30", 
                                "tmax30")) %>%
  filter(mon != "Annual") %>%
  mutate(month = as.integer(factor(mon, levels = month.name))) %>%
  select(-c(lon, lat, elev, mon))

# Summarize data by species ---------------------------------------------------#

# Species summaries
spp_summary <- sif %>%
  mutate(php_type = case_when(
    class_id %in% leaf_classes ~ "Leaves",
    class_id %in% flower_classes ~ "Flowers",
    class_id %in% fruit_classes ~ "Fruit"
  )) %>%
  group_by(common_name) %>%
  summarize(n_plants = n_distinct(individual_id),
            n_sites = n_distinct(site_name),
            first_yr = min(yr),
            last_yr = max(yr),
            n_yrs = n_distinct(yr),
            n_php = n_distinct(phenophase_id),
            leaf_php = n_distinct(phenophase_id[php_type == "Leaves"]),
            flower_php = n_distinct(phenophase_id[php_type == "Flowers"]),
            fruit_php = n_distinct(phenophase_id[php_type == "Fruit"]),
            .groups = "keep") %>%
  arrange(desc(n_sites), desc(n_plants), common_name) %>%
  data.frame()
spp_summary

# California barrel cactus only monitored at McDowell
count(filter(sif, common_name == "California barrel cactus"), site_id)
# Buck-horn cholla monitored at another 2 sites, but most obs at McDowell
count(filter(sif, common_name == "buck-horn cholla"), site_id)

# Summarize data by intensity category ----------------------------------------#

# Create table summarizing amount of information per intensity category
intensity_cats <- sif %>%
  group_by(class_id, intensity_label, intensity_type) %>%
  summarize(n = n(),
            n_spp = n_distinct(common_name),
            n_plants = n_distinct(individual_id),
            n_yrs = n_distinct(yr),
            values = paste(sort(unique(intensity_midpoint)), collapse = ", "),
            valuesq = paste(sort(unique(intensity_value)), collapse = ", "),
            .groups = "keep") %>%
  mutate(unique_values = ifelse(intensity_type == "qualitative", 
                                valuesq, values)) %>%
  select(-c(values, valuesq)) %>%
  data.frame()

# Create short name for intensity categories
intensity_cats <- intensity_cats %>%
  mutate(intensity_short = case_when(
    intensity_label == "No. young leaves" ~ "YoungLeaves",
    intensity_label == "No. breaking leaf buds" ~ "BreakingLeafBuds",
    intensity_label == "Leaf size (%)" ~ "LeafSize",
    intensity_label == "Leaf canopy fullness (%)" ~ "CanopyFullness",
    intensity_label == "Leaf canopy color (%)" ~ "CanopyColor",
    intensity_label == "No. flowers and flower buds" ~ "Flowers",
    intensity_label == "Open flowers (%)" ~ "OpenFlowers",
    intensity_label == "No. fruits" ~ "Fruits",
    intensity_label == "Ripe fruit (%)" ~ "RipeFruit",
    intensity_label == "No. fruit/seed drop" ~ "FruitDrop",
    .default = NA
  ))
intensity_cats

# Summarize filtered data for each plant, php, and year -----------------------#

# Recreating pl_ph_yr using rle(phenophase_status) to understand patterns:
  # length(rle$values) = number of 0/1 status sequences
  # first(rle$values) = state at first observation
  # last(rle$values) = state at last observation
pl_ph_yr <- sif %>%
  arrange(site_id, individual_id, phenophase_id, observation_date) %>%
  group_by(common_name, site_id, individual_id, class_id, 
           phenophase_description, intensity_label, yr) %>%
  summarize(nobs = n(),
            first_obs = min(day_of_year),
            last_obs = max(day_of_year),
            mean_int = round(mean(interval, na.rm = TRUE), 2),
            max_int = ifelse(nobs ==1, NA, max(interval, na.rm = TRUE)),
            n_inphase = sum(phenophase_status),
            prop_inphase = round(n_inphase / nobs, 2),
            n_states = length(rle(phenophase_status)$values),
            first_status = first(rle(phenophase_status)$values),
            last_status = last(rle(phenophase_status)$values),
            .groups = "keep") %>%
  data.frame()

# Summary across phenophases/intensity categories
php_summary <- pl_ph_yr %>%
  group_by(class_id, phenophase_description, intensity_label) %>%
  summarize(n_spp = n_distinct(common_name),
            n_plantyrs = n(),
            n_states_mn = round(mean(n_states), 2),
            n_states_max = max(n_states),
            prop_00 = round(sum(first_status == 0 & last_status == 0)/n_plantyrs, 2),
            .groups = "keep") %>%
  data.frame()
php_summary

# Quite a bit more fluctuation in phenophase status than there was in other 
# datasets (Kodiak, NEON), which is likely a function of the desert environment
# and how plants respond to weather changes. It's unclear the extent to which 
# this variation also reflects observer quality.

# Additional data filtering ---------------------------------------------------#

# Use a rule to exclude any species-phenophase combination where there's only 
# one year of data (if there are any)
sif <- sif %>%
  group_by(common_name, phenophase_description) %>%
  mutate(n_yrs = n_distinct(yr)) %>%
  ungroup() %>%
  filter(n_yrs > 1) %>%
  select(-n_yrs) %>%
  data.frame()

# Prep data for regression models ---------------------------------------------#

# Sort data and remove any observations with missing intensity values
sif_df <- sif %>%
  arrange(class_id, phenophase_description, common_name, site_id, 
          individual_id, observation_date) %>%
  filter(!is.na(intensity_midpoint)) %>%
  mutate(intensity_short = case_when(
    intensity_label == "No. breaking leaf buds" ~ "BreakingLeafBuds",
    intensity_label == "Leaf size (%)" ~ "LeafSize",
    intensity_label == "No. young leaves" ~ "YoungLeaves",
    intensity_label == "Leaf canopy fullness (%)" ~ "CanopyFullness",
    intensity_label == "Plant greeness (%)" ~ "Greenness",
    intensity_label == "Leaf canopy color (%)" ~ "CanopyColor",
    intensity_label == "No. flower heads" ~ "FlowerHeads",
    intensity_label == "No. flowers and flower buds" ~ "Flowers",
    intensity_label == "Open flowers (%)" ~ "OpenFlowers",
    intensity_label == "No. fruits" ~ "Fruits",
    intensity_label == "Ripe fruit (%)" ~ "RipeFruit",
    intensity_label == "No. fruit/seed drop" ~ "FruitDrop",
    .default = NA
  ))

# List of species-intensity category combinations
combos <- sif_df %>%
  mutate(plantyr = paste0(individual_id, "_", yr)) %>%
  group_by(class_id, intensity_label, intensity_type, common_name) %>%
  summarize(n_plantyrs = n_distinct(plantyr),
            n_plants = n_distinct(individual_id),
            n_yrs = n_distinct(yr),
            n_sites = n_distinct(site_name),
            .groups = "keep") %>%
  data.frame()

# What numeric intensity data are left?
combos_n <- filter(combos, intensity_type == "number")
combos_n
# Will focus on 3 intensity categories: flowers, fruit, fruit/seed drop and
# ignore leaf categories for now

# -----------------------------------------------------------------------------#
# Exploring variation in max counts for flowers phenophase --------------------#

flowers <- sif_df %>%
  filter(intensity_short == "Flowers") %>%
  select(common_name, individual_id, site_id, latitude, longitude, 
         phenophase_description, observation_date, yr, day_of_year,
         phenophase_status, intensity_label, intensity_midpoint) %>%
  rename(status = phenophase_status,
         phenophase = phenophase_description,
         id = individual_id,
         lat = latitude,
         lon = longitude,
         obsdate = observation_date, 
         doy = day_of_year)

# Distribution of intensity values
count(flowers, intensity_midpoint) %>%
  mutate(prop = round(n / sum(n), 3))
count(flowers, common_name, intensity_midpoint)

  # saguaro -------------------------------------------------------------------#
  # Aggregate data for each plant-year: saguaro
  flower_sag <- flowers %>%
    filter(common_name == "saguaro") %>%
    group_by(common_name, id, site_id, lat, lon, yr) %>%
    summarize(n_obs = n(),
              n_inphase = sum(status),
              max_count = max(intensity_midpoint),
              .groups = "keep") %>%
    data.frame()
  head(flower_sag)
  count(flower_sag, yr)
  # 2012-2015: 4-9 plants
  # 2016: 11 plants
  # 2017-2024: 25+ plants
  # Remove years with < 10 plants
  
  flower_sag <- flower_sag %>%
    group_by(yr) %>%
    mutate(n_plants = n()) %>%
    filter(n_plants >= 10) %>%
    select(-n_plants) %>%
    data.frame()
  
  count(flower_sag, max_count) # 324 total
  # 0:                        35
  # Less than 3 (1):           2
  # 3 to 10 (5):              17
  # 11 to 100 (50):          156
  # 101 to 1000 (500):       106
  # More than 1,000 (1001):    8
  
  # Create abundance categories
  # 0 (none), 1-5 (few), 50 (some), 500-1001 (many)
  # May need to combine none and few....
  flower_sag <- flower_sag %>%
    mutate(abund4 = case_when(
      max_count == 0 ~ "none",
      max_count %in% 1:5 ~ "few",
      max_count == 50 ~ "some",
      max_count >= 500 ~ "many",
    )) %>%
    mutate(abund4 = factor(abund4, 
                          levels = c("none", "few", "some", "many"),
                          ordered = TRUE)) %>%
    mutate(abund3 = case_when(
      max_count %in% 0:5 ~ "few",
      max_count == 50 ~ "some",
      max_count >= 500 ~ "many",
    )) %>%
    mutate(abund3 = factor(abund3, 
                          levels = c("few", "some", "many"),
                          ordered = TRUE))
  
  # Convert variables to factor
  flower_sag$fyr <- factor(flower_sag$yr)
  flower_sag$site <- factor(flower_sag$site_id)
  flower_sag$id <- factor(flower_sag$id)
  flower_sag$mcdo <- factor(ifelse(flower_sag$site_id %in% mcdo_sites, 1, 0))
  
  # Visualize for each year
  ggplot(flower_sag, aes(y = abund4, x = fyr)) +
    geom_jitter(aes(color = mcdo), width = 0.20, height = 0.20) +
    labs(x = "", y = "Abundance category (4)", title = "saguaro, flowers")
  ggplot(flower_sag, aes(y = abund3, x = fyr)) +
    geom_jitter(aes(color = mcdo), width = 0.20, height = 0.20) +
    labs(x = "", y = "Abundance category (3)", title = "saguaro, flowers")
  
  # Format and merge with weather data
    # Min temps in winter (DJF)
    # Freezing days in winter/spring (DJFMAM)
    # Precip in fall (SON)
    # Precip, 6 mo (fall + winter)
    # Precip, 9 mo (summer + fall + winter)
    # Max temps in summer (JJA)
    # Max temps in fall (SON)
    # Max temps in winter (DJF)  
    # Max temps in spring (MAM)
    # Max temps in summer + fall
    # Max temps, 12 mo
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
    group_by(site_id, elev, seasonyr) %>%
    summarize(tmin_wi = mean(tmin[season == "wi"]),
              freeze = sum(freezing),
              ppt_fa = sum(ppt[season == "fa"]),
              ppt_6 = sum(ppt[season %in% c("fa", "wi")]),
              ppt_9 = sum(ppt[season != "sp"]),
              tmax_su = mean(tmax[season == "su"]),
              tmax_fa = mean(tmax[season == "fa"]),
              tmax_wi = mean(tmax[season == "wi"]),
              tmax_sp = mean(tmax[season == "sp"]),
              tmax_6 = mean(tmax[season %in% c("su", "fa")]),
              tmax_12 = mean(tmax),
              .groups = "keep") %>%
    filter(seasonyr %in% unique(flower_sag$yr)) %>%
    data.frame()
  
  # Want to calculate anomalies, so, first calculate 30-year normals
  normvars <- norms %>%
    group_by(site_id) %>%
    summarize(tmin_wi_30 = mean(tmin30[month %in% c(12, 1:2)]),
              ppt_fa_30 = sum(ppt30[month %in% 9:11]),
              ppt_6_30 = sum(ppt30[month %in% c(9:12, 1:2)]),
              ppt_9_30 = sum(ppt30[month %in% c(6:12, 1:2)]),
              tmax_su_30 = mean(tmax30[month %in% 6:8]),
              tmax_fa_30 = mean(tmax30[month %in% 9:11]),
              tmax_wi_30 = mean(tmax30[month %in% c(12, 1:2)]),
              tmax_sp_30 = mean(tmax30[month %in% 3:5]),
              tmax_6_30 = mean(tmax30[month %in% 6:11]),
              tmax_12_30 = mean(tmax30),
              # Overall mean annual temperature, precipitation
              tmean_30 = mean(tmean30),
              ppt_30 = sum(ppt30)) %>%
    data.frame()
  # Merge normals with annual weather and calculate anomalies
  weathervars <- weathervars %>%
    left_join(normvars, by = "site_id") %>%
    mutate(tmin_wi_a = tmin_wi - tmin_wi_30,
           ppt_fa_a = ppt_fa - ppt_fa_30,
           ppt_6_a = ppt_6 - ppt_6_30,
           ppt_9_a = ppt_9 - ppt_9_30,
           tmax_su_a = tmax_su - tmax_su_30,
           tmax_fa_a = tmax_fa - tmax_fa_30,
           tmax_wi_a = tmax_wi - tmax_wi_30,
           tmax_sp_a = tmax_sp - tmax_sp_30,
           tmax_6_a = tmax_6 - tmax_6_30,
           tmax_12_a = tmax_12 - tmax_12_30)
  
  # Where are the monitored saguaros?
  sag_locs <- flower_sag %>%
    left_join(distinct(weathervars, site_id, elev), by = "site_id") %>%
    group_by(site_id, lat, lon, elev) %>%
    summarize(n_plants = n_distinct(id),
              n_yrs = n_distinct(yr),
              .groups = "keep") %>%
    left_join(normvars, by = "site_id") %>%
    data.frame()
  
  # Geographic location
  # ggplot(sag_locs, aes(x = lon, y = lat)) +
  #   geom_point(aes(size = n_plants)) +
  #   scale_size_continuous(breaks = c(1, 5, 10), range = c(2, 5))
  leaflet(sag_locs) %>% addTiles() %>%
    addCircleMarkers(~lon, ~lat, radius = ~n_plants, fillOpacity = 0.6) %>%
    addCircleMarkers(lng = -110.640, lat = 32.071, fillOpacity = 0.6, 
                     color = "red", radius = 5)
  
  # Elevation, Climate norms
  ggplot(sag_locs, aes(x = elev)) +
    geom_histogram(bins = 30, fill = "gray") +
    labs(x = "Elevation (m)", y = "No. sites") +
    geom_vline(aes(xintercept = mean(elev)), color = "steelblue4",
               linetype = "dashed") +
    geom_vline(aes(xintercept = 1074), color = "salmon4", linetype = "dashed") +
    annotate(geom = "text", label = "Mean", color = "steelblue4", 
             x = mean(sag_locs$elev) - 5, y = 9, hjust = 1) +
    annotate(geom = "text", label = "Renzi", color = "salmon4", 
             x = 1074 - 5, y = 9, hjust = 1) +
    theme_bw()
  ggplot(sag_locs) +
    geom_point(aes(x = tmean_30, y = ppt_30, color = elev),
               size = 2) +
    scale_color_viridis_c() +
    labs(x = "Mean annual temperature (degC)", 
         y = "Mean annual precipitation (mm)",
         color = "Elevation (m)") +
    theme_bw()
  # One much higher elevation (1195), colder, wetter than rest (SW of Tucson)
  # One of the lowest elevation sites (397) is much drier than expected (W of Phoenix)
    
  # Add weather data to flowering dataset (and standardize)
  flower_sag <- flower_sag %>%
    left_join(weathervars, by = c("site_id" = "site_id", "yr" = "seasonyr")) %>%
    mutate(tmin_wi_z = (tmin_wi_a - mean(tmin_wi_a)) / sd(tmin_wi_a),
           freeze_z = (freeze - mean(freeze)) / sd(freeze),
           ppt_fa_z = (ppt_fa_a - mean(ppt_fa_a)) / sd(ppt_fa_a),
           ppt_6_z = (ppt_6_a - mean(ppt_6_a)) / sd(ppt_6_a),
           ppt_9_z = (ppt_9_a - mean(ppt_9_a)) / sd(ppt_9_a),
           tmax_su_z = (tmax_su_a - mean(tmax_su_a)) / sd(tmax_su_a),
           tmax_fa_z = (tmax_fa_a - mean(tmax_fa_a)) / sd(tmax_fa_a),
           tmax_wi_z = (tmax_wi_a - mean(tmax_wi_a)) / sd(tmax_wi_a),
           tmax_sp_z = (tmax_sp_a - mean(tmax_sp_a)) / sd(tmax_sp_a),
           tmax_6_z = (tmax_6_a - mean(tmax_6_a)) / sd(tmax_6_a),
           tmax_12_z = (tmax_12_a - mean(tmax_12_a)) / sd(tmax_12_a),
           ppt_norm_z = (ppt_30 - mean(ppt_30)) / sd(ppt_30),
           tmean_norm_z = (tmean_30 - mean(tmean_30)) / sd(tmean_30))

  # Run ML ordinal models (4 cats) with plant random effects
  # Year with plant RE
  m_yr1 <- clmm(abund4 ~ fyr + (1|id), Hess = TRUE, nAGQ = 10,
               data = flower_sag)
  # Year with plant and site RE (can't use Quadrature methods with > 1 RE)
  m_yr2 <- clmm(abund4 ~ fyr + (1|id) + (1|site), Hess = TRUE,
                data = flower_sag)
  summary(m_yr1)
  summary(m_yr2)
  anova(m_yr1, m_yr2) 
  # Site RE is worth including (at least with few fixed effects)
  # Year effects similar in both models
  
  # Does it change things if we put 0 counts in the "few" category?
  m3_yr2 <- clmm(abund3 ~ fyr + (1|id) + (1|site), Hess = TRUE,
                data = flower_sag)
  summary(m3_yr2)
  # No, so we'll just stick with 4 categories for time being
  
  # Try different precip models
  m_ppt_fa <- clmm(abund4 ~ ppt_fa_z + (1|id) + (1|site), Hess = TRUE,
                   data = flower_sag)
  m_ppt_6 <- clmm(abund4 ~ ppt_6_z + (1|id) + (1|site), Hess = TRUE,
                  data = flower_sag)
  m_ppt_9 <- clmm(abund4 ~ ppt_9_z + (1|id) + (1|site), Hess = TRUE,
                  data = flower_sag)
  AIC(m_ppt_fa, m_ppt_6, m_ppt_9, m_yr2)
  summary(m_ppt_9)
  # Negative effect of precipitation (similar to Renzi)
  # 9-month precip much better than 3- or 6-mo, but an independent year model
  # is way better than all precip models
  
  # Does the effect of precipitation vary with mean annual rainfall?
  m_ppt_9b <- clmm(abund4 ~ ppt_9_z * ppt_norm_z + (1|id) + (1|site), 
                   Hess = TRUE, data = flower_sag)
  summary(m_ppt_9b)
  anova(m_ppt_9, m_ppt_9b) # Adding mean precip didn't help

  # Winter minimum temperatures
  m_tmin <- clmm(abund4 ~ tmin_wi_z + (1|id) + (1|site), Hess = TRUE,
                 data = flower_sag)
  summary(m_tmin) 
  AIC(m_ppt_9, m_tmin, m_yr2)
  # Negative effect of winter temperature (not similar to Renzi)
  # More support for a precipitation model than winter temperature, but 
  # year model still better
  
  # Does the effect of winter minimum temperatures vary with mean temperatures?
  # Having trouble with this with both REs and without nAGQ
  # m_tminb <- clmm(abund4 ~ tmin_wi_z * tmean_norm_z + (1|id) + (1|site), 
  #                 Hess = TRUE, data = flower_sag)
  m_tmin <- clmm(abund4 ~ tmin_wi_z + (1|id), nAGQ = 10, 
                 Hess = TRUE, data = flower_sag)
  m_tmini <- clmm(abund4 ~ tmin_wi_z * tmean_norm_z + (1|id), nAGQ = 10, 
                  Hess = TRUE, data = flower_sag)
  anova(m_tmin, m_tmini)
  # A model with an interaction is not a lot better (P = 0.16)
  summary(m_tmini)
  
  # Re-run model with clmm2 (for predict)
  m_tmini2 <- clmm2(abund4 ~ tmin_wi_z * tmean_norm_z, random = id, nAGQ = 10, 
                    Hess = TRUE, data = flower_sag)
  # Create dataframe for prediction
  newdat <- expand.grid(
    abund4 = unique(flower_sag$abund4),
    tmean_norm_z = c(min(flower_sag$tmean_norm_z), 0, 
                     max(flower_sag$tmean_norm_z)),
    tmin_wi_z = seq(min(flower_sag$tmin_wi_z), max(flower_sag$tmin_wi_z), 
                    length = 100),
    KEEP.OUT.ATTRS = FALSE
  )
  # Plot predictions (= probability that saguaro will have X flowers
  # given winter min temperatures for a cold/average/warm site)
  preds <- cbind(newdat, est = predict(m_tmini2, newdata = newdat)) %>%
    arrange(abund4, tmean_norm_z, tmin_wi_z) %>%
    mutate(mean_temp = case_when(
      tmean_norm_z == min(flower_sag$tmean_norm_z) ~ "cold",
      tmean_norm_z == 0 ~ "average",
      tmean_norm_z == max(flower_sag$tmean_norm_z) ~ "warm",
    )) %>%
    mutate(mean_temp = factor(mean_temp, levels = c("cold", "average", "warm")))
  ggplot(preds, aes(x = tmin_wi_z, y = est)) +
    geom_line(aes(color = abund4), linewidth = 1.3) +
    scale_color_brewer(palette = "BrBG") +
    facet_grid(mean_temp ~ .) +
    labs(x = "Winter min temp anomaly (z)", y = "Probability", color = "Flowers") +
    theme_bw()
  # Results are interesting: Cold site is very different than average/warm.
  # At average/warm site: 
  # Low probability of no/few flowers all the time though slightly higher 
  # probability in a warm winter. High probability of many flowers in a cold
  # winter and lower probability of many flowers in a warm year. Probability of
  # 11-100 (some flowers) always higher in warmer winter.
  # At cold site:
  # Low probability of many flowers all the time though slightly higher in warm
  # winter. Probability of no/few flowers higher in cold winter and lower in 
  # warm winter. Probability of some flowers increases with winter temps. 
  
  # Combine 9-mo precipitation and winter temperatures (additive)
  m_comb <- clmm(abund4 ~ ppt_9_z + tmin_wi_z * tmean_norm_z + (1|id), 
                 nAGQ = 10, Hess = TRUE, data = flower_sag)
  summary(m_comb)
  anova(m_comb, m_tmini)
  anova(m_comb, m_ppt_9)
  
  # Combine 9-mo precipitation and winter temperatures (3-way interaction)
  m_combi <- clmm(abund4 ~ ppt_9_z * tmin_wi_z * tmean_norm_z + (1|id), 
                  nAGQ = 10, Hess = TRUE, data = flower_sag)
  summary(m_combi)
  anova(m_comb, m_combi)
  anova(m_ppt_9, m_combi)  
  AIC(m_yr2, m_ppt_9, m_tmin, m_tmini, m_comb, m_combi)
  # Best is year model, then 3-way interaction, then ppt
  
  # Re-run model with clmm2 (for predict)
  m_combi2 <- clmm2(abund4 ~ ppt_9_z * tmin_wi_z * tmean_norm_z, random = id, 
                    nAGQ = 10, Hess = TRUE, data = flower_sag)
  # Create dataframe for prediction
  newdat_combi <- expand.grid(
    abund4 = unique(flower_sag$abund4),
    tmean_norm_z = c(min(flower_sag$tmean_norm_z), 0, 
                     max(flower_sag$tmean_norm_z)),
    tmin_wi_z = c(min(flower_sag$tmin_wi_z), 0, 
                     max(flower_sag$tmin_wi_z)),
    ppt_9_z = seq(min(flower_sag$ppt_9_z), max(flower_sag$ppt_9_z), 
                  length = 100),
    KEEP.OUT.ATTRS = FALSE
  )
  # Plot predictions (= probability that saguaro will have X flowers
  # given 9-month precipitation for cold/average/warm winter min temps at a 
  # cold/average/warm site)
  preds_combi <- cbind(newdat_combi, 
                       est = predict(m_combi2, newdata = newdat_combi)) %>%
    arrange(abund4, tmean_norm_z, tmin_wi_z, ppt_9_z) %>%
    mutate(mean_temp = case_when(
      tmean_norm_z == min(flower_sag$tmean_norm_z) ~ "cold",
      tmean_norm_z == 0 ~ "average",
      tmean_norm_z == max(flower_sag$tmean_norm_z) ~ "warm",
    )) %>%
    mutate(winter = case_when(
      tmin_wi_z == min(flower_sag$tmin_wi_z) ~ "cold winter",
      tmin_wi_z == 0 ~ "average winter",
      tmin_wi_z == max(flower_sag$tmin_wi_z) ~ "warm winter"
    )) %>%
    mutate(mean_temp = factor(mean_temp, 
                              levels = c("cold", "average", "warm"))) %>%
    mutate(winter = factor(winter, 
                           levels = c("cold winter", "average winter", 
                                      "warm winter")))
  ggplot(preds_combi, aes(x = ppt_9_z, y = est)) +
    geom_line(aes(color = abund4), linewidth = 1.3) +
    scale_color_brewer(palette = "BrBG") +
    facet_grid(mean_temp ~ winter) +
    labs(x = "Cumulative precipitation, 9 mo (z)", 
         y = "Probability", 
         color = "Flowers") +
    theme_bw()

  # Does it change things if we replace mean annual temps with elevation?
  flower_sag <- flower_sag %>%
    mutate(elev_z = (elev - mean(elev)) / sd(elev))
  m_tminielev <- clmm(abund4 ~ tmin_wi_z * elev_z + (1|id), nAGQ = 10,
                   Hess = TRUE, data =flower_sag)
  summary(m_tminielev) # Signif interaction
  m_combelev <- clmm(abund4 ~ ppt_9_z + tmin_wi_z * elev_z + (1|id), 
                   nAGQ = 10, Hess = TRUE, data = flower_sag)
  summary(m_combelev) # Precip helps, but now interaction P = 0.12
  m_combielev <- clmm(abund4 ~ ppt_9_z * tmin_wi_z * elev_z + (1|id), 
                      nAGQ = 10, Hess = TRUE, data = flower_sag)
  summary(m_combielev)
  AIC(m_yr2, m_ppt_9, m_tmin, m_tminielev, m_combelev, m_combielev)
  # Same as using mean temps:
  # yr best, followed by 3-way interaction, then ppt
  
  m_combielev2 <- clmm2(abund4 ~ ppt_9_z * tmin_wi_z * elev_z, random = id, 
                        nAGQ = 10, Hess = TRUE, data = flower_sag)

  # Create dataframe for prediction
  newdat_combie <- expand.grid(
    abund4 = unique(flower_sag$abund4),
    elev_z = c(min(flower_sag$elev_z), 0, max(flower_sag$elev_z)),
    tmin_wi_z = c(min(flower_sag$tmin_wi_z), 0, 
                  max(flower_sag$tmin_wi_z)),
    ppt_9_z = seq(min(flower_sag$ppt_9_z), max(flower_sag$ppt_9_z), 
                  length = 100),
    KEEP.OUT.ATTRS = FALSE
  )
  # Plot predictions (= probability that saguaro will have X flowers
  # given 9-month precipitation for cold/average/warm winter min temps at a 
  # cold/average/warm site)
  preds_combie <- cbind(newdat_combie, 
                       est = predict(m_combielev2, newdata = newdat_combie)) %>%
    arrange(abund4, elev_z, tmin_wi_z, ppt_9_z) %>%
    mutate(loc = case_when(
      elev_z == min(flower_sag$elev_z) ~ "low elev",
      elev_z == 0 ~ "average elev",
      elev_z == max(flower_sag$elev_z) ~ "high elev",
    )) %>%
    mutate(winter = case_when(
      tmin_wi_z == min(flower_sag$tmin_wi_z) ~ "cold winter",
      tmin_wi_z == 0 ~ "average winter",
      tmin_wi_z == max(flower_sag$tmin_wi_z) ~ "warm winter"
    )) %>%
    mutate(loc = factor(loc, 
                        levels = c("high elev", "average elev", "low elev"))) %>%
    mutate(winter = factor(winter, 
                           levels = c("cold winter", "average winter", 
                                      "warm winter")))
  ggplot(preds_combie, aes(x = ppt_9_z, y = est)) +
    geom_line(aes(color = abund4), linewidth = 1.3) +
    scale_color_brewer(palette = "BrBG") +
    facet_grid(loc ~ winter) +
    labs(x = "Cumulative precipitation, 9 mo (z)", 
         y = "Probability", 
         color = "Flowers") +
    theme_bw()
  
  # Does it matter if we throw out the observations of one saguaro that's way
  # higher elevation and was only observed in 2 years?
  flower_sag2 <- flower_sag %>%
    filter(elev < 1100) %>%
    # Re-calculate z-scores
    mutate(tmin_wi_z = (tmin_wi_a - mean(tmin_wi_a)) / sd(tmin_wi_a),
           freeze_z = (freeze - mean(freeze)) / sd(freeze),
           ppt_fa_z = (ppt_fa_a - mean(ppt_fa_a)) / sd(ppt_fa_a),
           ppt_6_z = (ppt_6_a - mean(ppt_6_a)) / sd(ppt_6_a),
           ppt_9_z = (ppt_9_a - mean(ppt_9_a)) / sd(ppt_9_a),
           tmax_su_z = (tmax_su_a - mean(tmax_su_a)) / sd(tmax_su_a),
           tmax_fa_z = (tmax_fa_a - mean(tmax_fa_a)) / sd(tmax_fa_a),
           tmax_wi_z = (tmax_wi_a - mean(tmax_wi_a)) / sd(tmax_wi_a),
           tmax_sp_z = (tmax_sp_a - mean(tmax_sp_a)) / sd(tmax_sp_a),
           tmax_6_z = (tmax_6_a - mean(tmax_6_a)) / sd(tmax_6_a),
           tmax_12_z = (tmax_12_a - mean(tmax_12_a)) / sd(tmax_12_a),
           ppt_norm_z = (ppt_30 - mean(ppt_30)) / sd(ppt_30),
           tmean_norm_z = (tmean_30 - mean(tmean_30)) / sd(tmean_30),
           elev_z = (elev - mean(elev)) / sd(elev))

  m2_yr <- clmm(abund4 ~ fyr + (1|id), nAGQ = 10, 
                Hess = TRUE, data = flower_sag2)
  m2_tmin <- clmm(abund4 ~ tmin_wi_z + (1|id), nAGQ = 10, 
                  Hess = TRUE, data = flower_sag2)
  m2_tmini <- clmm(abund4 ~ tmin_wi_z * tmean_norm_z + (1|id), nAGQ = 10, 
                   Hess = TRUE, data = flower_sag2)
  m2_tminielev <- clmm(abund4 ~ tmin_wi_z * elev_z + (1|id), nAGQ = 10, 
                       Hess = TRUE, data = flower_sag2)
  m2_ppt <- clmm(abund4 ~ ppt_9_z + (1|id), nAGQ = 10, 
                 Hess = TRUE, data = flower_sag2)
  m2_comb <- clmm(abund4 ~ ppt_9_z + tmin_wi_z * tmean_norm_z + (1|id), 
                  nAGQ = 10, Hess = TRUE, data = flower_sag2)
  m2_combelev <- clmm(abund4 ~ ppt_9_z + tmin_wi_z * elev_z + (1|id), 
                      nAGQ = 10, Hess = TRUE, data = flower_sag2)
  m2_combi <- clmm(abund4 ~ ppt_9_z * tmin_wi_z * tmean_norm_z + (1|id), 
                   nAGQ = 10, Hess = TRUE, data = flower_sag2)
  m2_combielev <- clmm(abund4 ~ ppt_9_z * tmin_wi_z * elev_z + (1|id), nAGQ = 10, 
                       Hess = TRUE, data = flower_sag2)
  AIC(m2_yr, m2_ppt, m2_tmin, m2_tmini, m2_comb, m2_combi)
  AIC(m2_yr, m2_ppt, m2_tmin, m2_tminielev, m2_combelev, m2_combielev)
  
  m2_combielev <- clmm2(abund4 ~ ppt_9_z * tmin_wi_z * elev_z, random = id, 
                        nAGQ = 10, Hess = TRUE, data = flower_sag2)
  
  # Create dataframe for prediction
  newdat_combie2 <- expand.grid(
    abund4 = unique(flower_sag2$abund4),
    elev_z = c(min(flower_sag2$elev_z), 0, max(flower_sag2$elev_z)),
    tmin_wi_z = c(min(flower_sag2$tmin_wi_z), 0, 
                  max(flower_sag2$tmin_wi_z)),
    ppt_9_z = seq(min(flower_sag2$ppt_9_z), max(flower_sag2$ppt_9_z), 
                  length = 100),
    KEEP.OUT.ATTRS = FALSE
  )
  # Plot predictions (= probability that saguaro will have X flowers
  # given 9-month precipitation for cold/average/warm winter min temps at a 
  # cold/average/warm site)
  preds_combie2 <- cbind(newdat_combie2, 
                        est = predict(m2_combielev, newdata = newdat_combie2)) %>%
    arrange(abund4, elev_z, tmin_wi_z, ppt_9_z) %>%
    mutate(loc = case_when(
      elev_z == min(flower_sag2$elev_z) ~ "low elev",
      elev_z == 0 ~ "average elev",
      elev_z == max(flower_sag2$elev_z) ~ "high elev",
    )) %>%
    mutate(winter = case_when(
      tmin_wi_z == min(flower_sag2$tmin_wi_z) ~ "cold winter",
      tmin_wi_z == 0 ~ "average winter",
      tmin_wi_z == max(flower_sag2$tmin_wi_z) ~ "warm winter"
    )) %>%
    mutate(loc = factor(loc, 
                        levels = c("high elev", "average elev", "low elev"))) %>%
    mutate(winter = factor(winter, 
                           levels = c("cold winter", "average winter", 
                                      "warm winter")))
  pp <- ggplot(preds_combie2, aes(x = ppt_9_z, y = est)) +
    geom_line(aes(color = abund4), linewidth = 1.3) +
    scale_color_brewer(palette = "BrBG", 
                       labels = c("0", "1-10", "11-100", '>100')) +
    facet_grid(loc ~ winter) +
    labs(x = "Cumulative 9-month precipitation, anomaly (standardized)", 
         y = "Probability", 
         color = "Flowers",
         title = "Saguaro flowers, excluding highest site (all sites < 1100 m)") +
    theme_bw() +
    theme(legend.position = "bottom")
  # ggsave("output/saguaro-flower-preds-elev.png",
  #        pp, height = 6.5, width = 6.5, units = "in")

  # velvet mesquite -----------------------------------------------------------#
  # Aggregate data for each plant-year: velvet mesquite
  flower_mesq <- flowers %>%
    filter(common_name == "velvet mesquite") %>%
    group_by(common_name, id, site_id, lat, lon, yr) %>%
    summarize(n_obs = n(),
              n_inphase = sum(status),
              max_count = max(intensity_midpoint),
              .groups = "keep") %>%
    data.frame()
  head(flower_mesq)
  count(flower_mesq, yr)
  # 2012-2018: 3, 8, 6, 1, 6, 11, 8
  # Remove years with < 8 plants?
  
  flower_mesq <- flower_mesq %>%
    group_by(yr) %>%
    mutate(n_plants = n()) %>%
    filter(n_plants >= 8) %>%
    select(-n_plants) %>%
    data.frame()
  
  count(flower_mesq, max_count) # 170 total
  # 0:                        3
  # Less than 3 (1):          1
  # 3 to 10 (5):              1
  # 11 to 100 (50):          20
  # 101 to 1000 (500):       51
  # 1001 to 10000 (5000):    78
  # More than 10000 (10001): 16
  
  # Create abundance categories
  # 0-50 (few/none), 500 (some), >= 5000 (many)
  flower_mesq <- flower_mesq %>%
    mutate(abund3 = case_when(
      max_count %in% 0:50 ~ "few",
      max_count == 500 ~ "some",
      max_count >= 5000 ~ "many"
    )) %>%
    mutate(abund3 = factor(abund3, 
                           levels = c("few", "some", "many"),
                           ordered = TRUE))
  
  # Convert variables to factor
  flower_mesq$fyr <- factor(flower_mesq$yr)
  flower_mesq$site <- factor(flower_mesq$site_id)
  flower_mesq$id <- factor(flower_mesq$id)
  flower_mesq$mcdo <- factor(ifelse(flower_mesq$site_id %in% mcdo_sites, 1, 0))
  
  # Visualize for each year
  ggplot(flower_mesq, aes(y = abund3, x = fyr)) +
    geom_jitter(aes(color = mcdo), width = 0.20, height = 0.20) +
    labs(x = "", y = "Abundance category (3)", title = "velvet mesquite, flowers")
  
  # Format and merge with weather data
  # Min temps in winter (DJF)
  # Freezing days in winter/spring (DJFMAM)
  # Precip in fall (SON)
  # Precip, 6 mo (fall + winter)
  # Precip, 9 mo (summer + fall + winter)
  # Max temps in summer (JJA)
  # Max temps in fall (SON)
  # Max temps in winter (DJF)  
  # Max temps in spring (MAM)
  # Max temps in summer + fall
  # Max temps, 12 mo
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
    group_by(site_id, elev, seasonyr) %>%
    summarize(tmin_wi = mean(tmin[season == "wi"]),
              freeze = sum(freezing),
              ppt_fa = sum(ppt[season == "fa"]),
              ppt_6 = sum(ppt[season %in% c("fa", "wi")]),
              ppt_9 = sum(ppt[season != "sp"]),
              tmax_su = mean(tmax[season == "su"]),
              tmax_fa = mean(tmax[season == "fa"]),
              tmax_wi = mean(tmax[season == "wi"]),
              tmax_sp = mean(tmax[season == "sp"]),
              tmax_6 = mean(tmax[season %in% c("su", "fa")]),
              tmax_12 = mean(tmax),
              .groups = "keep") %>%
    filter(seasonyr %in% unique(flower_mesq$yr)) %>%
    data.frame()
  
  # Want to calculate anomalies, so, first calculate 30-year normals
  normvars <- norms %>%
    group_by(site_id) %>%
    summarize(tmin_wi_30 = mean(tmin30[month %in% c(12, 1:2)]),
              ppt_fa_30 = sum(ppt30[month %in% 9:11]),
              ppt_6_30 = sum(ppt30[month %in% c(9:12, 1:2)]),
              ppt_9_30 = sum(ppt30[month %in% c(6:12, 1:2)]),
              tmax_su_30 = mean(tmax30[month %in% 6:8]),
              tmax_fa_30 = mean(tmax30[month %in% 9:11]),
              tmax_wi_30 = mean(tmax30[month %in% c(12, 1:2)]),
              tmax_sp_30 = mean(tmax30[month %in% 3:5]),
              tmax_6_30 = mean(tmax30[month %in% 6:11]),
              tmax_12_30 = mean(tmax30),
              # Overall mean annual temperature, precipitation
              tmean_30 = mean(tmean30),
              ppt_30 = sum(ppt30)) %>%
    data.frame()
  # Merge normals with annual weather and calculate anomalies
  weathervars <- weathervars %>%
    left_join(normvars, by = "site_id") %>%
    mutate(tmin_wi_a = tmin_wi - tmin_wi_30,
           ppt_fa_a = ppt_fa - ppt_fa_30,
           ppt_6_a = ppt_6 - ppt_6_30,
           ppt_9_a = ppt_9 - ppt_9_30,
           tmax_su_a = tmax_su - tmax_su_30,
           tmax_fa_a = tmax_fa - tmax_fa_30,
           tmax_wi_a = tmax_wi - tmax_wi_30,
           tmax_sp_a = tmax_sp - tmax_sp_30,
           tmax_6_a = tmax_6 - tmax_6_30,
           tmax_12_a = tmax_12 - tmax_12_30)
  
  # Where are the monitored velvet mesquites?
  mesq_locs <- flower_mesq %>%
    left_join(distinct(weathervars, site_id, elev), by = "site_id") %>%
    group_by(site_id, lat, lon, elev) %>%
    summarize(n_plants = n_distinct(id),
              n_yrs = n_distinct(yr),
              .groups = "keep") %>%
    left_join(normvars, by = "site_id") %>%
    data.frame()
  
  # Geographic location
  # ggplot(mesq_locs, aes(x = lon, y = lat)) +
  #   geom_point(aes(size = n_plants)) +
  #   scale_size_continuous(breaks = c(1, 5, 10), range = c(2, 5))
  leaflet(mesq_locs) %>% addTiles() %>%
    addCircleMarkers(~lon, ~lat, radius = ~n_plants, fillOpacity = 0.6)
  
  # Elevation, Climate norms
  ggplot(mesq_locs, aes(x = elev)) +
    geom_histogram(bins = 30, fill = "gray") +
    labs(x = "Elevation (m)", y = "No. sites") +
    geom_vline(aes(xintercept = mean(elev)), color = "steelblue4",
               linetype = "dashed") +
    annotate(geom = "text", label = "Mean", color = "steelblue4", 
             x = mean(mesq_locs$elev) - 5, y = 9.5, hjust = 1) +
    theme_bw()
  ggplot(mesq_locs) +
    geom_point(aes(x = tmean_30, y = ppt_30, color = elev),
               size = 2) +
    scale_color_viridis_c() +
    labs(x = "Mean annual temperature (degC)", 
         y = "Mean annual precipitation (mm)",
         color = "Elevation (m)") +
    theme_bw()
  
  # Add weather data to flowering dataset (and standardize)
  flower_mesq <- flower_mesq %>%
    left_join(weathervars, by = c("site_id" = "site_id", "yr" = "seasonyr")) %>%
    mutate(tmin_wi_z = (tmin_wi_a - mean(tmin_wi_a)) / sd(tmin_wi_a),
           freeze_z = (freeze - mean(freeze)) / sd(freeze),
           ppt_fa_z = (ppt_fa_a - mean(ppt_fa_a)) / sd(ppt_fa_a),
           ppt_6_z = (ppt_6_a - mean(ppt_6_a)) / sd(ppt_6_a),
           ppt_9_z = (ppt_9_a - mean(ppt_9_a)) / sd(ppt_9_a),
           tmax_su_z = (tmax_su_a - mean(tmax_su_a)) / sd(tmax_su_a),
           tmax_fa_z = (tmax_fa_a - mean(tmax_fa_a)) / sd(tmax_fa_a),
           tmax_wi_z = (tmax_wi_a - mean(tmax_wi_a)) / sd(tmax_wi_a),
           tmax_sp_z = (tmax_sp_a - mean(tmax_sp_a)) / sd(tmax_sp_a),
           tmax_6_z = (tmax_6_a - mean(tmax_6_a)) / sd(tmax_6_a),
           tmax_12_z = (tmax_12_a - mean(tmax_12_a)) / sd(tmax_12_a),
           ppt_norm_z = (ppt_30 - mean(ppt_30)) / sd(ppt_30),
           tmean_norm_z = (tmean_30 - mean(tmean_30)) / sd(tmean_30),
           elev_z = (elev - mean(elev)) / sd(elev))
  
  # Run ML ordinal models with plant random effects
  # Year with plant RE
  m_yr1 <- clmm(abund3 ~ fyr + (1|id), Hess = TRUE,
                data = flower_mesq)
  summary(m_yr1)
  # Year with plant and site REs
  m_yr2 <- clmm(abund3 ~ fyr + (1|id) + (1|site), Hess = TRUE,
                data = flower_mesq)
  summary(m_yr2) # Site effect is 0, so only plant RE moving forward
  # Re-run 1 RE model with quadrature methods
  m_yr <- clmm(abund3 ~ fyr + (1|id), Hess = TRUE, nAGQ = 10,
               data = flower_mesq)
  summary(m_yr)
  # Just RE, no year
  m_1 <- clmm(abund3 ~ 1 + (1|id), Hess = TRUE, nAGQ = 10,
               data = flower_mesq)
  anova(m_1, m_yr)
  # No evidence of annual variation
  
  # Elevation
  m_elev <- clmm(abund3 ~ elev_z + (1|id), Hess = TRUE, nAGQ = 10,
                 data = flower_mesq)
  summary(m_elev)
  anova(m_elev, m_1)
  # Negative elevation effect (P = 0.10)

  # Winter minimum temperatures
  m_tmin <- clmm(abund3 ~ tmin_wi_z + elev_z + (1|id),  
                 Hess = TRUE, nAGQ = 10, data = flower_mesq)
  summary(m_tmin)
  m_tmini <- clmm(abund3 ~ tmin_wi_z * elev_z + (1|id),
                  Hess = TRUE, nAGQ = 10, data = flower_mesq)
  summary(m_tmini)
  anova(m_elev, m_tmin, m_tmini) # P = 0.05, 0.11
  AIC(m_elev, m_tmin, m_tmini) # Interaction has lowest AIC
  # Negative winter min temperature effect, some evidence that it varies with loc
  
  # Precipitation
  m_ppt_fa <- clmm(abund3 ~ ppt_fa_z + elev_z + (1|id),
                   Hess = TRUE, nAGQ = 10, data = flower_mesq)
  m_ppt_6 <- clmm(abund3 ~ ppt_6_z + elev_z + (1|id),
                  Hess = TRUE, nAGQ = 10, data = flower_mesq)
  m_ppt_9 <- clmm(abund3 ~ ppt_9_z + elev_z + (1|id),
                  Hess = TRUE, nAGQ = 10, data = flower_mesq)
  AIC(m_elev, m_ppt_fa, m_ppt_6, m_ppt_9)
  summary(m_ppt_9)
  anova(m_elev, m_ppt_9) # Slight positive effect of precip
  
  # Does precip effect vary with location?
  m_ppt_9i <- clmm(abund3 ~ ppt_9_z * elev_z + (1|id),
                   Hess = TRUE, nAGQ = 10, data = flower_mesq)
  summary(m_ppt_9i)
  anova(m_ppt_9, m_ppt_9i)
  # No strong evidence of spatial variation
  
  # Maximum temperatures
  m_tmax_su <- clmm(abund3 ~ tmax_su_z + elev_z + (1|id), 
                    Hess = TRUE, nAGQ = 10, data = flower_mesq)
  m_tmax_fa <- clmm(abund3 ~ tmax_fa_z + elev_z + (1|id), 
                    Hess = TRUE, nAGQ = 10, data = flower_mesq)
  m_tmax_wi <- clmm(abund3 ~ tmax_wi_z + elev_z + (1|id), 
                    Hess = TRUE, nAGQ = 10, data = flower_mesq)
  m_tmax_6 <- clmm(abund3 ~ tmax_6_z + elev_z + (1|id), 
                   Hess = TRUE, nAGQ = 10, data = flower_mesq)
  AIC(m_elev, m_tmax_su, m_tmax_fa, m_tmax_wi, m_tmax_6)
  summary(m_tmax_fa)
  # Fall slightly better model, with temperatures having a negative effect
  
  # Do temperature effects vary with location?
  m_tmax_fai <- clmm(abund3 ~ tmax_fa_z * elev_z + (1|id), 
                     Hess = TRUE, nAGQ = 10, data = flower_mesq)
  summary(m_tmax_fai)
  anova(m_tmax_fa, m_tmax_fai)
  # No strong evidence of spatial variation
  
  # Full model
  m_full <- clmm(abund3 ~ elev_z + tmin_wi_z + elev_z:tmin_wi_z +
                   ppt_9_z + tmax_fa_z + (1|id),  
                 Hess = TRUE, nAGQ = 10, data = flower_mesq)
  summary(m_full)
  anova(m_full, m_tmini)
  AIC(m_elev, m_tmini, m_full)
  # Minimum temperatures model better than full
  
  # Re-run model with clmm2 (for predict)
  m_tmini2 <- clmm2(abund3 ~ tmin_wi_z * elev_z, random = id, nAGQ = 10, 
                    Hess = TRUE, data = flower_mesq)
  # Create dataframe for prediction
  newdat <- expand.grid(
    abund3 = unique(flower_mesq$abund3),
    elev_z = c(min(flower_mesq$elev_z), 0, 
               max(flower_mesq$elev_z)),
    tmin_wi_z = seq(min(flower_mesq$tmin_wi_z), max(flower_mesq$tmin_wi_z), 
                    length = 100),
    KEEP.OUT.ATTRS = FALSE
  )
  # Plot predictions (= probability that saguaro will have X flowers
  # given winter min temperatures for a cold/average/warm site)
  preds <- cbind(newdat, est = predict(m_tmini2, newdata = newdat)) %>%
    arrange(abund3, elev_z, tmin_wi_z) %>%
    mutate(loc = case_when(
      elev_z == min(flower_mesq$elev_z) ~ "low elev",
      elev_z == 0 ~ "average elev",
      elev_z == max(flower_mesq$elev_z) ~ "high elev",
    )) %>%
    mutate(loc = factor(loc, levels = c("high elev", "average elev", "low elev")))
  ppm <- ggplot(preds, aes(x = tmin_wi_z, y = est)) +
    geom_line(aes(color = abund3), linewidth = 1.3) +
    scale_color_manual(values = c("#a6611a", "#80cdc1", "#018571"), 
                       labels = c("0-100", "101-1000", ">1000")) +
    facet_grid(loc ~ .) +
    labs(x = "Winter minimum temperature anomaly (standardized)", 
         y = "Probability", 
         color = "Flowers",
         title = "Velvet mesquite flowers") +
    theme_bw() +
    theme(legend.position = "bottom")
  # ggsave("output/mesquite-flower-preds-elev.png",
  #        ppm, height = 6.5, width = 6.5, units = "in")
  
  # jojoba --------------------------------------------------------------------#
  # Tried some models for jojoba, but things were coming out odd.
  # Strong, negative effect of elevation, but that could be driven by one very
  # high-elevation site that monitored 3 plants in one year
  # Strong effect of spring temperatures, but jojoba can flower almost anytime, 
  # so in some plants those temperatures are occurring after flowering began.
  # Leaving jojoba for now.
  
# -----------------------------------------------------------------------------#
# Exploring correlations between flower and fruit counts ----------------------#

# Going back to si dataframe (before we did some filtering)
# Summarize amount and quality of information for each plant, phenophase, year
ppy <- si %>%
  left_join(select(timings, common_name, phenophase_description, 
                   q0.10, q0.25, q0.75, q0.90), 
            by = c("common_name", "phenophase_description")) %>%
  mutate(in50 = ifelse(day_of_year >= q0.25 & day_of_year <= q0.75, 1, 0)) %>%
  mutate(in80 = ifelse(day_of_year >= q0.10 & day_of_year <= q0.90, 1, 0)) %>%
  group_by(common_name, individual_id, site_id, phenophase_description, yr) %>%
  summarize(obs = n(),
            first_obs = min(day_of_year),
            last_obs = max(day_of_year),
            obs50 = sum(in50),
            obs80 = sum(in80),
            mean_int = round(mean(interval, na.rm = TRUE), 2),
            max_int = ifelse(obs == 1, NA, max(interval, na.rm = TRUE)),
            inphase = sum(phenophase_status),
            intvalue = sum(!is.na(intensity_value)),
            prop_inphase = round(inphase / obs, 2),
            prop_intvalue = round(intvalue / inphase, 2),
            max_intensity = ifelse(sum(is.na(intensity_midpoint)) == n(),
                                   NA, max(intensity_midpoint, na.rm = TRUE)),
            .groups = "keep") %>%
  data.frame()

# Summarize data on flowers/fruits for each plant and year
ff_ppy <- ppy %>%
  filter(phenophase_description %in% c("Flowers or flower buds", "Fruits")) %>%
  mutate(php = ifelse(phenophase_description == "Fruits", "fr", "fl")) %>%
  pivot_wider(
    id_cols = c(common_name, individual_id, yr),
    names_from = php,
    names_glue = "{php}_{.value}",
    values_from = c(obs, inphase, intvalue, obs50, obs80, max_intensity)
  ) %>%
  data.frame() %>%
  # Convert any NAs to 0s for number of observations (happened when there were
  # observations for one phenophase but not the other in a given year)
  mutate(across(fl_obs:fr_obs80, \(x) replace_na(x, 0)))

# Identify plant-years that have a sufficient number of observations of both 
# during appropriate periods that would make it worth evaluating correlations

# Try same filters that we used for ordinal regression models? (remove 
# plant-php-year combos when fewer than 5 observations were made between dates
# associated with the 10th and 90th percentile AND fewer than 2 observations
# were made between dates associated with the 25th and 75th percentiles)
ff_ppy <- ff_ppy %>%
  mutate(remove = ifelse(
    fl_obs50 < 2 | fr_obs50 < 2 | fl_obs50 < 5 | fr_obs50 < 5, 1, 0
  ))
count(ff_ppy, remove)
# This filter would remove 69% of plant-years...

# Filter out plant-years with too few observations, AND filter out plant-years
# where the max number of flowers is 0 but the max number of fruit is > 0 (which
# doesn't occur that often). Note: by excluding any plant-years with 0 flowers,
# we can get rid of those instances and also when flower and fruit = 0, which
# we don't care about
fff_ppy <- ff_ppy %>%
  filter(remove == 0) %>%
  filter(fl_max_intensity > 0) %>%
  rename(fl_max = fl_max_intensity,
         fr_max = fr_max_intensity)

# Plot logged max intensity values (need to add a small value to fruits, that
# can have max = 0)
corplot <- fff_ppy %>%
  filter(common_name != "California barrel cactus")  %>%
  mutate(fl_max_l = log(fl_max),
         fr_max_l = ifelse(fr_max == 0, log(fr_max + 0.1), log(fr_max))) %>%
  ggplot(aes(x = fl_max_l, y = fr_max_l)) +
  geom_jitter(width = 0.1, height = 0.1) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  facet_wrap(~ common_name, scales = "free_x") +
  labs(x = "log(max flowers)", y = "log(max fruit)") + 
  ggpubr::stat_cor(method = "pearson", cor.coef.name = "r")
# ggsave("output/flower-fruit-cor.png", corplot,
#        height = 5, width = 6.5, units = "in")

# High correlations between max number of flowers and fruit in cholla (0.71) 
# and saguaro (0.67), less so with mesquite (0.50) and almost none with 
# jojoba (0.19). 

# -----------------------------------------------------------------------------#
# How does timing of max #flowers relate to flower/open flower phenophases? ---#
# (need to have run though at least line 615)

# Might need to look at data near Phoenix and Tucson separately
tucson <- data.frame(lon = -110.9742, lat = 32.2540)
phoenix <- data.frame(lon = -112.0777, lat = 33.4482)

sag_locs <- sag_locs %>% 
  mutate(dist_tuc = geosphere::distHaversine(cbind(lon, lat),
                                             cbind(tucson$lon, tucson$lat))) %>%
  mutate(dist_tuc = dist_tuc / 1000) %>%
  mutate(tucson = ifelse(dist_tuc < 35, 1, 0)) %>%
  mutate(dist_phx = geosphere::distHaversine(cbind(lon, lat),
                                             cbind(phoenix$lon, phoenix$lat))) %>%
  mutate(dist_phx = dist_phx / 1000) %>%
  mutate(phoenix = ifelse(dist_phx < 50, 1, 0))

# Check that geographic filters look ok
leaflet(sag_locs) %>% addTiles() %>%
  addCircleMarkers(~lon, ~lat, radius = ~n_plants, fillOpacity = 0.6) %>%
  addCircleMarkers(~lon, ~lat, radius = ~n_plants, fillOpacity = 0.6,
                   data = filter(sag_locs, tucson == 1), color = "red") %>%
  addCircleMarkers(~lon, ~lat, radius = ~n_plants, fillOpacity = 0.6,
                   data = filter(sag_locs, phoenix == 1), color = "red")

# Will use flowers dataframe for max counts and will use si dataframe for
# phenophase status because we don't want as many strict filters for that 
# stuff.... Use data from 2012 on

# Get timing of max flower counts (in each region)
maxf_sag <- flowers %>%
  filter(common_name == "saguaro") %>%
  filter(yr >= 2012) %>%
  # Need to recalculate Tucson/Phoenix IDs since there are a few sites that
  # weren't included before (because there were < 10 plants monitored that year)
  # but there's no reason to exclude them for this.
  mutate(dist_tuc = geosphere::distHaversine(cbind(lon, lat),
                                             cbind(tucson$lon, tucson$lat))) %>%
  mutate(dist_tuc = dist_tuc / 1000) %>%
  mutate(tucson = ifelse(dist_tuc < 35, 1, 0)) %>%
  mutate(dist_phx = geosphere::distHaversine(cbind(lon, lat),
                                             cbind(phoenix$lon, phoenix$lat))) %>%
  mutate(dist_phx = dist_phx / 1000) %>%
  mutate(phoenix = ifelse(dist_phx < 50, 1, 0)) %>%
  rename(intensity = intensity_midpoint)

# Now, for each plant and year, identify the max count and then for all non-0s,
# the date(s) that the max count was observed
maxdates <- maxf_sag %>%
  group_by(id, site_id, lat, lon, tucson, phoenix, yr) %>%
  mutate(nobs = n(),
         max = ifelse(nobs == 1, NA, max(intensity, na.rm = TRUE))) %>%
  filter(!is.na(max) & max > 0) %>%
  select(-nobs) %>%
  group_by(id, site_id, lat, lon, tucson, phoenix, yr, max) %>%
  summarize(max_first = min(doy[intensity == max]),
            max_last = max(doy[intensity == max]),
            max_mn = mean(doy[intensity == max]),
            .groups = "keep") %>%
  data.frame()

# Remove observations of saguaros outside Tucson/Phoenix since they're often at
# different elevations or in different climates
maxdates <- maxdates %>%
  filter(tucson == 1 | phoenix == 1)

# Does it look like there's a difference between Tucson and Phoenix?
maxdates <- maxdates %>%
  mutate(fyr = as.factor(yr - min(yr)))
m_first <- lm(max_first ~ tucson + fyr, data = maxdates)
m_last <- lm(max_mn ~ tucson + fyr, data = maxdates)
m_mn <- lm(max_last ~ tucson + fyr, data = maxdates)
summary(m_first)$coefficients["tucson",]
summary(m_mn)$coefficients["tucson",]
summary(m_last)$coefficients["tucson",]
# Tucson does look to be 4-5 days earlier

# Look at distribution of dates across plants, years, location
maxdates %>%
  pivot_longer(cols = max_first:max_mn,
               names_to = "metric",
               values_to = "doy") %>%
  ggplot() +
    geom_histogram(aes(x = doy), bins = 50) +
    facet_grid(metric ~.)
# Just Tucson
maxdates %>%
  filter(tucson == 1) %>%
  pivot_longer(cols = max_first:max_mn,
               names_to = "metric",
               values_to = "doy") %>%
  ggplot() +
  geom_histogram(aes(x = doy), bins = 50) +
  facet_grid(metric ~.)
# As expected, but will probably want to exclude some outliers from analysis

# Look at distribution of mean dates for different years in Tucson
maxdates %>%
  filter(tucson == 1) %>%
  filter(yr %in% 2018:2024) %>%
  ggplot() +
    geom_histogram(aes(x = max_mn), bins = 50) +
    facet_grid(yr ~ .)
# 2024 looks like it was later, otherwise not a whole lot of difference

# Get phenophase status data (from Tucson area from 2012-2024)
sag_tuc <- si %>%
  filter(common_name == "saguaro",
         class_id %in% flower_classes, 
         yr %in% 2012:2024) %>%
  rename(lat = latitude, 
         lon = longitude,
         elev = elevation_in_meters,
         id = individual_id,
         doy = day_of_year,
         status = phenophase_status) %>%
  mutate(dist_tuc = geosphere::distHaversine(cbind(lon, lat),
                                             cbind(tucson$lon, tucson$lat))) %>%
  mutate(dist_tuc = dist_tuc / 1000) %>%
  mutate(tucson = ifelse(dist_tuc < 35, 1, 0)) %>%
  filter(tucson == 1) %>%
  select(site_id, id, lat, lon, phenophase_description, observation_date, yr, 
         doy, status, intensity_value, intensity_midpoint, interval)

# Add week of the year
sag_tuc <- sag_tuc %>%
  mutate(wk = week(observation_date)) %>%
  filter(wk < 53)

# Calculate the weekly proportion of individuals in each phenophase after
# removing all but one observation of a plant each week
wkprop <- sag_tuc %>%
  mutate(php = ifelse(phenophase_description == "Open flowers", 
                      "open", "flower")) %>%
  arrange(id, php, yr, wk, desc(status)) %>%
  distinct(id, php, yr, wk, .keep_all = TRUE) %>%
  group_by(php, yr, wk) %>%
  summarize(nobs = n(),
            inphase = sum(status),
            .groups = "keep") %>%
  data.frame() %>%
  mutate(prop = inphase / nobs) 

wkprop %>%
  group_by(yr) %>%
  summarize(min_nobs = min(nobs),
            mn_nobs = round(mean(nobs), 1))
# Mean number of observations per week much greater from 2017-2024

wkprop %>%
  group_by(wk) %>%
  summarize(min_nobs = min(nobs),
            mn_nobs = round(mean(nobs), 1)) %>%
  data.frame()
# More observations per week in spring, but not a huge dropoff in winter

# Exclude data before 2017, and then remove any weekly proportions that are 
# based on < 10 individuals per week
wkprop <- wkprop %>%
  filter(yr >= 2017) %>%
  filter(nobs >= 10) %>%
  mutate(fyr = factor(yr))

# Plot annual values
ggplot(wkprop, aes(x = wk, y = prop)) +
  geom_line(aes(color = fyr), alpha = 0.5) +
  geom_point(aes(color = fyr, size = nobs), alpha = 0.5) +
  facet_grid(php ~ .)

# Fit GAMs for phenophase flower status each year:
floweryr <- gam(prop ~ s(wk, by = fyr, bs = "cc", k = 20), 
                family = "binomial",
                weights = nobs,
                method = "REML",
                data = filter(wkprop, php == "flower"))
summary(floweryr)
gam.check(floweryr) 
# p-values are low, but do not improve (and predictions don't meaningfully 
# change) when k is higher
marginaleffects::plot_predictions(floweryr, by = c("wk", "fyr", "fyr"))

flowerpreds <- expand.grid(
  wk = 1:52,
  fyr = 2017:2024,
  KEEP.OUT.ATTRS = FALSE
  ) %>%
  mutate(fyr = factor(fyr))
predsfl <- predict(floweryr, newdata = flowerpreds, 
                   type = "link", se.fit = TRUE)

flowerpreds$est_l <- predsfl$fit
flowerpreds$se <- predsfl$se.fit
flowerpreds <- flowerpreds %>%
  mutate(lwr_l = est_l - 2 * se,
         upr_l = est_l + 2 * se) %>%
  mutate(estimate = exp(est_l) / (1 + exp(est_l)),
         lwr = exp(lwr_l) / (1 + exp(lwr_l)),
         upr = exp(upr_l) / (1 + exp(upr_l)))

ggplot(data = flowerpreds, aes(x = wk)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.5) +
  geom_line(aes(y = estimate)) +
  facet_wrap(~fyr, scales = "fixed")

# Fit GAMs for phenophase open open status each year:
openyr <- gam(prop ~ s(wk, by = fyr, bs = "cc", k = 20), 
                family = "binomial",
                weights = nobs,
                method = "REML",
                data = filter(wkprop, php == "open"))
summary(openyr)
gam.check(openyr) 
# p-values are low, but do not improve (and predictions don't meaningfully 
# change) when k is higher
marginaleffects::plot_predictions(openyr, by = c("wk", "fyr", "fyr"))

openpreds <- expand.grid(
  wk = 1:52,
  fyr = 2017:2024,
  KEEP.OUT.ATTRS = FALSE
) %>%
  mutate(fyr = factor(fyr))
predsop <- predict(openyr, newdata = openpreds, 
                   type = "link", se.fit = TRUE)

openpreds$est_l <- predsop$fit
openpreds$se <- predsop$se.fit
openpreds <- openpreds %>%
  mutate(lwr_l = est_l - 2 * se,
         upr_l = est_l + 2 * se) %>%
  mutate(estimate = exp(est_l) / (1 + exp(est_l)),
         lwr = exp(lwr_l) / (1 + exp(lwr_l)),
         upr = exp(upr_l) / (1 + exp(upr_l)))

ggplot(data = openpreds, aes(x = wk)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.5) +
  geom_line(aes(y = estimate)) +
  facet_wrap(~fyr, scales = "fixed")

# Next steps:
# Plot flower, open flower, and max dates together (for each year)