# Exploring intensity data from McDowell
# ER Zylstra

library(dplyr)
library(stringr)
library(lubridate)
library(ggplot2)
library(MASS)

# Plant phenophase classes
leaf_classes <- 1:5
flower_classes <- 6:9
fruit_classes <- 10:13

# List files with formatted intensity data from McDowell
intensity_files <- list.files("npn-data",
                              pattern = "intensity-siteMCDO",
                              full.names = TRUE)

# Load and format status-intensity data ---------------------------------------#

sites <- str_sub(intensity_files, 24, 28)

for (site in sites) {
  filename <- intensity_files[grepl(site, intensity_files)]
  si1 <- read.csv(filename)
  
  # Remove any true duplicates
  si1 <- si1[!duplicated(si1),]
  
  if (site == sites[1]) {
    si <- si1
  } else {
    si <- rbind(si, si1)
  }
  rm(si1)
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

# si %>% count(interval) %>% mutate(prop = n / 270135) %>% round(3)
# Vast majority (97%) of intervals are <= 7 days
# Most intervals (77%) <= 3 days

# Check that there's only one intensity category for each species-phenophase?
spil <- si %>%
  filter(!is.na(intensity_category_id)) %>%
  distinct(common_name, phenophase_description, intensity_name,
           intensity_type, intensity_label)
count(spil, common_name, phenophase_description) %>%
  pull(n) %>% 
  max()

# Make intensity midpoints = 0 if status = 0 and add intensity category labels
si <- si %>%
  mutate(intensity_midpoint = ifelse(phenophase_status == 0, 
                                     0, intensity_midpoint)) %>%
  dplyr::select(-c(intensity_name, intensity_type, intensity_label)) %>%
  left_join(spil, by = c("common_name", "phenophase_description"))
# Now the only NAs left in the intensity_midpoint column occur when the status 
# is yes but no intensity value was provided or for the Pollen release 
# phenophase

# Load plant details ----------------------------------------------------------#

# # Base URL for getting plant details through NPN API
# url_base <- "https://services.usanpn.org/npn_portal/individuals/getPlantDetails.xml?individual_id[0]=250&individual_id[1]=360"
# 
# # get vector of plant IDs
# plantids <- unique(si$individual_id)

# Was going to download plant details, but I had this file already from POP. 
# There is no information about size, age, or number of arms (though if we're
# going to use this as a case study in a paper, I'd definitely want to get this
# information through Mary F.). 

plantdetails <- read.csv("npn-data/mcdo-ancillary-individual-plant-data.csv") %>%
  dplyr::select(-c(Plant_Image_URL, Plant_Image_Upload_Date))

# One thing worth noting is when plants died. Will have entries in at least one of:
# Death_Reason, Death_Date_Observed, or Last_Date_Observed_Alive

# Two saguaros died:
# #117055: died from drought on 3/14/2024
# #117062: died of unknown cases, last seem alive on 11/16/2021

# Load weather data -----------------------------------------------------------#

# See mcdowell-weather-data.R for details about PRISM data download

weather <- read.csv("weather-data/mcdo-20162024.csv") %>%
  dplyr::select(-c(lat, lon))

# Filter data: plant-phenophase-year combinations -----------------------------#

# Want to:
  # remove any pollen release data (qualitative data)
  # remove plant-php-year combos with < 10 observations
  # remove plant-php-year combos with maximum interval >21 days 

  # For analyses of annual max counts, we want to include max annual counts of 
  # zero AS LONG AS observations were made within the period that the phenophase 
  # is typically observed. So, for each species and phenophase, we aggregated
  # all dates when plants at McDowell were considered in phase, and then 
  # elected to:
  # remove any plant-php-year combos where fewer than 5 observations were made 
  # between dates associated with the 10th and 90th percentile AND fewer than 2 
  # observations were made between dates associated with the 25th and 75th 
  # percentiles.

# Identify when species at McDowell are typically in various phenophases
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
# timings

# Summarize amount and quality of information for each plant, phenophase, year
pl_ph_yr <- si %>%
  left_join(dplyr::select(timings, common_name, phenophase_description, 
                          q0.10, q0.25, q0.75, q0.90), 
            by = c("common_name", "phenophase_description")) %>%
  mutate(in50 = ifelse(day_of_year >= q0.25 & day_of_year <= q0.75, 1, 0)) %>%
  mutate(in80 = ifelse(day_of_year >= q0.10 & day_of_year <= q0.90, 1, 0)) %>%
  group_by(common_name, individual_id, phenophase_description, yr) %>%
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

# Identify which plant-php-yr combinations could be removed 
pl_ph_yr <- pl_ph_yr %>%
  mutate(remove = case_when(
    phenophase_description == "Pollen release (flowers)" ~ 1,
    nobs < 10 ~ 1,
    max_int > 21 ~ 1,
    is.na(max_intensity) ~ 1,
    max_intensity == 0 & nobs80 < 5 ~ 1,
    max_intensity == 0 & nobs50 < 2 ~ 1,
    .default = 0
  ))
# check:
# count(filter(pl_ph_yr, remove == 1 & 
#                phenophase_description != "Pollen release (flowers)"), 
#       nobs < 10, max_int > 21, is.na(max_intensity),
#       (max_intensity == 0 & !is.na(max_intensity)), nobs80 < 5, nobs50 < 2)

si <- si %>%
  left_join(dplyr::select(pl_ph_yr, individual_id, phenophase_description, yr, 
                          remove),
            by = c("individual_id", "phenophase_description", "yr")) %>%
  filter(remove == 0) %>%
  dplyr::select(-remove)

# Summarize data by species ---------------------------------------------------#

# Species summaries
spp_summary <- si %>%
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
# write.table(spp_summary, "clipboard", sep = "\t", row.names = FALSE)

# Summarize data by intensity category ----------------------------------------#

# Create table summarizing amount of information per intensity category
intensity_cats <- si %>%
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
  dplyr::select(-c(values, valuesq)) %>%
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

# Summarize filtered data for each plant, php, and year -----------------------#

# Recreating pl_ph_yr using rle(phenophase_status) to understand patterns:
  # length(rle$values) = number of 0/1 status sequences
  # first(rle$values) = state at first observation
  # last(rle$values) = state at last observation
pl_ph_yr <- si %>%
  arrange(site_name, individual_id, phenophase_id, observation_date) %>%
  group_by(common_name, site_name, individual_id, class_id, 
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
# write.table(dplyr::select(php_summary, -class_id),
#             "clipboard", sep = "\t", row.names = FALSE)

# Quite a bit more fluctuation in phenophase status than there was in other 
# datasets (Kodiak, NEON), which likely reflects environment and how desert
# plants respond to weather changes. It's unclear the extent to which this 
# variation also reflects observer quality.

# Plot raw intensity data -----------------------------------------------------#

# Create column with logged intensity values
si <- si %>%
  mutate(intensity_midpoint_log = ifelse(intensity_midpoint == 0,
                                         log(intensity_midpoint + 0.01), 
                                         log(intensity_midpoint)))

# Loop through intensity categories and then species
for (i in 1:nrow(intensity_cats)) {
  
  si_int <- filter(si, intensity_label == intensity_cats$intensity_label[i])
  spps <- unique(si_int$common_name)
  
  # Remove observations where phenophase status is 1, but intensity value wasn't
  # provided (this creates breaks in plotted lines)
  si_int <- filter(si_int, !is.na(intensity_midpoint))
  
  # Plot logged values for numeric intensity categories
  if (si_int$intensity_type[1] == "number") {
    si_int <- si_int %>%
      mutate(yaxis = intensity_midpoint_log)
  } else {
    si_int <- si_int %>%
      mutate(yaxis = intensity_midpoint)
  }
  
  for (spp in spps) {
    
    # Create filename for png
    spp_nospace <- str_replace_all(spp, " ", "_")
    png_name <- paste0("output/mcdo-intensities-by-spp-php-yr/IntensityData-",
                       intensity_cats$intensity_short[i], "-", spp_nospace, 
                       ".png")
    
    # Create ggplot object and save to file (if file doesn't already exist)
    if (!file.exists(png_name)) {
      
      si_int_spp <- filter(si_int, common_name == spp)
      
      ylab <- ifelse(si_int$intensity_type[1] == "number",
                     paste0("log(", str_to_lower(si_int$intensity_label[1]), ")"),
                     si_int$intensity_label[1])
      
      # Figures with lines and points, with size proportional to no. of plants 
      agg <- si_int_spp %>%
        group_by(yr, day_of_year, yaxis) %>%
        summarize(n_indiv = n(), .groups = "keep") %>%
        data.frame()
      
      iplot <- ggplot(agg, aes(x = day_of_year, y = yaxis)) +
        facet_grid(yr ~ .) +
        geom_line(data = si_int_spp, show.legend = FALSE, alpha = 0.5,
                  aes(x = day_of_year, y = yaxis, color = factor(individual_id))) +
        geom_point(data = agg,
                   aes(x = day_of_year, y = yaxis, size = n_indiv),
                   shape = 16, alpha = 0.4) +
        scale_size_continuous(breaks = 1:4) +
        labs(title = paste0(str_to_sentence(spp), ",  ", intensity_cats$intensity_label[i]),
             y = ylab, x = "Day of year") +
        theme_bw()
      iplot
      
      ggsave(png_name,
             plot = iplot,
             width = 6.5,
             height = 9,
             units = "in")
    }
  }
}

# Additional data filtering, formatting for GAM models ------------------------#

# Removing data from 2016, since monitoring began late in the year
si <- filter(si, yr > 2016)

# Use a rule to exclude any species-phenophase combination where there's only 
# one year of data. 
si <- si %>%
  group_by(common_name, phenophase_description) %>%
  mutate(n_yrs = n_distinct(yr)) %>%
  ungroup() %>%
  filter(n_yrs > 1) %>%
  dplyr::select(-n_yrs) %>%
  data.frame()

# Prep data for regression models ---------------------------------------------#

# Sort data and remove any observations with missing intensity values
gamdf <- si %>%
  arrange(class_id, phenophase_description, common_name, site_name, 
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
combos <- gamdf %>%
  mutate(plantyr = paste0(individual_id, "_", yr)) %>%
  group_by(class_id, intensity_label, intensity_type, common_name) %>%
  summarize(n_plantyrs = n_distinct(plantyr),
            n_plants = n_distinct(individual_id),
            n_yrs = n_distinct(yr),
            n_sites = n_distinct(site_name),
            .groups = "keep") %>%
  data.frame()
# write.table(dplyr::select(combos, -class_id), "clipboard",
#             sep = "\t", row.names = FALSE)

# -----------------------------------------------------------------------------#
# Exploring variation in max counts for flowers phenophase --------------------#

flowers <- gamdf %>%
  filter(intensity_short == "Flowers") %>%
  dplyr::select(common_name, individual_id, site_name, latitude, longitude, 
                phenophase_description, observation_date, yr, day_of_year,
                phenophase_status, intensity_label, 
                intensity_midpoint) %>%
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

# Would like to model data for multiple species simultaneously, but they have 
# different intensity categories, so I'd have to create bins that work for all. 
# Max category == "More than 1,000": Barrel cacus, cassia, cholla, saguaro
# Max category == "More than 10,000": jojoba, soaptree yucca, velvet mesquite
# Might be able to combine anything over 1000 since there are relatively few
# max counts in that range for all species except jojoba, but will revisit this
# later

 # buck-horn cholla -----------------------------------------------------------#
  # Aggregate data for each plant-year: buck-horn cholla
  flowers_py_cholla <- flowers %>%
    filter(common_name == "buck-horn cholla") %>%
    group_by(common_name, id, site_name, lat, lon, yr) %>%
    summarize(n_obs = n(),
              n_inphase = sum(status),
              max_count = max(intensity_midpoint),
              .groups = "keep") %>%
    data.frame()
  flowers_py_cholla
  count(filter(gamdf, common_name == "buck-horn cholla" & 
                 intensity_short == "Flowers"),
        intensity_midpoint, intensity_value)
  count(flowers_py_cholla, max_count)
  # 0:                        1/70
  # Less than 3 (1):          5/70
  # 3 to 10 (5):              2/70
  # 11 to 100 (50):          27/70
  # 101 to 1000 (500):       31/70
  # More than 1,000 (1001):   4/70
  
  # Since there's only one zero, will have to group with few
  # 0-5 (few), 50 (some), 500-1001 (many)
  flowers_py_cholla <- flowers_py_cholla %>%
    mutate(abund = case_when(
      max_count == 0 ~ "few",
      max_count == 1 ~ "few",
      max_count == 5 ~ "few",
      max_count == 50 ~ "some",
      max_count >= 500 ~ "many",
    )) %>%
    mutate(abund = factor(abund, levels = c("few", "some", "many")))
  
  # Visualize for each year
  flowers_pya_cholla <- flowers_py_cholla %>%
    mutate(fyr = factor(yr)) %>%
    mutate(site = factor(site_name)) %>%
    group_by(common_name, fyr, site, abund) %>%
    summarize(n = n(), .groups = "keep") %>%
    data.frame()
  ggplot(flowers_pya_cholla, aes(y = abund, x = fyr)) +
    geom_point(aes(size = n)) +
    facet_grid(site ~ .) +
    scale_y_discrete(labels = c("0-5", "50", "500+")) +
    scale_size_continuous(breaks = 1:3, range = c(2, 4)) +
    labs(x = "", y = "Abundance category", size = "No. plants", 
         title = "buck-horn cholla, flowers")
  
  # Test run with MASS::polr
  flowers_py_cholla$site <- factor(flowers_py_cholla$site_name)
  flowers_py_cholla$fyr <- factor(flowers_py_cholla$yr)
  
  # By default, logistic regression (could change to probit)
  m1 <- polr(abund ~ site, Hess = TRUE, data = flowers_py_cholla)
    summary(m1)
    # p-values (really only valid for larger sample sizes)
    ctable <- coef(summary(m1))
    p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
    (ctable <- cbind(ctable, "p value" = round(p, 3)))
    # confidence intervals (profiling likelihood function)
    ci <- confint(m1)
    # confidence intervals (assuming normality)
    ci <- confint.default(m1)
    # Odds ratios
    exp(cbind(OR = coef(m1), ci))
  
  # Predicted probabilities (could include new values for continuous covs in model)
  newdat <- data.frame(
    site = sort(unique(flowers_py_cholla$site))
  )
  newdat <- cbind(newdat, round(predict(m1, newdat, type = "probs"), 3))
  newdat  
  
  # Mean rank/category for each species
  summary(lsmeans::lsmeans(m1, pairwise ~ site, mode = "mean"), type="response")$lsmeans
  
  m2 <- polr(abund ~ fyr, Hess = TRUE, data = flowers_py_cholla)
  summary(m2)
  m3 <- polr(abund ~ site + fyr, Hess = TRUE, data = flowers_py_cholla)
  summary(m3)
  m4 <- polr(abund ~ site * fyr, Hess = TRUE, data = flowers_py_cholla)
  summary(m4)
  # Warning about rank-deficient design
  
  AIC(m1, m2, m3)
  # Additive model slightly better than year alone
  
  # Strong evidence that max flower counts differed among years, mostly because
  # counts in 2021 were much lower at all sites (though 2018 was also low). 
  # Counts at Gateway slightly lower
  
  # These models could be much more useful if we replaced year with some 
  # measure(s) of weather at each site and year. Precipitation seems like an 
  # obvious target, particularly given that the 2020 drought was extreme...

  # saguaro -------------------------------------------------------------------#
  # Aggregate data for each plant-year: saguaro
  flowers_py_sag <- flowers %>%
    filter(common_name == "saguaro") %>%
    group_by(common_name, id, site_name, lat, lon, yr) %>%
    summarize(n_obs = n(),
              n_inphase = sum(status),
              max_count = max(intensity_midpoint),
              .groups = "keep") %>%
    data.frame()
  flowers_py_sag
  count(filter(gamdf, common_name == "saguaro" & intensity_short == "Flowers"),
        intensity_midpoint, intensity_value)
  count(flowers_py_sag, max_count)
  # 0:                        1/63
  # Less than 3 (1):          1/63
  # 3 to 10 (5):              7/63
  # 11 to 100 (50):          20/63
  # 101 to 1000 (500):       27/63
  # More than 1,000 (1001):   7/63
  
  # Will try same grouping as buck-horn cholla:
  # 0-5 (few), 50 (some), 500-1001 (many)
  flowers_py_sag <- flowers_py_sag %>%
    mutate(abund = case_when(
      max_count == 0 ~ "few",
      max_count == 1 ~ "few",
      max_count == 5 ~ "few",
      max_count == 50 ~ "some",
      max_count >= 500 ~ "many",
    )) %>%
    mutate(abund = factor(abund, levels = c("few", "some", "many")))
  
  # Visualize for each year
  flowers_pya_sag <- flowers_py_sag %>%
    mutate(fyr = factor(yr)) %>%
    mutate(site = factor(site_name)) %>%
    group_by(common_name, fyr, site, abund) %>%
    summarize(n = n(), .groups = "keep") %>%
    data.frame()
  ggplot(flowers_pya_sag, aes(y = abund, x = fyr)) +
    geom_point(aes(size = n)) +
    facet_grid(site ~ .) +
    scale_y_discrete(labels = c("0-5", "50", "500+")) +
    scale_size_continuous(breaks = 1:3, range = c(2, 4)) +
    labs(x = "", y = "Abundance category", size = "No. plants", 
         title = "saguaro, flowers")
  
  flowers_py_sag$site <- factor(flowers_py_sag$site_name)
  flowers_py_sag$fyr <- factor(flowers_py_sag$yr)
  
  # Models
  m1 <- polr(abund ~ site, Hess = TRUE, data = flowers_py_sag)
    summary(m1)
    # p-values (really only valid for larger sample sizes)
    ctable <- coef(summary(m1))
    p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
    (ctable <- cbind(ctable, "p value" = round(p, 3)))

  m2 <- polr(abund ~ fyr, Hess = TRUE, data = flowers_py_sag)
  summary(m2)
  m3 <- polr(abund ~ site + fyr, Hess = TRUE, data = flowers_py_sag)
  summary(m3)
  # Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred 
  m4 <- polr(abund ~ site * fyr, Hess = TRUE, data = flowers_py_sag)
  summary(m4)
  # Warning about rank-deficient design (didn't encounter this when 
  # 0 max counts were excluded). Could also group 0 with 1-5...
  
  AIC(m1, m2)
  # Site model much better than year (ignoring other 2 models)
  
  # Strong evidence that max flower counts differed among sites, with most 
  # saguaros at Brown's Ranch always having 500+ flowers (would be good to 
  # know something about size/age). Does look like some annual variation, but
  # maybe not consistent among sites and some estimation difficulties because
  # all plants monitored in 2018 had 500+ max counts.

  # jojoba --------------------------------------------------------------------#
  # Aggregate data for each plant-year: jojoba
  flowers_py_jojo <- flowers %>%
    filter(common_name == "jojoba") %>%
    group_by(common_name, id, site_name, lat, lon, yr) %>%
    summarize(n_obs = n(),
              n_inphase = sum(status),
              max_count = max(intensity_midpoint),
              .groups = "keep") %>%
    data.frame()
  flowers_py_jojo
  count(filter(gamdf, common_name == "jojoba" & intensity_short == "Flowers"),
        intensity_midpoint, intensity_value)
  count(flowers_py_jojo, max_count)
  # 0:                        0/82
  # Less than 3 (1):          0/82
  # 3 to 10 (5):              1/82
  # 11 to 100 (50):           7/82
  # 101 to 1000 (500):       34/82
  # 1001 to 10000 (5000):    35/82
  # More than 10,000 (10001): 5/82
  
  # Trying something new here:
  # 0-50 (few), 500 (some), 5000-10001 (many)
  flowers_py_jojo <- flowers_py_jojo %>%
    mutate(abund = case_when(
      max_count %in% 0:50 ~ "few",
      max_count == 500 ~ "some",
      max_count >= 5000 ~ "many",
    )) %>%
    mutate(abund = factor(abund, levels = c("few", "some", "many")))
  
  # Visualize for each year
  flowers_pya_jojo <- flowers_py_jojo %>%
    mutate(fyr = factor(yr)) %>%
    mutate(site = factor(site_name)) %>%
    group_by(common_name, fyr, site, abund) %>%
    summarize(n = n(), .groups = "keep") %>%
    data.frame()
  ggplot(flowers_pya_jojo, aes(y = abund, x = fyr)) +
    geom_point(aes(size = n)) +
    facet_grid(site ~ .) +
    scale_y_discrete(labels = c("0-50", "500", "5000+")) +
    scale_size_continuous(breaks = 1:3, range = c(2, 4)) +
    labs(x = "", y = "Abundance category", size = "No. plants", 
         title = "jojoba, flowers")
  
  flowers_py_jojo$site <- factor(flowers_py_jojo$site_name)
  flowers_py_jojo$fyr <- factor(flowers_py_jojo$yr)
  
  # Models
  m1 <- polr(abund ~ site, Hess = TRUE, data = flowers_py_jojo)
    summary(m1)
    # p-values (really only valid for larger sample sizes)
    ctable <- coef(summary(m1))
    p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
    (ctable <- cbind(ctable, "p value" = round(p, 3)))

  m2 <- polr(abund ~ fyr, Hess = TRUE, data = flowers_py_jojo)
  summary(m2)
  m3 <- polr(abund ~ site + fyr, Hess = TRUE, data = flowers_py_jojo)
  summary(m3)
  m4 <- polr(abund ~ site * fyr, Hess = TRUE, data = flowers_py_jojo)
  summary(m4)
  # Warning about rank-deficient design
  
  AIC(m1, m2, m3)
  # Additive model is best...
  
  # Strong evidence that max flower counts differed among sites, with lower 
  # number of flowers at Lost Dog. Does look like some annual variation, but
  # there are some estimation difficulties because all plants monitored in 2017
  # had 5000+ max counts.

# -----------------------------------------------------------------------------#
# Exploring variation in max counts for fruits phenophase ---------------------#
  
fruit <- gamdf %>%
  filter(intensity_short == "Fruits") %>%
  dplyr::select(common_name, individual_id, site_name, latitude, longitude, 
                phenophase_description, observation_date, yr, day_of_year,
                phenophase_status, intensity_label, 
                intensity_midpoint) %>%
  rename(status = phenophase_status,
         phenophase = phenophase_description,
         id = individual_id,
         lat = latitude,
         lon = longitude,
         obsdate = observation_date, 
         doy = day_of_year)

# Distribution of intensity values
count(fruit, intensity_midpoint) %>%
  mutate(prop = n / sum(n))

count(filter(si, intensity_label == "No. fruits"), 
      common_name, intensity_name, intensity_value)
# Would like to model data for multiple species simultaneously, but they have 
# different intensity categories, so I'd have to create bins that work for all. 
# Max category == "More than 1,000": Barrel cacus, cassia, cholla, saguaro
# Max category == "More than 10,000": jojoba, soaptree yucca, velvet mesquite
  
  # buck-horn cholla -----------------------------------------------------------#
  # Aggregate data for each plant-year: buck-horn cholla
  fruit_py_cholla <- fruit %>%
    filter(common_name == "buck-horn cholla") %>%
    group_by(common_name, id, site_name, lat, lon, yr) %>%
    summarize(n_obs = n(),
              n_inphase = sum(status),
              max_count = max(intensity_midpoint),
              .groups = "keep") %>%
    data.frame()
  fruit_py_cholla
  count(filter(gamdf, common_name == "buck-horn cholla" & 
                 intensity_short == "Fruits"),
        intensity_midpoint, intensity_value)
  count(fruit_py_cholla, max_count)
  # 0:                        2/73
  # Less than 3 (1):          4/73
  # 3 to 10 (5):              6/73
  # 11 to 100 (50):          27/73
  # 101 to 1000 (500):       33/73
  # More than 1,000 (1001):   1/73
  
  # 0-5 (few), 50 (some), 500-1001 (many)
  fruit_py_cholla <- fruit_py_cholla %>%
    mutate(abund = case_when(
      max_count %in% 0:5 ~ "few",
      max_count == 50 ~ "some",
      max_count >= 500 ~ "many",
    )) %>%
    mutate(abund = factor(abund, levels = c("few", "some", "many")))
  
  # Visualize for each year
  fruit_pya_cholla <- fruit_py_cholla %>%
    mutate(fyr = factor(yr)) %>%
    mutate(site = factor(site_name)) %>%
    group_by(common_name, fyr, site, abund) %>%
    summarize(n = n(), .groups = "keep") %>%
    data.frame()
  ggplot(fruit_pya_cholla, aes(y = abund, x = fyr)) +
    geom_point(aes(size = n)) +
    facet_grid(site ~ .) +
    scale_y_discrete(labels = c("0-5", "50", "500+")) +
    scale_size_continuous(breaks = 1:3, range = c(2, 4)) +
    labs(x = "", y = "Abundance category", size = "No. plants", 
         title = "buck-horn cholla, fruit")
  
  fruit_py_cholla$site <- factor(fruit_py_cholla$site_name)
  fruit_py_cholla$fyr <- factor(fruit_py_cholla$yr)
  
  # Models
  m1 <- polr(abund ~ site, Hess = TRUE, data = fruit_py_cholla)
    summary(m1)
    # p-values (really only valid for larger sample sizes)
    ctable <- coef(summary(m1))
    p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
    (ctable <- cbind(ctable, "p value" = round(p, 3)))

  m2 <- polr(abund ~ fyr, Hess = TRUE, data = fruit_py_cholla)
  summary(m2)
  m3 <- polr(abund ~ site + fyr, Hess = TRUE, data = fruit_py_cholla)
  summary(m3)
  m4 <- polr(abund ~ site * fyr, Hess = TRUE, data = fruit_py_cholla)
  summary(m4)
  # Warning about rank-deficient design
  
  AIC(m1, m2, m3)
  # Additive model is best...
  
  # Strong evidence that max flower counts differed among years (though all max 
  # counts in 2019 were 500+ (which causes some estimation problems). 
  # Counts in 2021 all low. Counts don't vary that much among sites

  # saguaro -------------------------------------------------------------------#
  # Aggregate data for each plant-year: buck-horn cholla
  fruit_py_sag <- fruit %>%
    filter(common_name == "saguaro") %>%
    group_by(common_name, id, site_name, lat, lon, yr) %>%
    summarize(n_obs = n(),
              n_inphase = sum(status),
              max_count = max(intensity_midpoint),
              .groups = "keep") %>%
    data.frame()
  fruit_py_sag
  count(filter(gamdf, common_name == "saguaro" & 
                 intensity_short == "Fruits"),
        intensity_midpoint, intensity_value)
  count(fruit_py_sag, max_count)
  # 0:                        6/63
  # Less than 3 (1):          2/63
  # 3 to 10 (5):             10/63
  # 11 to 100 (50):          17/63
  # 101 to 1000 (500):       27/63
  # More than 1,000 (1001):   1/63
  
  # 0-5 (few), 50 (some), 500-1001 (many)
  fruit_py_sag <- fruit_py_sag %>%
    mutate(abund = case_when(
      max_count %in% 0:5 ~ "few",
      max_count == 50 ~ "some",
      max_count >= 500 ~ "many",
    )) %>%
    mutate(abund = factor(abund, levels = c("few", "some", "many")))
  
  # Visualize for each year
  fruit_pya_sag <- fruit_py_sag %>%
    mutate(fyr = factor(yr)) %>%
    mutate(site = factor(site_name)) %>%
    group_by(common_name, fyr, site, abund) %>%
    summarize(n = n(), .groups = "keep") %>%
    data.frame()
  ggplot(fruit_pya_sag, aes(y = abund, x = fyr)) +
    geom_point(aes(size = n)) +
    facet_grid(site ~ .) +
    scale_y_discrete(labels = c("0-5", "50", "500+")) +
    scale_size_continuous(breaks = 1:3, range = c(2, 4)) +
    labs(x = "", y = "Abundance category", size = "No. plants", 
         title = "saguaro, fruit")
  
  fruit_py_sag$site <- factor(fruit_py_sag$site_name)
  fruit_py_sag$fyr <- factor(fruit_py_sag$yr)
  
  # Models
  m1 <- polr(abund ~ site, Hess = TRUE, data = fruit_py_sag)
    summary(m1)
    # p-values (really only valid for larger sample sizes)
    ctable <- coef(summary(m1))
    p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
    (ctable <- cbind(ctable, "p value" = round(p, 3)))
  
  m2 <- polr(abund ~ fyr, Hess = TRUE, data = fruit_py_sag)
  summary(m2)
  m3 <- polr(abund ~ site + fyr, Hess = TRUE, data = fruit_py_sag)
  summary(m3)
  m4 <- polr(abund ~ site * fyr, Hess = TRUE, data = fruit_py_sag)
  summary(m4)
  # Warning about rank-deficient design
  
  AIC(m1, m2, m3)
  # Additive model is best...
  
  # Similar to flower counts, there are big site differences (counts higher at 
  # Brown's Ranch). Some annual differences, but estimation issues for 2018
  # since all counts were 500+

  # jojoba --------------------------------------------------------------------#
  # Aggregate data for each plant-year: buck-horn cholla
  fruit_py_jojo <- fruit %>%
    filter(common_name == "jojoba") %>%
    group_by(common_name, id, site_name, lat, lon, yr) %>%
    summarize(n_obs = n(),
              n_inphase = sum(status),
              max_count = max(intensity_midpoint),
              .groups = "keep") %>%
    data.frame()
  fruit_py_jojo
  count(filter(gamdf, common_name == "jojoba" & 
                 intensity_short == "Fruits"),
        intensity_midpoint, intensity_value)
  count(fruit_py_jojo, max_count)
  # 0:                       20/63
  # Less than 3 (1):          3/63
  # 3 to 10 (5):              2/63
  # 11 to 100 (50):          14/63
  # 101 to 1000 (500):       13/63
  # 1001 to 10000 (5000):    11/63
  # More than 10000 (10001):  0/63
  
  # Try something different here since there are a lot of 0s
  # 0 (none), 1-50 (few), 500 (some), 5000-10001 (many)
  fruit_py_jojo <- fruit_py_jojo %>%
    mutate(abund = case_when(
      max_count == 0 ~ "none",
      max_count %in% 1:50 ~ "few",
      max_count == 500 ~ "some",
      max_count >= 5000 ~ "many",
    )) %>%
    mutate(abund = factor(abund, levels = c("none", "few", "some", "many")))

  # Visualize for each year
  fruit_pya_jojo <- fruit_py_jojo %>%
    mutate(fyr = factor(yr)) %>%
    mutate(site = factor(site_name)) %>%
    group_by(common_name, fyr, site, abund) %>%
    summarize(n = n(), .groups = "keep") %>%
    data.frame()
  ggplot(fruit_pya_jojo, aes(y = abund, x = fyr)) +
    geom_point(aes(size = n)) +
    facet_grid(site ~ .) +
    scale_y_discrete(labels = c("0", "1-50", "500", "5000+")) +
    scale_size_continuous(breaks = 1:3, range = c(2, 4)) +
    labs(x = "", y = "Abundance category", size = "No. plants", 
         title = "jojoba, fruit")
  
  fruit_py_jojo$site <- factor(fruit_py_jojo$site_name)
  fruit_py_jojo$fyr <- factor(fruit_py_jojo$yr)
  
  # Models
  m1 <- polr(abund ~ site, Hess = TRUE, data = fruit_py_jojo)
    summary(m1)
    # p-values (really only valid for larger sample sizes)
    ctable <- coef(summary(m1))
    p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
    (ctable <- cbind(ctable, "p value" = round(p, 3)))
  
  m2 <- polr(abund ~ fyr, Hess = TRUE, data = fruit_py_jojo)
  summary(m2)
  m3 <- polr(abund ~ site + fyr, Hess = TRUE, data = fruit_py_jojo)
  summary(m3)
  m4 <- polr(abund ~ site * fyr, Hess = TRUE, data = fruit_py_jojo)
  summary(m4)
  # Warning about rank-deficient design
  
  AIC(m1, m2, m3)
  # Year model is best
  
  # Counts in 2018 and 2021 were a little lower, counts in 2022 are a bit 
  # higher. In last few years, counts at Brown's ranch were 5000+ but counts at 
  # other sites were much lower.
  
# -----------------------------------------------------------------------------#
# Exploring variation in max counts for fruit drop phenophase -----------------#
  
drop <- gamdf %>%
  filter(intensity_short == "FruitDrop") %>%
  dplyr::select(common_name, individual_id, site_name, latitude, longitude, 
                phenophase_description, observation_date, yr, day_of_year,
                phenophase_status, intensity_label, 
                intensity_midpoint) %>%
  rename(status = phenophase_status,
         phenophase = phenophase_description,
         id = individual_id,
         lat = latitude,
         lon = longitude,
         obsdate = observation_date, 
         doy = day_of_year)

# Distribution of intensity values
count(drop, intensity_midpoint) %>%
  mutate(prop = n / sum(n))
# A lot more 0s than flowers/fruit phenophases

count(filter(si, intensity_label == "No. fruit/seed drop"), 
      common_name, intensity_name, intensity_value)
# Might be able to model multiple species together for this phenophase since
# the "More than X" intensity category was very rarely used.

  # Trying to model mutliple species together ---------------------------------#
  # Keeping it to the 3 species I used for flowers, fruit
  # Aggregate data for each plant-year and species
  drop_py <- drop %>%
    filter(common_name %in% c("buck-horn cholla", "saguaro", "jojoba")) %>%
    group_by(common_name, id, site_name, lat, lon, yr) %>%
    summarize(n_obs = n(),
              n_inphase = sum(status),
              max_count = max(intensity_midpoint),
              .groups = "keep") %>%
    data.frame()
  count(filter(gamdf, intensity_short == "FruitDrop"),
        intensity_midpoint, intensity_value)

  count(drop_py, max_count)
  # 0:                      63/179
  # Less than 3 (1):        34/179
  # 3 to 10 (5):            31/179
  # 11 to 100 (50):         42/179
  # 101 to 1000 (500):       6/179
  # More than 1000 (1001):   2/179
  # More than 10000 (10001): 1/179
  
  # Trying something, but not sure this is worth doing. I don't think there's 
  # any justification for splitting 1 and 5, and I'm not sure there are enough
  # observations to have the 500-10001 category....
  # 0 (none); 1-5 (few), 50 (some), 500-10001 (many)
  drop_py <- drop_py %>%
    mutate(abund = case_when(
      max_count == 0 ~ "none",
      max_count %in% 1:5 ~ "few",
      max_count == 50 ~ "some",
      max_count > 50 ~ "many",
    )) %>%
    mutate(abund = factor(abund, levels = c("none", "few", "some", "many")))
  
  # Visualize for each year
  drop_pya <- drop_py %>%
    mutate(fyr = factor(yr)) %>%
    mutate(site = factor(site_name)) %>%
    mutate(spp = factor(common_name)) %>%
    group_by(spp, fyr, site, abund) %>%
    summarize(n = n(), .groups = "keep") %>%
    data.frame()
  ggplot(drop_pya, aes(y = abund, x = fyr)) +
    geom_point(aes(size = n)) +
    facet_grid(site ~ spp) +
    scale_y_discrete(labels = c("0", "1-5", "50", "500+")) +
    scale_size_continuous(breaks = 1:3, range = c(2, 4)) +
    labs(x = "", y = "Abundance category", size = "No. plants", 
         title = "3 species, fruit drop")
  
  drop_py$site <- factor(drop_py$site_name)
  drop_py$fyr <- factor(drop_py$yr)
  drop_py$spp <- factor(drop_py$common_name)
  
  # Models
  m_null <- polr(abund ~ 1, Hess = TRUE, data = drop_py)
  m_site <- polr(abund ~ site, Hess = TRUE, data = drop_py)
  summary(m_site)
  m_spp <- polr(abund ~ spp, Hess = TRUE, data = drop_py)
  summary(m_spp)
  m_fyr <- polr(abund ~ fyr, Hess = TRUE, data = drop_py)
  summary(m_fyr)
  m_sitesppfyr <- polr(abund ~ site + spp + fyr, Hess = TRUE, data = drop_py)
  summary(m_sitesppfyr)

  AIC(m_null, m_site, m_fyr, m_spp, m_sitesppfyr)
  # Full additive model very slightly better than year model, but even that is
  # <5 points different than null model so no strong patterns.
  