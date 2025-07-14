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

# List files with formatted intensity data from Kodiak
intensity_files <- list.files("npn-data",
                              pattern = "intensity-siteMCDO",
                              full.names = TRUE)

# Load and format data --------------------------------------------------------#

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

# Summarize amount and quality of information for each plant, phenophase, year
pl_ph_yr <- si %>%
  group_by(common_name, individual_id, phenophase_description, yr) %>%
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

# Filter data: plant-phenophase-year combinations -----------------------------#

# Want to:
# remove any pollen release data (qualitative data)
# remove plant-php-year combos with < 5 observations
# remove plant-php-year combos with no observations in phase
# remove plant-php-year combos with no observations with intensity values
# remove plant-php-year combos with maximum interval >21 days #### Might be worth revisiting ####
pl_ph_yr <- pl_ph_yr %>%
  mutate(remove = case_when(
    phenophase_description == "Pollen release (flowers)" ~ 1,
    nobs < 5 ~ 1,
    n_inphase == 0 ~ 1,
    n_intvalue == 0 ~ 1,
    max_int > 21 ~ 1,
    .default = 0
  ))
si <- si %>%
  left_join(dplyr::select(pl_ph_yr, individual_id, phenophase_description, yr, remove),
            by = c("individual_id", "phenophase_description", "yr")) %>%
  filter(remove == 0) %>%
  dplyr::select(-remove)

# Check that there's only one intensity category for each species-phenophase?
spil <- si %>%
  filter(!is.na(intensity_category_id)) %>%
  distinct(common_name, phenophase_description, intensity_name,
           intensity_type, intensity_label)
spl <- si %>%
  distinct(common_name, phenophase_description)
# setdiff(spl, spil[, 1:2])

# Make intensity midpoints = 0 if status = 0 and add intensity category labels
si <- si %>%
  mutate(intensity_midpoint = ifelse(phenophase_status == 0, 
                                     0, intensity_midpoint)) %>%
  dplyr::select(-c(intensity_name, intensity_type, intensity_label)) %>%
  left_join(spil, by = c("common_name", "phenophase_description"))

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
    intensity_label == "Ripe fruit (%)" ~ "Ripe fruit",
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
  select(-n_yrs) %>%
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
    intensity_label == "Ripe fruit (%)" ~ "Ripe fruit",
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
  mutate(prop = n / sum(n))

# Would like to model data for multiple species simultaneously, but they have 
# different intensity categories, so I'd have to create bins that work for all. 
# Max category == "More than 1,000": Barrel cacus, cassia, cholla, saguaro
# Max category == "More than 10,000": jojoba, soaptree yucca, velvet mesquite

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
  # 0:                        0/69
  # Less than 3 (1):          5/69
  # 3 to 10 (5):              2/69
  # 11 to 100 (50):          27/69
  # 101 to 1000 (500):       31/69
  # More than 1,000 (1001):   4/69
  
  # Not sure about grouping, but for now will try:
  # 1-5 (few), 50 (some), 500-1001 (many) [Don't need "none" group]
  flowers_py_cholla <- flowers_py_cholla %>%
    mutate(abund = case_when(
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
    scale_y_discrete(labels = c("1-5", "50", "500+")) +
    scale_size_continuous(breaks = 1:3) +
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
  
  AIC(m1, m2, m3, m4)
  # Additive model is the best...
  
  # Strong evidence that max flower counts differed among years, mostly because
  # counts in 2021 were much lower at all sites. Some evidence that flower counts
  # wwere lower on average at the Gateway site than other two. 
  
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
  # 0:                        0/63
  # Less than 3 (1):          1/63
  # 3 to 10 (5):              8/63
  # 11 to 100 (50):          20/63
  # 101 to 1000 (500):       27/63
  # More than 1,000 (1001):   7/63
  
  # Will try same grouping as buck-horn cholla:
  # 1-5 (few), 50 (some), 500-1001 (many) [Don't need "none" group]
  flowers_py_sag <- flowers_py_sag %>%
    mutate(abund = case_when(
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
    scale_y_discrete(labels = c("1-5", "50", "500+")) +
    scale_size_continuous(breaks = 1:3) +
    labs(x = "", y = "Abundance category", size = "No. plants", 
         title = "saguaro, flowers")
  
  # Test run with MASS::polr
  flowers_py_sag$site <- factor(flowers_py_sag$site_name)
  flowers_py_sag$fyr <- factor(flowers_py_sag$yr)
  
  # By default, logistic regression (could change to probit)
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
  
  m4 <- polr(abund ~ site * fyr, Hess = TRUE, data = flowers_py_sag)
  summary(m4)
  
  AIC(m1, m2, m3, m4)
  # Additive model is the best...
  
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
  # 1-50 (few), 500 (some), 5000-10001 (many) [Don't need "none" group]
  flowers_py_jojo <- flowers_py_jojo %>%
    mutate(abund = case_when(
      max_count %in% 1:50 ~ "few",
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
    scale_y_discrete(labels = c("1-50", "500", "5000+")) +
    scale_size_continuous(breaks = 1:3) +
    labs(x = "", y = "Abundance category", size = "No. plants", 
         title = "jojoba, flowers")
  
  # Test run with MASS::polr
  flowers_py_jojo$site <- factor(flowers_py_jojo$site_name)
  flowers_py_jojo$fyr <- factor(flowers_py_jojo$yr)
  
  # By default, logistic regression (could change to probit)
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
  
  AIC(m1, m2, m3, m4)
  # Additive model is the best...
  
  # Strong evidence that max flower counts differed among sites, with lower 
  # number of flowers at Lost Dog. Does look like some annual variation, but
  # there are some estimation difficulties because all plants monitored in 2017
  # had 5000+ max counts.
  
# -----------------------------------------------------------------------------#
# Exploring variation in max counts for flowers phenophase --------------------#
  
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
  # 0:                        0/71
  # Less than 3 (1):          4/71
  # 3 to 10 (5):              6/71
  # 11 to 100 (50):          27/71
  # 101 to 1000 (500):       33/71
  # More than 1,000 (1001):   1/71
  
  # Not sure about grouping, but for now will try:
  # 1-5 (few), 50 (some), 500-1001 (many) [Don't need "none" group]
  fruit_py_cholla <- fruit_py_cholla %>%
    mutate(abund = case_when(
      max_count == 1 ~ "few",
      max_count == 5 ~ "few",
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
    scale_y_discrete(labels = c("1-5", "50", "500+")) +
    scale_size_continuous(breaks = 1:3) +
    labs(x = "", y = "Abundance category", size = "No. plants", 
         title = "buck-horn cholla, fruit")
  
  # Test run with MASS::polr
  fruit_py_cholla$site <- factor(fruit_py_cholla$site_name)
  fruit_py_cholla$fyr <- factor(fruit_py_cholla$yr)
  
  # By default, logistic regression (could change to probit)
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
  
  AIC(m1, m2, m3, m4)
  # Additive model is the best...
  
  # Strong evidence that max flower counts differed among years, similarly
  # among sites. All max counts in 2019 were 500+ (which causes some
  # estimation problems). Counts in 2021 all low.

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
  # 0:                        0/57
  # Less than 3 (1):          2/57
  # 3 to 10 (5):             10/57
  # 11 to 100 (50):          17/57
  # 101 to 1000 (500):       27/57
  # More than 1,000 (1001):   1/57
  
  # 1-5 (few), 50 (some), 500-1001 (many) [Don't need "none" group]
  fruit_py_sag <- fruit_py_sag %>%
    mutate(abund = case_when(
      max_count == 1 ~ "few",
      max_count == 5 ~ "few",
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
    scale_y_discrete(labels = c("1-5", "50", "500+")) +
    scale_size_continuous(breaks = 1:3) +
    labs(x = "", y = "Abundance category", size = "No. plants", 
         title = "saguaro, fruit")
  
  # Test run with MASS::polr
  fruit_py_sag$site <- factor(fruit_py_sag$site_name)
  fruit_py_sag$fyr <- factor(fruit_py_sag$yr)
  
  # By default, logistic regression (could change to probit)
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
  
  AIC(m1, m2, m3, m4)
  # Additive model is the best...
  
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
  # 0:                        0/43
  # Less than 3 (1):          3/43
  # 3 to 10 (5):              2/43
  # 11 to 100 (50):          14/43
  # 101 to 1000 (500):       13/43
  # 1001 to 10000 (5000):    11/43
  # More than 10000 (10001):  0/43
  
  # Use flowers grouping
  # 1-50 (few), 500 (some), 5000-10001 (many) [Don't need "none" group]
  fruit_py_jojo <- fruit_py_jojo %>%
    mutate(abund = case_when(
      max_count %in% 1:50 ~ "few",
      max_count == 500 ~ "some",
      max_count >= 5000 ~ "many",
    )) %>%
    mutate(abund = factor(abund, levels = c("few", "some", "many")))

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
    scale_y_discrete(labels = c("1-5", "50", "500+")) +
    scale_size_continuous(breaks = 1:3) +
    labs(x = "", y = "Abundance category", size = "No. plants", 
         title = "jojoba, fruit")
  
  # Test run with MASS::polr
  fruit_py_jojo$site <- factor(fruit_py_jojo$site_name)
  fruit_py_jojo$fyr <- factor(fruit_py_jojo$yr)
  
  # By default, logistic regression (could change to probit)
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
  # Warning about rank-deficient design (like flowers)
  
  AIC(m1, m2, m3)
  # Site model is the best...
  
  # Less data for jojoba fruit, especially in early years. I think this is why
  # the simplest (site) model has the lowest AIC. In last couple years, counts
  # at Brown's ranch were 500+ but counts at other sites were much lower.
  
  
  
    