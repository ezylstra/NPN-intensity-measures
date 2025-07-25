# Exploring intensity data
# ER Zylstra

library(dplyr)
library(stringr)
library(lubridate)
library(ggplot2)
library(cowplot)
library(mgcv)
library(marginaleffects)
library(glmmTMB)

# List of sites
sites <- c("bona", "deju", "srer", "ellen")

# Plant phenophase classes
leaf_classes <- 1:5
flower_classes <- 6:9
fruit_classes <- 10:13

# List files with formatted intensity data
intensity_files <- list.files("npn-data",
                              pattern = "intensity-site",
                              full.names = TRUE)

# Load and format data --------------------------------------------------------#

# Summarize data available by species, individual, phenophase
for (site in sites) {
  
  sitecap <- str_to_upper(site)
  filename <- intensity_files[grepl(sitecap, intensity_files)]
  si <- read.csv(filename)
  
  # Remove any true duplicates
  si <- si[!duplicated(si),]
  
  # Occasionally there are two records for a plant-phenophase in one day with 
  # different intensity or status values. Removing these observations too (not
  # trying to resolve conflicts, just deleting all observations).
  si <- si %>%
    group_by(individual_id, phenophase_id, observation_date) %>%
    mutate(n_obs = n()) %>%
    ungroup() %>%
    filter(n_obs == 1) %>%
    select(-n_obs) %>%
    data.frame()
  
  # Doublecheck that there's only one observation of each plant-phenophase per day
  if (nrow(si) != nrow(distinct(si, individual_id, phenophase_id, observation_date))) {
    warning("There is more than one observation of some plant phenophases at ",
            sitecap, " in a day")
  }
  
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
    select(-c(interval_raw, same_ind, same_php, same_yr))
  
  assign(paste0("si_", site), si)
}

# Create intensity dataset combined for all sites
si_files <- ls()[ls() %in% paste0("si_", sites)]
si <- do.call(rbind, mget(si_files))
rownames(si) <- NULL
rm(list = si_files)

# Add short site name
si <- si %>%
  mutate(site = str_to_lower(str_sub(site_name, 1, 4))) %>%
  mutate(site = str_replace_all(site, "home", "ellen"))

# Summarize amount and quality of information for each plant, phenophase, year
pl_ph_yr <- si %>%
  group_by(site, common_name, individual_id, phenophase_description, yr) %>%
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
  # remove plant-php-year combos with no observations in phase
  # remove plant-php-year combos with no observations with intensity values
  # remove plant-php-year combos with < 5 observations
  # remove plant-php-year combos when max interval > 14 days
  # (using 14-day max interval cutoff rather than mean interval doesn't remove
  # species/plants from dataset but reduces the number of years for some
  # plants at Ellen's site; most now 3-4 years)

# Side note: Checked that individual_id's for plants at NEON sites weren't 
# changing annually like the site_id's do. Luckily, that doesn't seem to be the 
# case, as most plants are monitored for multiple years.

pl_ph_yr <- pl_ph_yr %>%
  mutate(remove = case_when(
    n_inphase == 0 ~ 1,
    n_intvalue == 0 ~ 1,
    nobs < 5 ~ 1,
    # mean_int > 14 ~ 1,
    max_int > 14 ~ 1,
    .default = 0
  ))
si <- si %>%
  left_join(select(pl_ph_yr, individual_id, phenophase_description, yr, remove),
            by = c("individual_id", "phenophase_description", "yr")) %>%
  filter(remove == 0) %>%
  select(-remove)

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
  select(-c(intensity_name, intensity_type, intensity_label)) %>%
  left_join(spil, by = c("common_name", "phenophase_description"))

# Filter data: species --------------------------------------------------------#

# Site-species summaries
site_spp <- si %>%
  mutate(php_type = case_when(
    class_id %in% leaf_classes ~ "Leaves",
    class_id %in% flower_classes ~ "Flowers",
    class_id %in% fruit_classes ~ "Fruit"
  )) %>%
  group_by(site, common_name, species_functional_type) %>%
  summarize(n_plants = n_distinct(individual_id),
            first_yr = min(yr),
            last_yr = max(yr),
            n_yrs = n_distinct(yr),
            n_php = n_distinct(phenophase_id),
            leaf_php = n_distinct(phenophase_id[php_type == "Leaves"]),
            flower_php = n_distinct(phenophase_id[php_type == "Flowers"]),
            fruit_php = n_distinct(phenophase_id[php_type == "Fruit"]),
            .groups = "keep") %>%
  arrange(site, species_functional_type, desc(n_plants), desc(n_yrs),
          desc(n_php)) %>%
  data.frame()

# Identify which species to focus on at Ellen's site: 1/2 species per functional
# type, prioritizing no. plants, no. years, and then no. phenophases
site_spp <- site_spp %>%
  group_by(site, species_functional_type) %>%
  mutate(priority = row_number()) %>%
  # Keep all species that have at least 3 plants
  mutate(priority = ifelse(n_plants > 2, 1, priority)) %>%
  data.frame()

# Remove any species that are lower priority
si <- si %>%
  left_join(select(site_spp, site, common_name, priority),
            by = c("site", "common_name")) %>%
  filter(priority == 1) %>%
  select(-priority)
site_spp <- site_spp %>%
  filter(priority == 1) %>%
  select(-priority)

# Filter data: intensity categories -------------------------------------------#

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
  select(-c(values, valuesq)) %>%
  data.frame()

# Identify categories that we want to exclude because:
# Only one or two intensity values represented in dataset
# Fewer than 50 observations across species, plants, and years
intensity_cats <- intensity_cats %>%
  mutate(exclude = case_when(
    str_count(unique_values, ",") < 2 ~ 1,
    n < 50 ~ 1, 
    .default = 0
  ))

# Remove categories from dataset and intensity category table
si <- si %>%
  left_join(select(intensity_cats, intensity_label, exclude), 
            by = "intensity_label") %>%
  filter(exclude == 0) %>%
  select(-exclude)
intensity_cats <- intensity_cats %>%
  filter(exclude == 0) %>%
  select(-exclude)

# Create short name for intensity categories
intensity_cats <- intensity_cats %>%
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

# Plot raw intensity data -----------------------------------------------------#

# Create column with logged intensity values
si <- si %>%
  mutate(intensity_midpoint_log = ifelse(intensity_midpoint == 0,
                                         log(intensity_midpoint + 0.01), 
                                         log(intensity_midpoint)))

# Identify the min/max number of plants per species, intensity category, and 
# year to keep point size consistent among plots
nplants_size <- si %>%
  group_by(intensity_label, common_name, yr) %>%
  summarize(n_plants = n_distinct(individual_id), .groups = "keep") %>%
  data.frame()

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
    png_name <- paste0("output/intensities-by-spp-php-yr/IntensityData-",
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
        scale_size_continuous(limits = c(1, max(nplants_size$n_plants))) +
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

# NEON sites had some years where monitoring started late (but note that this 
# may not matter for phenophases that occur later in the year)
  # BONA [quaking aspen]: 2017 and 2020 
  # DEJU [dwarf birch]: 2016 and 2020
  # SRER [velvet mesquite, creosote bush, desert zinnia]: 2016
si <- si %>%
  filter(!(site == "bona" & yr %in% c(2017, 2020))) %>%
  filter(!(site == "deju" & yr %in% c(2016, 2020))) %>%
  filter(!(site == "srer" & yr == 2016))

# Use a rule to exclude any species-phenophase combination where there's only 
# one year of data. For now, will keep the few combinations where we have one 
# plant observed over multiple years.
si <- si %>%
  group_by(common_name, phenophase_description) %>%
  mutate(n_yrs = n_distinct(yr)) %>%
  ungroup() %>%
  filter(n_yrs > 1) %>%
  select(-n_yrs) %>%
  data.frame()

# Add short name for intensity categories
si <- si %>%
  left_join(select(intensity_cats, intensity_label, intensity_short),
            by = "intensity_label")

# Sort data and remove any observations with missing intensity values
gamdf <- si %>%
  arrange(class_id, phenophase_description, site, common_name, individual_id, 
          observation_date) %>%
  filter(!is.na(intensity_midpoint))

# Summarize filtered data -----------------------------------------------------#

# Site-species summary
site_spp <- gamdf %>%
  mutate(php_type = case_when(
    class_id %in% leaf_classes ~ "Leaves",
    class_id %in% flower_classes ~ "Flowers",
    class_id %in% fruit_classes ~ "Fruit"
  )) %>%
  group_by(site, common_name, species_functional_type) %>%
  summarize(n_plants = n_distinct(individual_id),
            first_yr = min(yr),
            last_yr = max(yr),
            n_yrs = n_distinct(yr),
            n_php = n_distinct(phenophase_id),
            leaf_php = n_distinct(phenophase_id[php_type == "Leaves"]),
            flower_php = n_distinct(phenophase_id[php_type == "Flowers"]),
            fruit_php = n_distinct(phenophase_id[php_type == "Fruit"]),
            .groups = "keep") %>%
  arrange(site, species_functional_type, desc(n_plants), desc(n_yrs),
          desc(n_php)) %>%
  mutate(site = str_to_upper(site),
         yr_range = paste(first_yr, "-", last_yr),
         phps = paste0(n_php, " (", leaf_php, ",", flower_php, 
                       ",", fruit_php, ")")) %>%
  relocate(yr_range, .before = "n_yrs") %>%
  select(-c(first_yr, last_yr, n_php, leaf_php, flower_php, fruit_php)) %>%
  data.frame()
# site_spp
# write.table(site_spp, "clipboard", sep ="\t", row.names = FALSE)

# Summarize data available for each plant, phenophase, and year
  # Going to use rle(phenophase_status) to understand the patterns:
  # length(rle$values) = number of 0/1 status sequences
  # first(rle$values) = state at first observation
  # last(rle$values) = state at last observation
pl_ph_yr <- gamdf %>%
  arrange(site, individual_id, phenophase_id, observation_date) %>%
  group_by(site, common_name, individual_id, phenophase_description,
           intensity_label, class_id, yr) %>%
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
            n_transitions_mn = n_states_mn - 1,
            n_transitions_max = n_states_max - 1,
            prop_00 = round(sum(first_status == 0 & last_status == 0) / n_plantyrs, 2),
            .groups = "keep") %>%
  data.frame() %>%
  select(-c(class_id, n_states_mn, n_states_max))
# php_summary
# write.table(php_summary, "clipboard", sep ="\t", row.names = FALSE)

# Summary for each species-intensity category combination
spp_int <- gamdf %>%
  mutate(plantyr = paste0(individual_id, "_", yr)) %>%
  group_by(class_id, intensity_label, intensity_type, common_name) %>%
  summarize(n_plantyrs = n_distinct(plantyr),
            n_plants = n_distinct(individual_id),
            n_yrs = n_distinct(yr),
            .groups = "keep") %>%
  data.frame()
# head(spp_int)
# write.table(spp_int, "clipboard", sep = "\t", row.names = FALSE)

# -----------------------------------------------------------------------------#
# Exploring variation in max annual counts ------------------------------------#

# Exploring with no. of young leaves on velvet mesquite
veme <- gamdf %>%
  filter(intensity_short == "YoungLeaves") %>%
  filter(common_name == "velvet mesquite") %>%
  dplyr::select(common_name, individual_id, latitude, longitude, 
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
count(veme, intensity_midpoint) %>%
  mutate(prop = n / sum(n))

# Aggregate data for each plant-year
veme_py <- veme %>%
  group_by(common_name, id, yr) %>%
  summarize(n_obs = n(),
            n_inphase = sum(status),
            max_count = max(intensity_midpoint),
            .groups = "keep") %>%
  data.frame()
veme_py
count(veme_py, max_count)
# 0:                     0/100
# Less than 3 (1):       0/100
# 3 to 10 (5):           0/100
# 11 to 100 (50):        1/100
# 101 to 1000 (500):    61/100
# 1001 to 10000 (5000): 38/100

# Not really worth looking into too much given that there isn't too much
# variation and there's only one site

# Leaf size -------------------------------------------------------------------#

leafs <- gamdf %>%
  filter(intensity_label == "Leaf size (%)") %>%
  mutate(prop = intensity_midpoint / 100)

(leafs_spp <- unique(leafs$common_name))
# quaking aspen (4 yrs; 3141 obs)
# dwarf birch (4 yrs; 2701 obs)
# forsythia (4 yrs; 625 obs)
# sugar maple (4 yrs; 527 obs)

# Theoretically, leaf size (for deciduous plants with a single flush of leaves) 
# is supposed to increase monotonically during phenophase. Does this seem to be
# the case?

leafs_byind <- leafs %>%
  group_by(common_name, yr, individual_id) %>%
  summarize(n_obs = n(),
            min_doy = min(day_of_year),
            max_doy = max(day_of_year),
            n_yes = sum(phenophase_status),
            first_yes = min(day_of_year[phenophase_status == 1]),
            last_yes = max(day_of_year[phenophase_status == 1]),
            # How many series of yeses (ideally, just one)
            n_yes_series = sum(rle(phenophase_status)$values == 1),
            vals = paste(rle(intensity_midpoint)$values[-length(rle(intensity_midpoint)$values)], collapse = ","),
            # Do intensity values (except for the last 0) monotonically increase?
            monotone_incr = all(as.integer(str_split(vals, ",")[[1]]) ==
                                  cummax(as.integer(str_split(vals, ",")[[1]]))),
            .groups = "keep") %>%
  select(-vals) %>%
  data.frame()

leafs_sppyr <- leafs_byind %>%
  group_by(common_name, yr) %>%
  summarize(n_indiv = n(),
            n_yesdates = round(mean(n_yes), 1),
            single_yes_series = sum(n_yes_series == 1),
            single_yes_monotone = sum(n_yes_series == 1 & monotone_incr == TRUE),
            .groups = "keep") %>%
  data.frame() %>%
  mutate(percent_single_yes = round(single_yes_series / n_indiv, 2),
         percent_single_mono = round(single_yes_monotone / n_indiv, 2))
leafs_sppyr

# Plots for one species
leafss <- filter(leafs, common_name == leafs_spp[1])
yrs <- unique(leafss$yr)
for (j in 1:length(yrs)) {
  leafss1 <- filter(leafss, yr == yrs[j])
  ls_plot <- ggplot(data = leafss1, aes(x = day_of_year, y = intensity_midpoint)) +
    geom_line() +
    facet_wrap(~factor(individual_id)) +
    labs(x = "Day of year", y = "% Leaf size", 
         title = paste0(str_to_sentence(leafss1$common_name[1]),
                        ", ", leafss1$yr[1])) +
    theme_bw()
  print(ls_plot)
}

# Run models, 
# Just for species in years when at least 50% of individuals had a single series of yeses
# Just when the mean number of dates observed in phase is > 5
leafs_combos <- leafs_sppyr %>%
  filter(percent_single_yes > 0.5 & n_yesdates > 5)

for (i in 1:nrow(leafs_combos)) {
  spp <- leafs_combos$common_name[i]
  sppyr <- leafs_combos$yr[i]
  leafs1 <- leafs %>%
    filter(common_name == spp & yr == sppyr) %>%
    filter(prop != 0) %>%
    mutate(individual_id = factor(individual_id))

  agg <- leafs1 %>%
    group_by(day_of_year, prop) %>%
    summarize(n_indiv = n_distinct(individual_id), .groups = "keep") %>%
    data.frame()
  
  m <- gam(prop ~ s(day_of_year, k = 5, bs = "cr") + s(individual_id, bs = "re"),
           data = leafs1, method = "REML", family = betar(link="logit"), select = TRUE)
  # summary(m)
  
  # plot_predictions(m, by = c("day_of_year", "individual_id"))
  # plot_predictions(m, by = "day_of_year") +
  #   geom_point(data = agg, aes(x = day_of_year, y = prop, size = n_indiv), color = "gray") +
  #   theme_bw()

  plot_dat <- data.frame(
    day_of_year = seq(min(leafs1$day_of_year), max(leafs1$day_of_year)),
    individual_id = leafs1$individual_id[1]
  )
  
  ilink <- family(m)$linkinv
  preds <- predict(m, newdata = plot_dat, type = "link", 
                   se.fit = TRUE, exclude = "s(individual_id)")
  preds <- cbind(plot_dat, preds)
  preds <- preds %>%
    mutate(lwr_ci = ilink(fit - (2 * se.fit)),
           upr_ci = ilink(fit + (2 * se.fit)),
           fitted = ilink(fit))
  
  pred_plot <- ggplot(data = preds, aes(x = day_of_year, y = fitted)) +
    geom_ribbon(aes(ymin = lwr_ci, ymax = upr_ci), alpha = 0.2) +
    # geom_line(data = leafs1, 
    #           aes(x = day_of_year, y = prop, color = individual_id),
    #           alpha = 0.3,
    #           show.legend = FALSE) +
    geom_point(data = agg,
               aes(x = day_of_year, y = prop, size = n_indiv),
               shape = 16, alpha = 0.7, color = "forestgreen") +
    geom_line() +
    labs(x = "Day of year",
         y = "Proportion of full leaf size",
         size = "No. plants",
         title = paste0(str_to_sentence(leafs1$common_name[1]),
                        ", ", leafs1$yr[1])) +
    theme_bw()
  
  print(pred_plot)
  
}

# Model multiple years at once? Try forythia, which has 4 years of data
forsyth <- leafs %>%
  filter(common_name == "forsythia") %>%
  filter(prop != 0) %>%
  mutate(individual_id = factor(individual_id)) %>%
  mutate(fyr = factor(yr))

agg_forsyth <- forsyth %>%
  group_by(day_of_year, prop) %>%
  summarize(n_obs = n(), .groups = "keep") %>%
  data.frame()

m2 <- gam(prop ~ s(day_of_year, k = 5, bs = "cr") + 
            s(day_of_year, fyr, bs = "fs", xt = list(bs = "cr")) +
            s(individual_id, bs = "re"),
          data = forsyth, method = "REML", family = betar(link="logit"), 
          select = TRUE)
summary(m2)

plot_dat2 <- data.frame(
  day_of_year = rep(seq(min(forsyth$day_of_year), max(forsyth$day_of_year)), 4),
  fyr = rep(c(2020, 2021, 2023, 2024), each = max(forsyth$day_of_year) - min(forsyth$day_of_year) + 1),
  individual_id = forsyth$individual_id[1]
)

ilink <- family(m2)$linkinv
preds <- predict(m2, newdata = plot_dat2, type = "link", 
                 se.fit = TRUE, exclude = "s(individual_id)")
preds <- cbind(plot_dat2, preds)
preds <- preds %>%
  mutate(lwr_ci = ilink(fit - (2 * se.fit)),
         upr_ci = ilink(fit + (2 * se.fit)),
         fitted = ilink(fit)) %>%
  mutate(fyr = factor(fyr))

ggplot(preds, aes(x = day_of_year, y = fitted)) +
  geom_ribbon(aes(ymin = lwr_ci, ymax = upr_ci, fill = fyr), alpha = 0.3) +
  geom_line(aes(color = fyr))

m3 <- gam(prop ~ s(day_of_year, k = 5, bs = "cr") + s(fyr, individual_id, bs = "re"),
          data = forsyth, method = "REML", family = betar(link="logit"), 
          select = TRUE)
summary(m3)

predsG <- predict(m3, newdata = filter(plot_dat2, fyr == 2020), type = "link", 
                 se.fit = TRUE, exclude = "s(fyr,individual_id)")
predsG <- cbind(filter(plot_dat2, fyr == 2020), predsG)
predsG <- predsG %>%
  mutate(lwr_ci = ilink(fit - (2 * se.fit)),
         upr_ci = ilink(fit + (2 * se.fit)),
         fitted = ilink(fit)) %>%
  mutate(fyr = factor(fyr))

ggplot(preds, aes(x = day_of_year, y = fitted)) +
  # geom_ribbon(aes(ymin = lwr_ci, ymax = upr_ci, fill = fyr), alpha = 0.3) +
  geom_line(aes(color = fyr)) +
  geom_ribbon(data = predsG, aes(ymin = lwr_ci, ymax = upr_ci), alpha = 0.3) +
  geom_line(data = predsG, aes(y = fitted)) +
  labs(x = "Day of year", y = "Proportion of full leaf size",
       color = "Year", title = "Forsythia, leaf size (%)") +
  theme_bw()
