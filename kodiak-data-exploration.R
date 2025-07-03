# Exploring intensity data from Kodiak NWR
# ER Zylstra

library(dplyr)
library(stringr)
library(lubridate)
library(ggplot2)
library(cowplot)
library(mgcv)
library(marginaleffects)
library(glmmTMB)


# Plant phenophase classes
leaf_classes <- 1:5
flower_classes <- 6:9
fruit_classes <- 10:13

# List files with formatted intensity data from Kodiak
intensity_files <- list.files("npn-data",
                              pattern = "intensity-siteKNWR",
                              full.names = TRUE)

# Load and format data --------------------------------------------------------#

sites <- str_sub(intensity_files, 24, 29)

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
  select(-n_obs)

# Doublecheck that there's only one observation of each plant-phenophase per day
if (nrow(si) != nrow(distinct(si, individual_id, phenophase_id, observation_date))) {
  warning("There is more than one observation of some plant phenophases at ",
          sitecap, " in a day")
}

# Are there observations where the plant is out of phase but an intensity value
# is reported? If so, change the status but add a column to note that the data
# were amended.
si <- si %>%
  mutate(amended_status = ifelse(phenophase_status == 0 & !is.na(intensity_midpoint), 
                                 1, 0)) %>%
  mutate(phenophase_status = ifelse(phenophase_status == 0 & !is.na(intensity_midpoint),
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
    select(-c(interval_raw, same_ind, same_php, same_yr))
  
# si %>% count(interval) %>% mutate(prop = n / 40686)
# Vast majority (98%) of intervals are 3 days

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
  # remove plant-php-year combos with no observations in phase
  # remove plant-php-year combos with no observations with intensity values
  # remove plant-php-year combos with < 5 observations
  # remove plant-php-year combos with maximum interval >21 days
pl_ph_yr <- pl_ph_yr %>%
  mutate(remove = case_when(
    n_inphase == 0 ~ 1,
    n_intvalue == 0 ~ 1,
    nobs < 5 ~ 1,
    max_int > 21 ~ 1,
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
  select(-c(values, valuesq)) %>%
  data.frame()

# Create short name for intensity categories
intensity_cats <- intensity_cats %>%
  mutate(intensity_short = case_when(
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
            prop_00 = round(sum(first_status == 0 & last_status == 0) / n_plantyrs, 2),
            .groups = "keep") %>%
  data.frame()
# write.table(select(php_summary, -class_id), 
#             "clipboard", sep = "\t", row.names = FALSE)

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
    png_name <- paste0("output/kodiak-intensities-by-spp-php-yr/IntensityData-",
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
        # geom_point(data = agg,
        #            aes(x = day_of_year, y = yaxis, size = n_indiv),
        #            shape = 16, alpha = 0.4) +
        # scale_size_continuous(breaks = 1:4) +
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

# Weird pattern observed for one of the blueberry plants in 2019 --------------#

# Look at blueberry data (Road 1) in 2019
blue19 <- si %>%
  filter(common_name == "oval-leaf blueberry" & yr == 2019) %>%
  select(site_name, observation_date, day_of_year, interval, 
         class_id, phenophase_status, 
         intensity_label, intensity_type, intensity_midpoint, 
         intensity_midpoint_log) %>%
  mutate(name = ifelse(site_name == "Blueberry-Road System 1", 
                       "Road1", "Road2")) %>%
  mutate(name = factor(name)) %>%
  mutate(y = ifelse(intensity_type == "number", 
                    intensity_midpoint_log, intensity_midpoint))

ggplot(blue19, aes(x = day_of_year, y = y)) +
  geom_line(aes(color = name)) +
  facet_wrap(~intensity_label, scales = "free_y")
ggplot(filter(blue19, class_id %in% fruit_classes), 
       aes(x = day_of_year, y = y)) +
  geom_line(aes(color = name)) +
  facet_wrap(~intensity_label, scales = "free_y", ncol = 1)

blue19 %>%
  filter(name == "Road1" & day_of_year %in% 160:244) %>%
  select(intensity_label, observation_date, day_of_year, phenophase_status) %>%
  arrange(observation_date)

# Blueberry plant (Road 1) data in 2019 probably has some errors. Two periods
# in June (days 160-175) and August (days 232-244) where status of all 
# phenophases was zero, including leaves, fruits, and ripe fruits, that were
# all in phase (and often had very high intensity values) before and after.

# Also, there I'm not sure which observations of fruit drop are valid and which
# aren't, since there's a super high intensity value on day 223, but nothing
# before or after.

# Remove these questionable observations for now.
si <- si %>%
  filter(!(site_name == "Blueberry-Road System 1" &
             yr == 2019 &
             day_of_year %in% c(160:175, 232:244) &
             phenophase_status == 0)) %>%
  filter(!(site_name == "Blueberry-Road System 1" &
             yr == 2019 &
             day_of_year == 223 &
             intensity_label == "No. fruit/seed drop"))

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
            .groups = "keep") %>%
  data.frame()

# Create dataframe to hold estimates from GAM models:
combos_u <- combos %>%
  select(class_id, intensity_label, intensity_type, common_name, n_plants) %>%
  mutate(yr = NA)
gams <- gamdf %>%
  group_by(class_id, intensity_label, intensity_type, common_name, yr) %>%
  summarize(n_plants = n_distinct(individual_id),
            .groups = "keep") %>%
  data.frame() %>%
  rbind(combos_u) %>%
  arrange(class_id, common_name, yr) %>%
  mutate(k = NA,
         smooth = NA,
         model = NA,
         fyrs_p = NA,
         aic = NA,
         devexpl = NA,
         kcheck = NA,
         peak_date = NA, 
         peak_value = NA,
         convergence = NA)
# write.table(select(combos, -class_id), "clipboard", 
#             sep = "\t", row.names = FALSE)

# No. fruits ------------------------------------------------------------------#

int_label <- "No. fruits"
spps <- sort(unique(si$common_name))

# Extract data for that intensity category
counts_i <- filter(gamdf, intensity_label == int_label)

# Extract data for one species
spp <- spps[3]
intdf <- filter(counts_i, common_name == spp)

  message("Running models for ", spp, ": ", intdf$intensity_label[1])
  
  # Identify the rows of gams dataframe associated with species & intensity
  row1 <- which(gams$intensity_label == intdf$intensity_label[1] &
                  gams$common_name == spp & is.na(gams$yr))
  rowyrs <- which(gams$intensity_label == intdf$intensity_label[1] &
                    gams$common_name == spp & !is.na(gams$yr))

  # Make variables into factors
  intdf$fyr <- factor(intdf$yr)
  intdf$individual_id <- factor(intdf$individual_id)
  
  # Nudge intensity midpoint 0 values
  intdf <- intdf %>%
    mutate(midpoint_adj = case_when(
      intensity_midpoint == 0 ~ 0.01,
      .default = intensity_midpoint
    ))  
  
  # Aggregate data across individual plants within a year
  counts_yragg <- intdf %>%
    group_by(fyr, yr, day_of_year, intensity_midpoint, midpoint_adj) %>%
    summarize(n_indiv = n(), .groups = "keep") %>%
    data.frame()
  
  # From help for glmmTMB: If your response variable is defined on the 
  # closed interval [a,b], transform it to [0,1] via y_scaled <- (y-a)/(b-a)
  # Don't know whether this will work or not....
  
  intdf <- intdf %>%
    mutate(y_scaled = intensity_midpoint / max(intensity_midpoint))
  counts_yragg <- counts_yragg %>%
    mutate(y_scaled = intensity_midpoint / max(intensity_midpoint))
  # Tried to use y_scaled in lots of different forms, but nothing was working 
  # consistently (see bottom of script)
  
  
  
  library(ordbetareg)
  mobr <- ordbetareg(y_scaled ~ s(day_of_year),
                     data = filter(intdf, yr == 2015),
                     true_bounds = c(0, 1))
  plot_predictions(mobr, by = "day_of_year")
  
  mobr2 <- ordbetareg(intensity_midpoint ~ s(day_of_year, by = fyr),
                      data = intdf,
                      true_bounds = c(0, 5000))
  mobr2
  plot_predictions(mobr2, by = c("day_of_year", "fyr", "fyr"))
  # Getting warnings about divergent transitions, but on the whole, these
  # don't look bad.
  
  mobr3 <- ordbetareg(intensity_midpoint ~ s(day_of_year, by = fyr) + (1|individual_id),
                      data = intdf,
                      true_bounds = c(0, 5000),
                      cores = 2, chains = 2, 
                      control=list(adapt_delta=0.95))
  mobr3
  plot_predictions(mobr3, by = c("day_of_year", "fyr", "fyr"))
  # Come back to this. Not showing curves correctly in years with more than
  # one individual....
  
  # Probably want to specify cores and chains
  # Also, want to see about adding a random effect
  # Need to look more into how visualize/summarize results 
  # (though for now, I can at least use marginaleffects::plot_predictions)
  
  
  
  
  
  # Trying glmmTMB with ordbeta... --------------------------------------------#  
  mtest <- glmmTMB(y_scaled ~ s(day_of_year, fyr, bs = "fs"), 
                   data = intdf, 
                   REML = TRUE, family = ordbeta)
  plot_predictions(mtest, by = c("day_of_year", "fyr"), vcov = FALSE)
  
  myr <- glmmTMB(y_scaled ~ s(day_of_year, k = 30),
                 weights = n_indiv,
                 data = filter(counts_yragg, yr == 2022), 
                 REML = TRUE, family = ordbeta)
  plot_predictions(myr, by = c("day_of_year"), vcov = FALSE)
  
  glmmTMB(prop ~ s(day_of_year, k = kval, bs = cubic_bs) +
            s(day_of_year, fyr, k = kval, bs = "fs", xt = list(bs = cubic_bs), m = 1), 
          weights = n_indiv/mean(n_indiv),
          data = props_yragg, REML = TRUE, family = ordbeta)
  
  
  # Extract and save model settings
  gams$k[c(rowyrs, row1)] <- kval
  # gams$smooth[c(rowyrs, row1)] <- cubic_bs
  
  
  
  # Model with annual and global smooths ("GS" model in Pedersen et al. 2019).
  # Wrapping call in suppressWarnings since these models always spit out a 
  # warning about repeated 1-d smooths that we can safely ignore (according
  # to Gavin Simpson)
  suppressWarnings(
    m_gs <- gam(midpoint_adj ~ s(day_of_year, k = kval) +
                  s(day_of_year, fyr, k = kval, bs = "fs", m = 1), 
                weights = n_indiv,
                data = counts_yragg, method = "REML", 
                family = Gamma(link = "log"), select = TRUE)
  )
  summary(m_gs)
  plot(m_gs, pages = 1)
  plot_predictions(m_gs, by = c("day_of_year", "fyr"))
  
  
  
  
  
  # Extract and save model results
  gams$model[rowyrs] <- "GS"
  gams$fyrs_p[rowyrs] <- summary(m_gs)$s.table[2,"p-value"]
  gams$aic[rowyrs] <- round(AIC(m_gs), 2)
  gams$devexpl[rowyrs] <- round(summary(m_gs)$dev.expl * 100, 1)
  gams$kcheck[rowyrs] <- ifelse(any(k.check(m_gs)[, "p-value"] < 0.10),
                                "problem", "ok")
  
  # Plot global smooth:
  # plot(m_gs, select = 1)
  # plot(m_gs, select = 1, trans = exp, shift = coef(m_gs)[1], 
  #      rug = FALSE, se = FALSE)
  # plot_predictions(m_gs, by = "day_of_year", exclude = "s(day_of_year,fyr)", 
  #                  type = "link", transform = exp) 
  # Note that the CI around the global effect is huge. Similar to: 
  # https://stats.stackexchange.com/questions/645096/interpretation-help-of-summary-from-basic-gam-models-with-random-smooths 
  # If we reduce both k's a lot, then the CI gets smaller, but the annual 
  # smooths don't seem like they fit the data well at all and the gam.check 
  # p-values are much smaller.
  
  # Using marginaleffects package to look at annual smooths. 
  # (need to use link and transform arguments to calculate properly and 
  # prevent CI's from extending below 0.)
  # plot_predictions(m_gs, by = c("day_of_year", "fyr"), type = "link", transform = exp)
  # plot_predictions(m_gs, by = c("day_of_year", "fyr", "fyr"), type = "link", transform = exp)
  # Note: can get dataframe with estimates by including draw = FALSE for 
  # original data or newdata
  
  # Model with annual smooths and NO global smooth ("S" model in Pedersen)
  # m_s <- gam(midpoint_adj ~ s(day_of_year, fyr, k = kval, bs = "fs",
  #                             xt = list(bs = cubic_bs)),
  #            weights = n_indiv/mean(n_indiv),
  #            data = counts_yragg, method = "REML",
  #            family = Gamma(link = "log"), select = TRUE)
  # summary(m_s)
  # plot_predictions(m_s, by = c("day_of_year", "fyr"))
  
  # Model with only a global smooth ("G" model in Pedersen; but keeping a 
  # yearly random effect in the model).
  m_g <- gam(midpoint_adj ~ s(day_of_year, k = kval, bs = cubic_bs) + s(fyr, bs = "re"), 
             weights = n_indiv/mean(n_indiv),
             data = counts_yragg, method = "REML", 
             family = Gamma(link = "log"), select = TRUE)
  # summary(m_g)
  # plot_predictions(m_g, by = c("day_of_year"))
  # gam.check(m_g)
  
  # Extract and save model results
  gams$model[row1] <- "G"
  gams$aic[row1] <- round(AIC(m_g), 2)
  gams$devexpl[row1] <- round(summary(m_g)$dev.expl * 100, 1)
  gamcheck_pvalues <- k.check(m_g)[, "p-value"]
  gamcheck_pvalues <- gamcheck_pvalues[!is.na(gamcheck_pvalues)]
  gams$kcheck[row1] <- ifelse(any(gamcheck_pvalues < 0.10), "problem", "ok")
  
  if (gams$kcheck[row1] == "problem") {
    warning("gam.check for ", spp, ":", counts$intensity_label[i],
            " (G) indicates potential problems with model fit.")
  }
  if (any(gams$kcheck[rowyrs] == "problem")) {
    warning("gam.check for ", spp, ":", counts$intensity_label[i],
            " (GS) indicates potential problems with model fit.")
  }
  

