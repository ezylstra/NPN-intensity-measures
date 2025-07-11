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
library(ordbetareg)


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
            prop_00 = round(sum(first_status == 0 & last_status == 0)/n_plantyrs, 2),
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

# Model exploration for no. fruits dataset ------------------------------------#

int_label <- "No. fruits"
spps <- sort(unique(si$common_name))

# Extract data for that intensity category
counts_i <- filter(gamdf, intensity_label == int_label)

# Extract data for one species
spp <- spps[3]
intdf <- filter(counts_i, common_name == spp)
intdf$fyr <- factor(intdf$yr)

# Nudge intensity midpoint 0 values
intdf <- intdf %>%
  mutate(midpoint_adj = case_when(
    intensity_midpoint == 0 ~ 0.01,
    .default = intensity_midpoint
  ))  

message("Running models for ", spp, ": ", intdf$intensity_label[1])

# Pick one year (year with 3 plants, different max values)?
intdf1 <- filter(intdf, yr == 2021)
intdf1$individual_id <- factor(intdf1$individual_id)

# GAM (Gamma family after nudging 0s) -----------------------------------------#

# Model with random effects for individual plants
# Note: need to have sufficiently large k or predicted mean values get silly high
mgam <- gam(midpoint_adj ~ s(day_of_year, bs = "cr", k = 20) + 
              s(individual_id, bs = "re"), 
            data = intdf1, method = "REML", 
            family = Gamma(link = "log"), select = TRUE)
  
summary(mgam)
gam.check(mgam)
plot(mgam, pages = 1)

# Plots -- why are they CIs so different (plot_predictions bigger)?
plot(mgam, select = 1, trans = exp, shift = coef(mgam)[1],
     rug = FALSE, se = TRUE, seWithMean = TRUE, unconditional = TRUE)
plot_predictions(mgam, by = "day_of_year", exclude = "s(individual_id)",
                 type = "link", transform = exp, points = 0.5)

# Predictions by hand (CIs match up with plot_predictions, not plot.gam)
ilink <- family(mgam)$linkinv
newdata_g <- data.frame(
  day_of_year = seq(min(intdf1$day_of_year), max(intdf1$day_of_year)),
  individual_id = intdf1$individual_id[1])
preds_g <- predict(mgam, newdata_g, type = "link", se.fit = TRUE, 
                   exclude = "s(individual_id)")
preds_g <- cbind(newdata_g, preds_g)
preds_g <- preds_g %>%
  mutate(lwr_ci = ilink(fit - (2 * se.fit)),
         upr_ci = ilink(fit + (2 * se.fit)),
         fitted = ilink(fit))
ggplot(preds_g, aes(x = day_of_year, y = fitted)) +
  geom_ribbon(aes(ymin = lwr_ci, ymax = upr_ci), alpha = 0.3) +
  geom_point(data = intdf1, 
             aes(x = day_of_year, y = intensity_midpoint, color = individual_id),
             alpha = 0.3) +
  geom_line() +
  theme_bw()

# Still not a great fit to the data, with a steeper peak that tops off at values
# that seem too high (and way more so for CIs)

# GAM for ordered categorical data (mgcv::gam with ocat family) ---------------#

# Try using each unique values as a category (though we get similar results
# if we group 1-50 values together)
intdf1 <- intdf1 %>%
  mutate(y = case_when(
    intensity_midpoint == 0 ~ 1,
    intensity_midpoint == 1 ~ 2,
    intensity_midpoint == 5 ~ 3,
    intensity_midpoint == 50 ~ 4,
    intensity_midpoint == 500 ~ 5,
    intensity_midpoint == 5000 ~ 6,
  ))
R <- length(unique(intdf1$y))

mocat <- gam(y ~ s(day_of_year, bs = "cr", k = 20) + 
               s(individual_id, bs = "re"), 
             data = intdf1, method = "REML", 
             family = ocat(R = R), select = TRUE)
summary(mocat)
gam.check(mocat)
plot(mocat, pages = 1)

# Plot predictions
ilink <- family(mocat)$linkinv # Don't really need this (identity link)
newdata_g <- data.frame(
  day_of_year = seq(min(intdf1$day_of_year), max(intdf1$day_of_year)),
  individual_id = intdf1$individual_id[1])
preds_g <- predict(mocat, newdata = newdata_g, type = "link", se.fit = TRUE,
                   exclude = "s(individual_id)")
preds_g <- cbind(newdata_g, preds_g)
preds_g <- preds_g %>%
  mutate(lwr_ci = ilink(fit - (2 * se.fit)),
         upr_ci = ilink(fit + (2 * se.fit)),
         fitted = ilink(fit))

cutpoints <- mocat$family$getTheta(TRUE)
ggplot(preds_g, aes(x = day_of_year, y = fitted)) +
  geom_ribbon(aes(ymin = lwr_ci, ymax = upr_ci), alpha = 0.3) +
  geom_line() +
  geom_hline(yintercept = cutpoints, linetype = "dashed") +
  theme_bw()
ggplot(preds_g, aes(x = day_of_year, y = fitted)) +
  geom_ribbon(aes(ymin = lwr_ci, ymax = upr_ci), alpha = 0.3) +
  geom_line() +
  geom_hline(yintercept = cutpoints, linetype = "dashed") +
  scale_y_continuous(limits = c(-2, 15)) +
  theme_bw()
# Interpretting this is tricky. Can predict dates that you're likely to
# transition from one abundance category to the next, but that's not very 
# intuitive.
  
# Ordered beta regression (glmmTMB with ordbeta family) -----------------------#

# Ordered beta regression is like a mixture model with cutpoints determining 
# extreme categories (0ish, 1ish) and assuming a continuous response on (0, 1). 
# Unlike zero-one inflated betas, don't need to specify separate models and 
# parameters for 0, (0-1), and 1.

# From help for glmmTMB: If your response variable is defined on the 
# closed interval [a,b], need to first transform it to [0,1] via:
# y_scaled <- (y-a)/(b-a)

intdf1 <- intdf1 %>%
  mutate(y_scaled = intensity_midpoint / max(intensity_midpoint))

mob <- glmmTMB(y_scaled ~ s(day_of_year, k = 20) + (1|individual_id), 
               data = intdf1, 
               REML = TRUE, family = ordbeta)
mob 
# Note that upper cutpoint is 0.125, but that is between two highest 
# categories (500 = 0.1; 5000 = 1)
count(intdf1, intensity_midpoint, y_scaled) %>% format(scientific = FALSE)

# Extract inverse link function
ilink <- family(mob)$linkinv
newdata_g <- data.frame(
  day_of_year = seq(min(intdf1$day_of_year), max(intdf1$day_of_year)),
  individual_id = NA)
preds_g <- predict(mob, newdata_g, type = "link", se.fit = TRUE)
preds_g <- cbind(newdata_g, preds_g)
preds_g <- preds_g %>%
  mutate(lwr_ci = ilink(fit - (2 * se.fit)),
         upr_ci = ilink(fit + (2 * se.fit)),
         fitted = ilink(fit))
ggplot(preds_g, aes(x = day_of_year, y = fitted)) +
  geom_ribbon(aes(ymin = lwr_ci, ymax = upr_ci), alpha = 0.3) +
  geom_point(data = intdf1, 
             aes(x = day_of_year, y = y_scaled, color = individual_id),
             alpha = 0.3) +
  geom_line() +
  theme_bw()
# Predicted values don't approach upper boundary, even when lots of observations
# there (Doesn't seem to just be a 2021 issue. See similar things if I model 
# data for all years)

# Ordered beta regression (ordbetareg package) --------------------------------#

# Done in a Bayesian framework using brms functions and Stan

# Settings that we can tweak:
  # Number of chains: 4 is default
  # Number of cores: 1 is default
  # Number of iterations per chain including warmup (iter): default is 2000
  # Number of warmup iterations (warmup): default is iter/2
  # Thinning rate (thin): default is 1
  # Priors (eg, coef_prior_mean = 0, coef_prior_SD = 5)
  # Knots (list containing specific values to be used for basis construction of smooothing terms)
  # Initial values (init default is "0", but could use "random")
  # Can set a seed
  # Sampling parameters (used to eliminate or decrease number of divergent
    # transitions). Best option is to use control = list(adapt_delta = X) 
    # where X is something between 0.8 and 1. If get warning about tree depth 
    # use control = list(max_treedepth = X) where X is > 10.

  # Insufficient sampling, just to get a sense of what we've got
  mobr <- ordbetareg(intensity_midpoint ~ s(day_of_year) + (1|individual_id),
                     data = intdf1,
                     true_bounds = c(0, 5000),
                     cores = 3, 
                     chains = 3,
                     iter = 1000,
                     control = list(adapt_delta = 0.95))
  mobr
  conditional_effects(mobr, effects = "day_of_year")    # This makes sense to me
  conditional_smooths(mobr, smooths = "s(day_of_year)") # This does not
  
  posterior_check <- pp_check_ordbeta(mobr, ndraws = 100)
  posterior_check$discrete
  posterior_check$continuous
  
  # All years:
  # Here allowing different wiggliness each year. If wanted the same wiggliness, 
  # then use s(doy, fyr, bs = 'fs'). See Pederson et al. 2019
  mobr2 <- ordbetareg(intensity_midpoint ~ fyr + s(day_of_year, by = fyr) +
                        (1|individual_id),
                      data = intdf,
                      true_bounds = c(0, 5000),
                      cores = 3, 
                      chains = 3,
                      iter = 1000,
                      control = list(adapt_delta = 0.99))
  mobr2
  conditional_effects(mobr2, effects = "day_of_year:fyr")
  plot_predictions(mobr2, condition = c("day_of_year", "fyr", "fyr")) 

  # Just including year as a random effect?
  mobr3 <- ordbetareg(intensity_midpoint ~ 
                        s(day_of_year) +
                        (1|fyr) +
                        (1|individual_id),
                      data = intdf,
                      true_bounds = c(0, 5000),
                      cores = 3, 
                      chains = 3,
                      iter = 1000,
                      control = list(adapt_delta = 0.99,
                                     max_treedepth = 15))
  mobr3
  conditional_effects(mobr3, effects = "day_of_year")
  
  # Stuff below is from ordbetareg vignette:
  all_draws <- prepare_predictions(mobr)
  cutzero <- plogis(all_draws$dpars$cutzero)
  cutone <- plogis(all_draws$dpars$cutzero + exp(all_draws$dpars$cutone))
  
  intdf1 %>%
    ggplot(aes(x = intensity_midpoint)) +
    geom_histogram(bins = 100) +
    theme_minimal() +
    geom_vline(xintercept = mean(cutzero)*100,linetype=2) +
    geom_vline(xintercept = mean(cutone)*100,linetype=2)
  
# ordbetareg package seems to give better results than the glmmTMB model with
# ordbeta family. However, the ordbetareg/brms approach takes a lot longer and
# can require some fussing to avoid warnings about diverent transitions and/or
# max treedepth. Even with a "better" fit than glmmTMB, I still don't think 
# that these models fit our overdispersed count data well. 
