# Exploring intensity data for canopy trees in southern Appalachian Trail area
# ER Zylstra

library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(brms)
library(leaflet)
# library(geosphere)
# library(mgcv)
# library(cowplot)


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
    max_int > 30 ~ 1,
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
  group_by(site_id, latitude, longitude) %>%
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



  
