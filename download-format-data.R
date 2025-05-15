# Exploring intensity data from sites that have collected it consistently
# ER Zylstra
# 14 May 2025

library(dplyr)
library(stringr)
library(lubridate)
library(rnpn)
library(ggplot2)

rm(list = ls())

# Set parameters --------------------------------------------------------------#

# Identify potential site(s) of interest (can call rnpn::npn_stations() if needed)

# Ellen's house:
  # bbox <- c(-70.70, 43.05, -70.69, 43.10)
  # wkt <- spocc::bbox2wkt(bbox)
  # npn_stations_by_location(wkt) 
  # Site 2 (this shouldn't change)
  
# NEON sites in AK and AZ:
  # All associated with network_id == 77; network_name = National Ecological Observatory Network (NEON)
  # Note that unlike other sites, NEON site IDs change annually when new data is imported.
  neons <- npn_stations(state_code = c("AK", "AZ")) %>%
    filter(network_id == 77) %>%
    data.frame()
  # neons
  
  # Names of sites we want: 
  # BONA_092.phenology.phe - primary: Bonanza Creek in central AK (quaking aspen). Site ID = 57121
  # DEJU_065.phenology.phe - primary: Delta Junction in central AK (dwarf birch). Site ID = 57123
  # SRER_060.phenology.phe - phenocam: Santa Rita Exp Range (creosote, velvet mesquite, desert zinnia) Site ID = 57104
  
  # Not sure if we also want: 
  # BONA_092.phenology.phe - phenocam: Site ID = 57122
  # DEJU_065.phenology.phe - phenocam: Site ID = 57124
  # SRER_060.phenology.phe - primary: Site ID = 57103

site_table <- data.frame(
  site_short = c("Ellen", rep("BONA", 2), rep("DEJU", 2), rep("SRER", 2)),
  site_name = c("Home", 
                "BONA_092.phenology.phe - primary",
                "BONA_092.phenology.phe - phenocam",
                "DEJU_065.phenology.phe - primary",
                "DEJU_065.phenology.phe - phenocam",
                "SRER_060.phenology.phe - primary",
                "SRER_060.phenology.phe - phenocam"),
  site_id = c(2, 57121, 57122, 57123, 57124, 57103, 57014))

# Identify years of interest (through last calendar year)
yrs <- 2009:(year(Sys.Date()) - 1)

# Name of person requesting NPN data
requestor <- "erinz"

# Download and format NPN data ------------------------------------------------#

# For now, need to download each year separately and then combine, because if
# we request the entire yrs range in the function call, the download will abort
# with errors when it encounters a year with no data.

for (site in site_table$site_short) {
  
  si_csv_name <- paste0("npn-data/si-site", site, "-",
                        min(yrs), "-", max(yrs), ".csv")

  if (!file.exists(si_csv_name)) {
    
    for (yr in yrs) {
  
      # Download status/intensity data (one row for each observation of a plant
      # or animal species and phenophase)
      
      tryCatch({
        assign(paste0("si_", yr), 
          npn_download_status_data(
          request_source = requestor,
          years = yr,
          station_ids = c(site_table$site_id[site_table$site_short == site]),
          climate_data = FALSE,
          additional_fields = c("site_name", 
                                "observedby_person_id",
                                "species_functional_type"))
        )
      }, 
      error = function(e) {
        cat("Note: No data to download for ", site, " in ", yr, "\n")
      })
    }
    
    # Combine data for all years
    si_files <- ls()[ls() %in% paste0("si_", yrs)]
    si_orig <- do.call(rbind, mget(si_files))
    
    # Remove animal observations and unnecessary fields
    si_orig <- si_orig %>%
      filter(kingdom == "Plantae") %>%
      select(-c(update_datetime, kingdom, abundance_value))
    
    # Write to file
    write.csv(si_orig, si_csv_name, row.names = FALSE)
    rm(list = c("si_orig", si_files))
    
  }  
} 

# Load information about phenophases, intensity categories --------------------#

# Phenophase-intensity associations (as of 2024)
ph_int <- read.csv("npn-data/phenophases-intensities-2024.csv")

# Values/midpoints for each intensity category
ivalues <- read.csv("npn-data/intensity-values-2024.csv") %>%
  rename(intensity_value = intensity_value_name) %>%
  select(-intensity_value_id)

# Extract list of phenophases that were used in 2024
ph_2024 <- sort(unique(ph_int$phenophase_id))

# Extract list of intensity categories that were used in 2024
int_2024 <- sort(unique(ivalues$intensity_category_id))

# Load status-intensity data for one site -------------------------------------#

# Starting with Ellen's site for now
si_file <- paste0("npn-data/si-siteEllen-", 
                  min(yrs), "-", max(yrs), ".csv")

si <- read.csv(si_file) %>%
  # Create yr column
  mutate(yr = year(observation_date)) %>%
  # Convert html coding of ">" to symbol
  mutate(phenophase_description = str_replace_all(phenophase_description,
                                                  "&gt;", ">")) %>%
  # Remove whitespaces in common_name column
  mutate(common_name = str_trim(common_name)) %>%
  # Create phenophase_descrip that's the same as phenophase_description except 
  # all parenthetical references to taxonomic groups or locations are removed
  # mutate(phenophase_descrip = str_replace(phenophase_description,
  #                                         " \\s*\\([^\\)]+\\)", "")) %>%
  # Remove any records with unknown phenophase status
  filter(phenophase_status != -1) 

# Filter out phenophases or intensity categories that weren't used in 2024 
# (keeping NAs for now)
si <- si %>%
  filter(phenophase_id %in% ph_2024) %>%
  filter(intensity_category_id %in% int_2024 | is.na(intensity_category_id))

# count(si, yr, is.na(intensity_category_id))
# Starting in 2013, there was always an intensity category provided.
# Will exclude few early years when intensity category wasn't always specified
si <- si %>%
  filter(yr > 2012)

# Merge phenophase and intensity information with si 
ph_merge <- ph_int %>%
  distinct(phenophase_id, class_id, class_name)
ivalues_missing <- ivalues %>%
  distinct(intensity_category_id, intensity_name, intensity_type) %>%
  mutate(intensity_value = NA,
         intensity_midpoint = NA) %>%
  relocate(intensity_type, .after = "intensity_value")
ivalues <- rbind(ivalues, ivalues_missing)

si <- si %>%
  select(-c(observation_id, observedby_person_id)) %>%
  left_join(ph_merge, by = "phenophase_id") %>%
  left_join(ivalues, by = c("intensity_category_id", "intensity_value"))

# Create a new column with intensity labels (factor)
si <- si %>%
  # Remove anything in parentheses in intensity name
  mutate(intensity_label = str_replace(intensity_name, " \\s*\\([^\\)]+\\)", "")) %>%
  # Remove the word " present" from intensity name
  mutate(intensity_label = str_remove(intensity_label, " present")) %>%
  # Remove the word "Potential" from intensity name
  mutate(intensity_label = str_remove(intensity_label, "Potential ")) %>%
  # Remove the word " percentage" from intensity name
  mutate(intensity_label = str_remove(intensity_label, " percentage")) %>%
  # Remove "Recent " from intensity name
  mutate(intensity_label = str_remove(intensity_label, "Recent ")) %>%
  # Replace " or " with "/" in intensity name
  mutate(intensity_label = str_replace(intensity_label, " or ", "/")) %>%
  # Add "No. " in front or "(%)" at the end
  mutate(intensity_label = case_when(
    intensity_type == "number" ~ paste0("No. ", str_to_lower(intensity_label)),
    intensity_type == "percent" ~ paste0(str_to_sentence(intensity_label), " (%)"),
    .default = intensity_label
  )) %>%
  arrange(class_id) %>%
  mutate(intensity_label = factor(intensity_label, 
                                  levels = unique(intensity_label)))

leaf_classes <- 1:4
flower_classes <- 6:8
fruit_classes <- 10:13

# Create plant list for Ellen's site (after filtering by monitoring freq, intensity data)
site_plants <- si %>%
  group_by(site_id, site_name, latitude, longitude,
           elevation_in_meters, state, species_id, genus, species, common_name,
           species_functional_type, individual_id) %>%
  summarize(n_obs = n(),
            first_year = min(yr),
            last_year = max(yr),
            .groups = "keep") %>%
  arrange(species_functional_type, common_name, .locale = "en") %>%
  data.frame()
# Remove site and species information we don't need to simplify dataframes
si <- si %>%
  select(-c(site_name, latitude, longitude, elevation_in_meters, state,
            species_id, genus, species, species_functional_type))

# Look at data for one species ------------------------------------------------#

spp <- "sugar maple"

# Notes based on activity curves (prop of individuals with yes) for sugar maple 
# in NE US in 2022: 
  # Sugar maples flower before leafing out. 
  # Fruits appear mid-year but may not ripen until winter?

sispp <- filter(si, common_name == spp)
count(sispp, individual_id) # 3 individuals

# Look at amount of observations in phase, with intensity values by phenophase
sispp %>%
  group_by(class_id, phenophase_description, intensity_label) %>%
  summarize(nobs = n(),
            yr_first = first(yr),
            yr_last = last(yr),
            # Number of observations in phenophase
            n_inphase = sum(phenophase_status == 1),
            # Number of observations with an intensity value
            n_intvalue = sum(!is.na(intensity_value)),
            # Proportion of observations in phenophase
            prop_inphase = round(n_inphase / nobs, 2),
            # Proportion of in-phase observations with intensity values
            prop_intvalue = round(n_intvalue / n_inphase, 2),
            .groups = "keep") %>%
  data.frame()
# No observations of Pollen release

# First check that there's only one observation of each plant-phenophase per day
  # nrow(sispp)
  # nrow(distinct(sispp, individual_id, phenophase_description, day_of_year, yr))
# If these have the same number of rows, we're good

# Calculate the interval between observations of the same plant, phenophase
sispp <- sispp %>%
  arrange(individual_id, phenophase_description, yr, day_of_year)
sispp$interval <- NA
for (i in 2:nrow(sispp)) {
  sispp$interval[i] <- ifelse(
    sispp$individual_id[i] == sispp$individual_id[i - 1] &
      sispp$phenophase_description[i] == sispp$phenophase_description[i - 1] &
      sispp$yr[i] == sispp$yr[i - 1], 
    sispp$day_of_year[i] - sispp$day_of_year[i - 1], 
    NA
  )
}

# Summarize amount and quality of information for each plant, phenophase, year
pl_ph_yr <- sispp %>%
  group_by(common_name, individual_id, phenophase_description, yr) %>%
  summarize(nobs = n(),
            first_obs = min(day_of_year),
            last_obs = max(day_of_year),
            mean_int = round(mean(interval, na.rm = TRUE), 2),
            max_int = max(interval, na.rm = TRUE),
            n_inphase = sum(phenophase_status),
            n_intvalue = sum(!is.na(intensity_value)),
            prop_inphase = round(n_inphase / nobs, 2),
            prop_intvalue = round(n_intvalue / n_inphase, 2),
            .groups = "keep") %>%
  data.frame()

# Yearly differences: Less frequent monitoring in 2017, all other years similar
pl_ph_yr %>%
  group_by(yr) %>%
  summarize(mean_int = round(mean(mean_int), 1))
# Less frequent monitoring in 2017.
# Slightly more frequent monitoring in 2018-2024 than 2013-2016

# Phenophase differences:
pl_ph_yr %>%
  group_by(phenophase_description) %>%
  summarize(mean_int = round(mean(mean_int), 1))
# Not much difference among monitoring frequency among phenophases

# So for exploring intensity data for sugar maple at Ellen's site, we will:
  # remove plant-phenophase-year combos with no observations in phase (includes all pollen release)
  # remove plant-phenophase-year combos with no observations with intensity values
  # remove plant-phenophase-year combos when mean interval > 14 days (2017 only)

pl_ph_yr <- pl_ph_yr %>%
  mutate(remove = case_when(
    n_inphase == 0 ~ 1,
    n_intvalue == 0 ~ 1,
    mean_int > 14 ~ 1,
    .default = 0
  ))
sispp2 <- sispp %>%
  left_join(select(pl_ph_yr, individual_id, phenophase_description, yr, remove),
            by = c("individual_id", "phenophase_description", "yr")) %>%
  filter(remove == 0) %>%
  select(-remove)

# Things to explore: ###########################################################
# Status data for each phenophase: yes/no by year, individual (how much variation?)
# Intensity data for each phenophase: Changes over time for an individual, by year
# What are we gaining with intensity data? By comparing status and intensity value
  # might get a more nuanced understanding of phenology (eg, how quickly does 
  # a tree go from no leaves to full canopy?)
# Combining intensity values for flowers (#), open flowers (%) = # open flowers
# Combining intensity values for fruits (#), ripe fruits (%) = # ripe fruits

# Simplified dataframe for plotting
sispp_plot <- sispp2 %>%
  select(common_name, individual_id, class_id, phenophase_description, 
         yr, observation_date, day_of_year, phenophase_status, intensity_label,
         intensity_midpoint) %>%
  rename(id = individual_id, 
         class = class_id, 
         phenophase = phenophase_description,
         obsdate = observation_date,
         doy = day_of_year,
         status = phenophase_status,
         intensity = intensity_midpoint) %>%
  mutate(obsdate = ymd(obsdate), 
         id = as.factor(id)) %>%
  # Change intensity values to 0 (from NA), when status = 0
  mutate(intensity = ifelse(status == 0, 0, intensity))

# Plot intensity values for each phenophase, colored by plant -----------------#
# Leaves
ggplot(filter(sispp_plot, yr == 2021 & class %in% leaf_classes)) +
  geom_line(aes(x = doy, y = intensity, color = id)) +
  geom_point(aes(x = doy, y = intensity, color = id)) +
  facet_grid(intensity_label ~ ., scales = "free_y", 
             labeller = label_wrap_gen(14)) + 
  labs(title = paste0(str_to_sentence(spp), " - Ellen's site,",
                      " - 2021: Leaf phenophases"), 
       x = "Day of year", y = "Estimate", color = "Plant ID") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

# Flowers
ggplot(filter(sispp_plot, yr == 2021 & class %in% flower_classes)) +
  geom_line(aes(x = doy, y = intensity, color = id)) +
  geom_point(aes(x = doy, y = intensity, color = id)) +
  facet_grid(intensity_label ~ ., scales = "free_y", 
             labeller = label_wrap_gen(14)) + 
  labs(title = paste0(str_to_sentence(spp), " - Ellen's site,",
                      " - 2021: Flower phenophases"), 
       x = "Day of year", y = "Estimate", color = "Plant ID") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

# Fruits
ggplot(filter(sispp_plot, yr == 2021 & class %in% fruit_classes)) +
  geom_line(aes(x = doy, y = intensity, color = id)) +
  geom_point(aes(x = doy, y = intensity, color = id)) +
  facet_grid(intensity_label ~ ., scales = "free_y", 
             labeller = label_wrap_gen(14)) + 
  labs(title = paste0(str_to_sentence(spp), " - Ellen's site,",
                      " - 2021: Fruit phenophases"), 
       x = "Day of year", y = "Estimate", color = "Plant ID") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

# How much annual variation is there in intensity values, over years
# (within a plant, aggregating across plants)
ggplot(filter(sispp_plot, class == 3)) +
  geom_line(aes(x = doy, y = intensity, color = factor(yr))) +
  geom_point(aes(x = doy, y = intensity, color = factor(yr))) +
  facet_grid(id ~ .) + 
  labs(title = paste0(str_to_sentence(spp), " - Ellen's site,",
                      ": Leaf canopy fullness (%)"), 
       x = "Day of year", y = "Estimate", color = "Year") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

# Looks like there's more variation among individuals than years, which is
# interesting.

