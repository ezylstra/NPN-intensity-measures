# Exploring intensity data from sites that have collected it consistently
# ER Zylstra
# 28 Feb 2025

library(dplyr)
library(stringr)
library(lubridate)
library(rnpn)
library(ggplot2)

rm(list = ls())

# Set parameters --------------------------------------------------------------#

# Identify potential site(s) of interest (can call rnpn::npn_stations() if needed)
site_ids <- c(2, 55876, 55878, 55859)
  # 2: Ellen's house
  # 55876: NEON: Bonanza Creek in central AK (quaking aspen)
  # 55878: NEON: Delta Junction in central AK (dwarf birch)
  # 55859: NEON: Santa Rita Exp Range (creosote, velvet mesquite, desert zinnia)

# Erin P had mentioned one other site, but it has very few records so ignoring
# for now
  # 47097: ABQ BioPark Botanic Garden (multiple spp: desert shrubs might have >1 
  #   flowering or fruiting periods per year)

# Identify years of interest (through last calendar year)
yrs <- 2009:(year(Sys.Date()) - 1)

# Name of person requesting NPN data
requestor <- "erinz"

# Download and format NPN data ------------------------------------------------#

for (site_id in site_ids) {
  
  si_csv_name <- paste0("npn-data/si-site", site_id, "-",
                        min(yrs), "-", max(yrs), ".csv")

  if (!file.exists(si_csv_name)) {
  
    # Download status/intensity data (one row for each observation of a plant or 
    # animal species and phenophase)
    si_orig <- npn_download_status_data(
      request_source = requestor,
      years = yrs,
      station_ids = site_id,
      climate_data = FALSE,
      additional_fields = c("site_name", 
                            "observedby_person_id",
                            "species_functional_type")
    )
    
    # Remove animal observations and unnecessary fields
    si_orig <- si_orig %>%
      filter(kingdom == "Plantae") %>%
      select(-c(update_datetime, kingdom, abundance_value))
    
    # Write to file
    write.csv(si_orig, si_csv_name, row.names = FALSE)
    rm(si_orig)
    
  }  
} 

# Load information about phenophases, intensity categories --------------------#

# Phenophase-intensity associations (as of 2024)
ph_int <- read.csv("npn-data/phenophases-intensities-2024.csv")

# Values/midpoints for each intensity category
ivalues <- read.csv("npn-data/intensity-values-2024.csv") %>%
  rename(intensity_value = intensity_value_name)

# Extract list of phenophases that were used in 2024
ph_2024 <- sort(unique(ph_int$phenophase_id))

# Extract list of intensity categories that were used in 2024
int_2024 <- sort(unique(ivalues$intensity_category_id))

# Load status-intensity data for one site -------------------------------------#

# Starting with Ellen's site for now
site_id <- 2
si_file <- paste0("npn-data/si-site", site_id, "-", 
                  min(yrs), "-", max(yrs), ".csv")

si <- read.csv(si_file) %>%
  # Create yr column
  mutate(yr = year(observation_date)) %>%
  # Convert html coding of ">" to symbol
  mutate(phenophase_description = str_replace_all(phenophase_description,
                                                  "&gt;", ">")) %>%
  # Create phenophase_descrip that's the same as phenophase_description except 
  # all parenthetical references to taxonomic groups or locations are removed
  mutate(phenophase_descrip = str_replace(phenophase_description,
                                          " \\s*\\([^\\)]+\\)", "")) %>%
  # Remove any records with unknown phenophase status
  filter(phenophase_status != -1)

# Filter out phenophases or intensity categories that weren't used in 2024 
# (keeping -9999s for now)
si <- si %>%
  filter(phenophase_id %in% ph_2024) %>%
  filter(intensity_category_id %in% c(int_2024, -9999))

count(si, yr, intensity_category_id == -9999)
# Starting in 2013, there was always an intensity category provided.
# Will exclude few early years when intensity category wasn't always specified
si <- si %>%
  filter(yr > 2012)

# Merge phenophase and intensity information with si 
ph_merge <- ph_int %>%
  distinct(phenophase_id, class_id, class_name)
ivalues_missing <- ivalues %>%
  distinct(intensity_category_id, intensity_name, intensity_type) %>%
  mutate(intensity_value_id = NA,
         intensity_value = "-9999",
         intensity_midpoint = NA) %>%
  relocate(intensity_type, .after = "intensity_value")
ivalues <- rbind(ivalues, ivalues_missing)


si <- si %>%
  select(-c(observation_id, observedby_person_id, site_name, latitude,
            longitude, elevation_in_meters, state, genus, species,
            phenophase_descrip)) %>%
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

# Look at data for one species ------------------------------------------------#

spp <- "sugar maple"

# Notes based on activity curves (prop of individuals with yes) for sugar maple 
# in NE US in 2022: 
# Sugar maples flower before leafing out. 
# Fruits appear mid-year but may not ripen until winter? (Not as sure about this)

sispp <- filter(si, common_name == spp)
count(sispp, individual_id) # 3 individuals

sispp %>%
  group_by(class_id, phenophase_description, intensity_label) %>%
  summarize(nobs = n(),
            yr_first = first(yr),
            yr_last = last(yr),
            # Number of observations in phenophase (yeses)
            n_inphase = sum(phenophase_status == 1),
            # Number of observations with an intensity value
            n_intvalue = sum(intensity_value != "-9999"),
            # Proportion of observations
            prop_inphase = round(n_inphase / nobs, 2),
            prop_intvalue = round(n_intvalue / n_inphase, 2),
            .groups = "keep") %>%
  data.frame()

# Calculate the interval between observations of the same plant, phenophase
sispp_obs <- sispp %>%
  distinct(common_name, individual_id, phenophase_description, day_of_year, yr) %>%
  arrange(individual_id, phenophase_description, yr, day_of_year) 
sispp_obs$interval <- NA
for (i in 2:nrow(sispp_obs)) {
  sispp_obs$interval[i] <- ifelse(
    sispp_obs$individual_id[i] == sispp_obs$individual_id[i - 1] &
      sispp_obs$phenophase_description[i] == sispp_obs$phenophase_description[i - 1] &
      sispp_obs$yr[i] == sispp_obs$yr[i - 1], 
    sispp_obs$day_of_year[i] - sispp_obs$day_of_year[i - 1], 
    NA
  )
}

# NEXT UP:
# Summarize amount and quality of information for each plant, phenophase, year
# Number of observations
# Mean interval between observations (or maybe largest gap?)
# Any in phase (status == yes)?
# Number of observations with intensity values when status == 1





# Things to explore: ###########################################################
# Status data for each phenophase: yes/no by year, individual (how much variation?)
# Intensity data for each phenophase: Changes over time for an individual, by year
# What are we gaining with intensity data? By comparing status and intensity value
  # might get a more nuanced understanding fo phenology (eg, how quickly does 
  # a tree go from no leaves to full canopy?)
# Combining intensity values for flowers (#), open flowers (%) = # open flowers
# Combining intensity values for fruits (#), ripe fruits (%) = # ripe fruits

# Simplified dataframe for plotting
sispp1 <- sispp %>%
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
  mutate(intensity = ifelse(status == 0, 0, intensity)) %>%
  # Remove phenophases in a given year if the status is always no
  group_by(phenophase, yr) %>%
  mutate(observed = ifelse(sum(status) > 0, 1, 0)) %>%
  mutate(values = ifelse(sum(intensity) > 0, 1, 0)) %>%
  ungroup() %>%
  data.frame()

# Plot intensity values for each phenophase, colored by plant -----------------#
# Leaves
ggplot(filter(sispp1, yr == 2021 & class %in% leaf_classes & values == 1)) +
  geom_line(aes(x = doy, y = intensity, color = id)) +
  geom_point(aes(x = doy, y = intensity, color = id)) +
  facet_grid(intensity_label ~ ., scales = "free_y", 
             labeller = label_wrap_gen(14)) + 
  labs(title = paste0(str_to_sentence(spp), " - site ", site_id,
                      " - 2021: Leaf phenophases"), 
       x = "Day of year", y = "Estimate", color = "Plant ID") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

# Flowers
ggplot(filter(sispp1, yr == 2021 & class %in% flower_classes & values == 1)) +
  geom_line(aes(x = doy, y = intensity, color = id)) +
  geom_point(aes(x = doy, y = intensity, color = id)) +
  facet_grid(intensity_label ~ ., scales = "free_y", 
             labeller = label_wrap_gen(14)) + 
  labs(title = paste0(str_to_sentence(spp), " - site ", site_id,
                      " - 2021: Flower phenophases"), 
       x = "Day of year", y = "Estimate", color = "Plant ID") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

# Fruits
ggplot(filter(sispp1, yr == 2021 & class %in% fruit_classes & values == 1)) +
  geom_line(aes(x = doy, y = intensity, color = id)) +
  geom_point(aes(x = doy, y = intensity, color = id)) +
  facet_grid(intensity_label ~ ., scales = "free_y", 
             labeller = label_wrap_gen(14)) + 
  labs(title = paste0(str_to_sentence(spp), " - site ", site_id, 
                      " - 2021: Fruit phenophases"), 
       x = "Day of year", y = "Estimate", color = "Plant ID") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

# How much annual variation is there in intensity values, over years
# (within a plant, aggregating across plants)
ggplot(filter(sispp1, class == 3 & values == 1)) +
  geom_line(aes(x = doy, y = intensity, color = factor(yr))) +
  geom_point(aes(x = doy, y = intensity, color = factor(yr))) +
  facet_grid(id ~ .) + 
  labs(title = paste0(str_to_sentence(spp), " - site ", site_id, 
                      ": Leaf canopy fullness (%)"), 
       x = "Day of year", y = "Estimate", color = "Year") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")
# Need to remove a year if there's large gaps between observations #######

# Looks like there's more variation among individuals than years, which is
# interesting.

