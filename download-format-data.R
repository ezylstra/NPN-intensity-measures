# Exploring intensity data from sites that have collected it consistently
# ER Zylstra
# 27 Feb 2025

library(dplyr)
library(stringr)
library(lubridate)
library(rnpn)

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
  left_join(select(ph_merge, phenophase_id, class_id, class_name), 
            by = "phenophase_id") %>%
  left_join(ivalues, by = c("intensity_category_id", "intensity_value"))

# Look at data for one species ------------------------------------------------#

spp <- "sugar maple"

sispp <- filter(si, common_name == spp)
sispp %>%
  group_by(class_id, phenophase_description, intensity_name, intensity_type) %>%
  summarize(nobs = n(),
            yr_first = first(yr),
            yr_last = last(yr),
            n_inphase = sum(phenophase_status == 1),
            n_intvalue = sum(intensity_value != "-9999"),
            prop_inphase = round(n_inphase / nobs, 2),
            prop_intvalue = round(n_intvalue / n_inphase, 2),
            .groups = "keep") %>%
  data.frame()



