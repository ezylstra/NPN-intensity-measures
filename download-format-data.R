# Downloading and formatting data for sites that have collected intensity data 
# consistently

# ER Zylstra
# 14 July 2025

library(dplyr)
library(stringr)
library(lubridate)
library(rnpn)

# Identify sites --------------------------------------------------------------#

# Ellen's house has site_id = 2 (which shouldn't change)
  # bbox <- c(-70.70, 43.05, -70.69, 43.10)
  # wkt <- spocc::bbox2wkt(bbox)
  # npn_stations_by_location(wkt) 
  
# Kodiak NWR (network_id == 1257)
  # Collected data at 12 sites to monitor berry (fruit) resources
  # I think the data were collected via remote cameras that were checked every 3 or so days
  # Species: Site IDs
  # devilsclub: 53455, 53456
  # oval-leaf blueberry: 53453, 53454
  # red elderberry: 53457, 53458, 53459, 53460
  # salmonberry: 53461, 53462, 53463, 53464
  kodiak_sites <- npn_stations() %>%
    filter(network_id == 1257) %>%
    data.frame() %>%
    rename(site_name = station_name,
           site_id = station_id) %>%
    mutate(site_short = paste0("KNWR", str_pad(1:12, 2, pad = "0"))) %>%
    mutate(species = case_when(
      str_detect(site_name, "Blueberry") ~ "oval-leaf blueberry",
      str_detect(site_name, "Devilsclub") ~ "devilsclub",
      str_detect(site_name, "Elderberry") ~ "red elderberry",
      str_detect(site_name, "Salmonberry") ~ "salmonberry"
    )) %>%
    select(site_name, site_id, site_short, species)

# NEON sites
  # All associated with network_id == 77, network_name = National Ecological Observatory Network (NEON)
  # Note that unlike other sites, NEON site IDs change annually when new data is imported.

  # Names of sites we want: 
  # BONA_092.phenology.phe - primary: Bonanza Creek in central AK (quaking aspen). Current site ID = 57121
  # DEJU_065.phenology.phe - primary: Delta Junction in central AK (dwarf birch). Current site ID = 57123
  # SRER_060.phenology.phe - phenocam: Santa Rita Exp Range (creosote, velvet mesquite, desert zinnia) Current site ID = 57104
  
  # Not sure if we also want: 
  # BONA_092.phenology.phe - phenocam: Current site ID = 57122
  # DEJU_065.phenology.phe - phenocam: Current site ID = 57124
  # SRER_060.phenology.phe - primary: Current site ID = 57103

  neons <- data.frame(site_short = c("BONA", "DEJU", "SRER"),
                           species = c("quaking aspen", 
                                       "dwarf birch",
                                       "creosote, velvet mesquite, desert zinnia"))
  
  neon_sites <- npn_stations() %>%
    filter(network_id == 77) %>%
    data.frame() %>%
    filter(str_detect(station_name, paste(neons$site_short, collapse = "|"))) %>%
    select(station_name, station_id) %>%
    arrange(station_name) %>%
    rename(site_name = station_name,
           site_id = station_id) %>%
    mutate(site_short = str_sub(site_name, 1, 4)) %>%
    left_join(neons, by = "site_short")
  
# McDowell Sonoran Conservancy (network_id == 622)
  # Collected data at 5 sites (though two are "Training", "Beta Testing" that
  # we'll ignore). Will download data from other sites separately since file 
  # sizes are large

  mcdo_sites <- npn_stations() %>%
    filter(network_id == 622) %>%
    data.frame() %>%
    filter(station_name != "Training") %>%
    filter(station_name != "Beta Testing") %>%
    select(station_name, station_id) %>%
    arrange(station_id) %>%
    rename(site_name = station_name,
           site_id = station_id) %>%
    mutate(site_short = paste0("MCDO", 1:3), 
           species = "various")

# Combine all
  all_sites <- neon_sites %>%
    rbind(kodiak_sites) %>%
    rbind(mcdo_sites) %>%
    rbind(data.frame(site_name = "Home",
                     site_id = 2,
                     site_short = "ELLEN",
                     species = "various"))

# Download NPN status-intensity data and save to file -------------------------#

# Identify years of interest (through last calendar year)
yrs <- 2009:(year(Sys.Date()) - 1)
  
# Name of person requesting NPN data
requestor <- "erinz"

# For now, need to download each year separately and then combine because if
# we request the entire yrs range in the function call, the download will abort
# with errors when it encounters a year with no data.

for (site in unique(all_sites$site_short)) {
  
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
          station_ids = c(all_sites$site_id[all_sites$site_short == site]),
          climate_data = FALSE,
          additional_fields = c("site_name", 
                                "observedby_person_id",
                                "species_functional_type"))
        )
      }, 
      error = function(e) {
        cat("Note: No data to download for", site, "in", yr, "\n")
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
# (created in phenophases-intensities.R)

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

# Create phenophase and intensity dataframes that can be merged with 
# status-intensity data
ph_merge <- ph_int %>%
  distinct(phenophase_id, class_id, class_name)
ivalues_missing <- ivalues %>%
  distinct(intensity_category_id, intensity_name, intensity_type) %>%
  mutate(intensity_value = NA,
         intensity_midpoint = NA) %>%
  relocate(intensity_type, .after = "intensity_value")
ivalues <- rbind(ivalues, ivalues_missing)

# Load status-intensity data for each site and format -------------------------#

for (site in unique(all_sites$site_short)) {

  si_file <- paste0("npn-data/si-site", site, "-", 
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
  
  # For all NEON and KWR sites, extract data for species of interest. 
  # (Not sure what species will be most useful at Ellen's site and will 
  # keep data for all species at McDowell for now)
  if (!site %in% c("ELLEN", paste0("MCDO", 1:3))) {
    si <- si %>%
      filter(str_detect(common_name, 
                        str_replace_all(all_sites$species[all_sites$site_short == site][1], 
                                        ", ", "|")))
  }

  # Filter out phenophases or intensity categories that weren't used in 2024 
  # or phenophases that don't have intensity categories associated with them
  # (for now, keeping phenophases that could have intensity categories, but 
  # observers didn't report intensity values)
  si <- si %>%
    filter(phenophase_id %in% ph_2024) %>%
    filter(intensity_category_id %in% int_2024 | is.na(intensity_category_id))

  # Merge phenophase and intensity information with si 
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
  
  min_yr <- min(si$yr)
  max_yr <- max(si$yr)
  
  # Save filtered and formatted data to file
  new_filename <- paste0("npn-data/intensity-site", site, "-", 
                         min_yr, "-", max_yr, ".csv")
  write.csv(si, new_filename, row.names = FALSE)
 
}

# Download desert plants data to complement McDowell --------------------------#

# Species: saguaro, buck-horn cholla, jojoba, velvet mesquite, California 
# barrel cactus
dspp <- npn_species() %>% data.frame()
dspp_id <- dspp %>%
  filter(common_name %in% c("saguaro", "buck-horn cholla", "jojoba",
                            "velvet mesquite", "California barrel cactus")) %>%
  pull(species_id)

for (i in 1:length(dspp_id)) {
  
  si_csv_name <- paste0("npn-data/si-spp", dspp_id[i], "-",
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
                 species_ids = dspp_id[i],
                 climate_data = FALSE,
                 additional_fields = c("site_name", 
                                       "observedby_person_id",
                                       "species_functional_type"))
        )
      }, 
      error = function(e) {
        cat("Note: No data to download for", 
            dspp$common_name[dspp$species_id == dspp_id[i]], "in", yr, "\n")
      })
    }
    
    # Combine data for all years
    si_files <- ls()[ls() %in% paste0("si_", yrs)]
    si_orig <- do.call(rbind, mget(si_files))
    
    # Remove unnecessary fields
    si_orig <- si_orig %>%
      select(-c(update_datetime, kingdom, abundance_value))
    
    # Write to file
    write.csv(si_orig, si_csv_name, row.names = FALSE)
    rm(list = c("si_orig", si_files))
    
  }  
}

# Load status-intensity data for each species and format ----------------------#

for (i in 1:length(dspp_id)) {
  
  si_file <- paste0("npn-data/si-spp", dspp_id[i],  "-", 
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
  # or phenophases that don't have intensity categories associated with them
  # (for now, keeping phenophases that could have intensity categories, but 
  # observers didn't report intensity values)
  si <- si %>%
    filter(phenophase_id %in% ph_2024) %>%
    filter(intensity_category_id %in% int_2024 | is.na(intensity_category_id))
  
  # Merge phenophase and intensity information with si 
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
  
  min_yr <- min(si$yr)
  max_yr <- max(si$yr)
  
  # Save filtered and formatted data to file
  new_filename <- paste0("npn-data/intensity-spp", dspp_id[i], "-", 
                         min_yr, "-", max_yr, ".csv")
  write.csv(si, new_filename, row.names = FALSE)
  
}
