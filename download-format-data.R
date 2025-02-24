# Exploring intensity data from sites that have collected it consistently
# ER Zylstra
# 24 Feb 2025

library(dplyr)
library(stringr)
library(lubridate)
library(rnpn)


# Set parameters --------------------------------------------------------------#

# Identify site(s) of interest (start with Ellen's site for now)
site_ids <- 2

# Identify years of interest (through last calendar year)
yrs <- 2009:(year(Sys.Date()) - 1)

# Name of person requesting NPN data
requestor <- "erinz"

# Use NPN-LPP-reports/phenophases-intensities.R to clean up or aggregate classes?

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
      station_ids = site_ids,
      climate_data = FALSE,
      additional_fields = c("site_name", "observedby_person_id")
    )
    
    # Remove animal observations and unnecessary fields
    si_orig <- si_orig %>%
      filter(kingdom == "Plantae") %>%
      select(-c(update_datetime, kingdom))
    
    # Write to file
    write.csv(si_orig, si_csv_name, row.names = FALSE)
    rm(si_orig)
    
  }  
} 

# Explore data for a plant or two ---------------------------------------------#

# Start with oak trees at Ellen's site (one black oak, one northern red oak)
site_id <- 2

si_file <- paste0("npn-data/si-site", site_id, "-", 
                  min(yrs), "-", max(yrs), ".csv")

si <- read.csv(si_file) %>%
  filter(genus == "Quercus") %>%
  mutate(yr = year(observation_date))

# Look at when each phenophase was recorded, since I think some were only used
# in first year or two
php_yr <- si %>%
  group_by(phenophase_description) %>%
  summarize(nobs = n(),
            yr_first = min(yr),
            yr_last = max(yr))
php_yr

# Remove observations of phenophases that weren't used after 2010
php_keep <- php_yr$phenophase_description[php_yr$yr_last > 2010]
si <- si %>%
  filter(phenophase_description %in% php_keep)

