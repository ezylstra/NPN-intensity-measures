# Exploring intensity data from sites that have collected it consistently
# ER Zylstra
# 25 Feb 2025

library(dplyr)
library(stringr)
library(lubridate)
library(rnpn)

rm(list = ls())

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

#- Download info on phenophases and intensities -------------------------------#

# Phenophase classes
phc <- npn_pheno_classes() %>% 
  rename(class_id = id,
         class_name = name,
         class_description = description) %>%
  # Restricting to plants
  filter(class_id %in% 1:13) %>%
  data.frame()

ph <- npn_phenophases() %>% 
  select(-color) %>%
  rename(class_id = pheno_class_id,
         category = phenophase_category) %>% 
  filter(class_id %in% 1:13) %>%
  arrange(class_id, phenophase_id) %>%
  data.frame() %>%
  left_join(select(phc, class_id, class_name), by = "class_id")


intensity_file <- "npn-data/intensities.csv"

  # Extract list of intensity categories in NPN database
  ia_orig <- npn_abundance_categories() 
    # Returns tibble with category ids, each with a dataframe that lists the value
    # ids and value names (which is what appears in intensity_value column in a
    # status-intensity dataset)
  
  # Remove any blank entries 
  ia_orig <- ia_orig %>%
    filter(!(is.na(category_name) | category_name == ""))
  
  # Create a 2-dimensional dataframe that contains all the data
  cat_ids <- select(ia_orig, category_id, category_name) %>% data.frame()
  cat_values <- ia_orig$category_values
  cat_values <- mapply(cbind, cat_values, "category_id" = cat_ids$category_id, 
                       SIMPLIFY = FALSE)
  cat_values <- mapply(cbind, cat_values, "category_name" = cat_ids$category_name, 
                       SIMPLIFY = FALSE)
  ia_orig <- bind_rows(cat_values) %>%
    select(category_id, category_name, value_id, value_name, value_description)
  
  # NOTE: categories 13 and 41 (Leaf size, in %) seem to identical #############
  
  # Identify whether values are a number/count, percent, or if they're qualitative. 
  # If not qualitative, extract bounding values to calculate an appoximate
  # midpoint.
  ia <- ia_orig %>%
    select(-value_description) %>%
    mutate(value1 = NA,
           value2 = NA,
           type = case_when(
             str_detect(value_name, "%") ~ "percent",
             str_detect(value_name, "[0-9]") ~ "number",
             .default = "qualitative"
           ))
  val12 <- which(colnames(ia) %in% c("value1", "value2"))
  for (i in 1:nrow(ia)) {
    if (str_detect(ia$value_name[i], " to ")) {
      ia[i, val12] <- str_split_fixed(ia$value_name[i], " to ", 2)
      ia[i, val12] <- as.numeric(str_remove(ia[i, val12], ","))
    } else if (str_detect(ia$value_name[i], "-")) {
      ia[i, val12] <- str_split_fixed(ia$value_name[i], "-", 2)
      ia[i, val12[2]] <- str_remove(ia[i, val12[2]], "%")
    } else if (str_detect(ia$value_name[i], "% or more")) {
      ia[i, val12] <- str_remove(ia$value_name[i], "% or more")
    } else if (str_detect(ia$value_name[i], "Less than ")) {
      ia[i, val12[1]] <- 0
      ia[i, val12[2]] <- str_remove(ia$value_name[i], "Less than ")
      ia[i, val12[2]] <- str_remove(ia[i, val12[2]], "%")
    } else if (str_detect(ia$value_name[i], "More than ")) {
      ia[i, val12] <- str_remove(ia$value_name[i], "More than ")
      ia[i, val12[1]] <- as.numeric(str_remove(ia[i, val12[1]], ",")) + 1
      ia[i, val12[2]] <- as.numeric(str_remove(ia[i, val12[2]], ",")) + 1
    }
  }
  ia <- ia %>%
    mutate_at(c("value1", "value2"), as.numeric)
  
  # Assigning a middle-ish value for each range (For count ranges, keeping it to 
  # nice numbers like 5, 50, 500, and 5000)
  ia <- ia %>%
    mutate(mag = nchar(value1) - 1) %>%
    mutate(value = case_when(
      value1 == value2 ~ round(value1),
      type == "number" & value1 == 0 ~ 1,
      type == "number" & value1 != 0 ~ 
        plyr::round_any(rowMeans(across(value1:value2)), 5 * (10 ^ mag)),
      type == "percent" ~ round(rowMeans(across(value1:value2))),
      .default = NA
    )) %>%
    select(-mag)
  
  # Some of these peak intensity categories are confusing (e.g., cats 31-34)
  # Is there an easy way to match up phenophases with intensity categories?
  


# Explore data for one site ---------------------------------------------------#

# Start with Ellen's site
site_id <- 2
si_file <- paste0("npn-data/si-site", site_id, "-", 
                  min(yrs), "-", max(yrs), ".csv")

si <- read.csv(si_file) %>%
  mutate(yr = year(observation_date))

# Look at when each phenophase was recorded, since I think some were only used
# in first year or two
php_yr <- si %>%
  group_by(phenophase_description) %>%
  summarize(nobs = n(),
            yr_first = min(yr),
            yr_last = max(yr)) %>%
  data.frame()
php_yr %>% arrange(yr_last, desc(yr_first))
# A bunch not used after 2010

# A bunch of conifer related phenophases not used after 2019, but I don't see
# new conifer phenophases popping up after that. Stopped monitoring some conifers?
count(filter(si, grepl("conifers", phenophase_description)), common_name)
si %>%
  group_by(common_name) %>%
  summarize(nobs = n(),
            yr_first = min(yr),
            yr_last = max(yr)) %>%
  data.frame() %>%
  arrange(yr_last, yr_first)
# Looks like just balsam fir, which wasn't monitored after 2019

# Some grasses/sedges phenophases just start being used in 2020
count(filter(si, grepl("grasses/sedges", phenophase_description)), common_name)
# This is just Pennsylvania sedge, that was only monitored in 2020-2024

# Remove observations of phenophases that weren't used after 2010. 
# Keep everything else.
php_keep <- php_yr$phenophase_description[php_yr$yr_last > 2010]
si <- si %>%
  filter(phenophase_description %in% php_keep)

# Look at when each intensity category was recorded, since I think some were 
# only used in a year or two
intensity_yr <- si %>%
  group_by(intensity_category_id) %>%
  summarize(nobs = n(),
            yr_first = min(yr),
            yr_last = max(yr)) %>%
  data.frame()
intensity_yr %>% arrange(yr_last, desc(yr_first))
# A bunch of categories just used in 2011 
# Two categories (40, 42) just used in 2012-2015
# For now, will remove any categories that have not been used since 2015
int_keep <- intensity_yr$intensity_category_id[intensity_yr$yr_last > 2015]
si <- si %>%
  filter(intensity_category_id %in% int_keep)

# Extract phenophase and intensity categories in status-intensity dataset, and
# merge with intensity info
intensity_info <- read.csv(intensity_file)
si_phpi <- si %>%
  filter(intensity_value != -9999) %>%
  group_by(phenophase_id, phenophase_description, intensity_category_id,
           intensity_value) %>%
  summarize(nobs = n(),
            nspp = n_distinct(common_name),
            .groups = "keep") %>%
  left_join(select(intensity_info, intensity_category_id, intensity_name,
                   intensity_value, intensity_type, intensity_mid),
            by = c("intensity_category_id", "intensity_value")) %>%
  data.frame()
si_phpi

# Look into intensity categories (and associated phenophases) that are qualitative
phq <- unique(si_phpi$phenophase_description[si_phpi$intensity_type == "qualitative"])
si %>%
  filter(intensity_value != "-9999") %>%
  filter(phenophase_description %in% phq) %>%
  group_by(common_name, species_functional_type, 
           phenophase_description, intensity_value) %>%
  summarize(nobs = n(),
            yr_first = min(yr),
            yr_last = max(yr),
            .groups = "keep") %>%
  data.frame()

si %>% 
  filter(common_name %in% c("Pennsylvania sedge", "paper birch"), 
         phenophase_description == "Pollen release (flowers)",
         intensity_value != "-9999") %>%
  select(common_name, individual_id, phenophase_description, observation_date,
         yr, phenophase_status, intensity_value) %>%
  arrange(common_name, observation_date)


