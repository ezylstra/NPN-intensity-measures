# Exploring relationships between plant phenophases and intensity categories
# ER Zylstra
# 12 May 2025

library(dplyr)
library(stringr)
library(lubridate)
library(rnpn)

rm(list = ls())

# Create custom rounding function ---------------------------------------------#

# This is the same as plyr::round_any(), but loading the plyr package often 
# causes conflicts with function names when dplyr is loaded

round_any <- function(x, accuracy) {
  round(x / accuracy) * accuracy
}

# Download info on phenophases and phenoclasses -------------------------------#

phc <- npn_pheno_classes() %>% 
  rename(class_id = id,
         class_name = name,
         class_description = description) %>%
  # Restricting to plants
  filter(class_id %in% 1:13) %>%
  # Removing unnecessary column
  select(-sequence) %>%
  data.frame()

ph <- npn_phenophases() %>% 
  select(-color) %>%
  # Remove trailing whitespaces in the phenophase_name column
  mutate(phenophase_name = str_trim(phenophase_name)) %>%
  # Rename columns
  rename(class_id = pheno_class_id,
         category = phenophase_category) %>% 
  # Restricting to plants
  filter(class_id %in% 1:13) %>%
  # Merge with phenophase class information
  left_join(select(phc, class_id, class_name), by = "class_id") %>%
  arrange(class_id, phenophase_id) %>%
  data.frame()

# Download info on intensity categories ---------------------------------------#

ia_orig <- npn_abundance_categories() %>%
  select(category_id, category_name, value_id, value_name) %>%
  data.frame()

# Create new dataframe, that identifies whether values are a number/count, 
# percent, or if they're qualitative. If not qualitative, extract bounding 
# values to calculate an appoximate midpoint.
ia <- ia_orig %>%
  mutate(value1 = NA,
         value2 = NA,
         intensity_type = case_when(
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

# Assigning a middle-ish value for each range (keeping it to nice numbers like 
# 5, 50, 500, and 5000)
ia <- ia %>%
  mutate(mag = nchar(value1) - 1) %>%
  mutate(value = case_when(
    value1 == value2 ~ round(value1),
    intensity_type == "number" & value1 == 0 ~ 1,
    intensity_type == "number" & value1 != 0 ~ 
      round_any(rowMeans(across(value1:value2)), 5 * (10 ^ mag)),
    intensity_type == "percent" ~ round(rowMeans(across(value1:value2))),
    .default = NA
  )) %>%
  select(-c(mag, value1, value2)) %>%
  # Rename columns
  rename(intensity_category_id = category_id,
         intensity_name = category_name,
         intensity_value_id = value_id,
         intensity_value_name = value_name,
         intensity_midpoint = value)

#########################
# Few notes/questions on intensity categories

# Categories 13 and 41 (Leaf size, in %) seem to be identical
filter(ia, intensity_category_id %in% c(13, 41))
  # I think category 13 was replaced by 41

# Intensity categories labled as "peak" (categories 31-34) are confusing 
# because they can have numeric or qualitative responses.
filter(ia, grepl("peak", intensity_name))
filter(ia_orig, grepl("peak", category_name))
  # Doesn't look like these were used in recent years 
#########################

# IMPORTANT NOTE: ##############################################################
# Want to match up phenophases and intensity categories. Originally I did this
# using status intensity data (downloaded all plant data in 2024!). However, I
# later learned about the Secondary Phenophase Details table available in the 
# API that lists all species-phenophase-abundance/intensity categories (they
# call this the SSPI data). Import the SSPI data and compare at the bottom of 
# this script. 
################################################################################

# Download status-intensity data ----------------------------------------------#

# Commenting out much of the code below since I don't want to accidentally start 
# a new download or overwrite the existing file (there are a lot of records and 
# it takes a while!)

# Limit download to plants by finding relevant taxonomic orders first
  all_spp <- npn_species()
  plant_spp <- all_spp %>%
    filter(kingdom == "Plantae") %>%
    select(species_id, common_name, functional_type, class_id, class_name)
  plant_classes <- unique(plant_spp$class_id)
  
# Set year of interest
  yr <- 2024

# # Download status-intensity data (this takes a long time! millions of records)
#   si_orig <- npn_download_status_data(
#     request_source = "erinz",
#     years = yr,
#     class_ids = plant_classes,
#     climate_data = FALSE,
#     additional_fields = c("species_functional_type")
#   )
# 
# # Remove unnecessary columns 
#   si_orig <- si_orig %>%
#     select(-c(update_datetime, elevation_in_meters, genus, species,
#               kingdom, day_of_year, abundance_value))
# 
# # Don't need all the observations. Just want a table that provides combinations
# # of species, phenophase, and intensity category (when specified). 
#   si_combos <- si_orig %>%
#     filter(intensity_category_id != -9999) %>%
#     distinct(species_id, common_name, species_functional_type, phenophase_id,
#              phenophase_description, intensity_category_id, intensity_value)
# # Write this to file
#   write.csv(si_combos,
#             paste0("npn-data/si-plants-", yr, "-spp-ph-int-combos.csv"),
#             row.names = FALSE)

# Load si data and combine with phenophase, intensity information -------------#

# Read in the species-phenophase-intensity category combination file
si_combos <- read.csv(paste0("npn-data/si-plants-", yr, 
                             "-spp-ph-int-combos.csv"))

# Check that each intensity category has multiple values (other than -9999)
si_combos %>%
  filter(intensity_value != "-9999") %>%
  group_by(intensity_category_id) %>%
  summarize(nlevels = n_distinct(intensity_value)) %>%
  data.frame()

# Create better intensity "type" classification to account for "peak" categories
ia <- ia %>%
  group_by(intensity_category_id, intensity_name) %>%
  mutate(type1 = first(intensity_type),
         type2 = last(intensity_type)) %>%
  ungroup() %>%
  mutate(intensity_type = ifelse(type1 == type2, intensity_type, 
                                 paste0(type1, "/", type2))) %>%
  select(-c(type1, type2)) %>%
  data.frame()

# Summarize information about each intensity category
ia_cats <- ia %>%
  group_by(intensity_category_id, intensity_name, intensity_type) %>%
  summarize(nlevels = n(), .groups = "keep") %>%
  data.frame()

# Create dataframe with unique combinations of phenophase and intensity category
# that were used in 2024
ph_int <- si_combos %>%
  group_by(phenophase_id, phenophase_description, intensity_category_id) %>%
  summarize(n_spp = n_distinct(common_name), .groups = "keep") %>%
  left_join(select(ph, phenophase_id, category, class_id, class_name),
            by = "phenophase_id") %>%
  left_join(ia_cats, by = "intensity_category_id") %>%
  data.frame() %>%
  arrange(class_id, phenophase_id)

# Write this to file
# write.csv(ph_int,
#           paste0("npn-data/phenophases-intensities-", yr, ".csv"),
#           row.names = FALSE)

# Write dataframe with intensity values/midpoints to file
ia_yr <- ia %>%
  filter(intensity_category_id %in% unique(ph_int$intensity_category_id))
# write.csv(ia_yr,
#           paste0("npn-data/intensity-values-", yr, ".csv"),
#           row.names = FALSE)

#########################
# Few notes/questions on phenophase-intensity associations

select(ph_int, -nlevels)

# Two phenophase classes are missing from these phenophase-intensity summaries: 
# 5 (Falling leaves or needles) and 9 (End of flowering). 
# I think this is because there are no intensity categories associated with 
# phenophases in these classes. 

# Most phenophases associated with a single intensity category. 
# Exceptions have 2 intensity categories with different max values (1000, 10000)
  # phenophase 500 (Flowers and flower buds); intensity categories 48 & 49
  # phenophase 516 (Fruits present); intensity categories 56 & 57
  # phenophase 504 (Recent fruit or seed drop); intensity categories 59 & 60

# Most intensity categories associated with a single phenophase.
# Exceptions are the following:
  # category 39: Breaking leaf buds (10000); phenophases 371, 480, 481
  # category 38: Plant greenness; phenophases 489, 497, 509 
  # category 51: Pollen release; phenophases 502, 503

# Only one qualitative intensity category used: Pollen release.
# None of the "peak" intensity categories were used.
#########################

# SSPI data (from API) --------------------------------------------------------#

# Read in csv with species-phenophase-abundance/intensity category combinations
# Got this (SSPI data) from the NPN API see:
# OneDrive/NPN/IntensityMeasures/get-phenophase-intensity-pairs.R

spi <- read.csv("npn-data/spps-phenophases-intensities.csv")

head(select(spi, -additional_definition))
head(ia)
head(ia_cats)
head(ph_int)

spi <- spi %>%
  rename(intensity_category_id = abundance_category)

# Do all intensity categories in ia/ia_cats appear in SSPI data?
spi_cats <- unique(spi$intensity_category_id)
filter(ia_cats, !intensity_category_id %in% spi_cats)
# No. Missing the following in the SSPI data: 
# 18 (Young needle bundles present)
# 21 (Initial growth of buds or shoots)
# 22 (Initial growth of shoots)
filter(ph_int, intensity_category_id %in% c(18, 21, 22))
# These categories must have been phased out because they weren't used in 2024
# But why not in the table from the API, that should include categories that 
# have been deactivated?

# Do all intensity categories in ph_int (2024 combinations of phenophases and 
# intensity categories in status data) appear in SSPI data?
filter(ph_int, !intensity_category_id %in% spi_cats)
# Yes

# Merge SSPI data with info on intensity categories I created before:
spi2 <- spi %>%
  # Get rid of categories that had been deactivated
  filter(is.na(stopdate)) %>%
  # Get rid of spp-phenophase combinations that don't have intensity category
  filter(!is.na(intensity_category_id)) %>%
  # Merge with intensity category info
  left_join(ia_cats, by = "intensity_category_id")
count(spi2, nlevels) 
# Everything matches up - no intensity categories in SSPI data that I didn't
# already have from rnpn downloads


# If the following phenophases don't appear in current SSPI table, what have 
# they been replaced with? (18 looks like it was replaced with a similar 
# category and different max value)
# 21 (Initial growth of buds or shoots)
# 22 (Initial growth of shoots)

# Look for the following phenophases in SSPI table:
# 71 (Emergence above ground)
# 482/492/508 (Initial growth)
spi %>% 
  filter(phenophase_id %in% c(71, 482, 492, 508)) %>% 
  count(intensity_category_id)
# No intensity categories associated with these phenophases (always NA)
