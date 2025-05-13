# Exploring relationships between plant phenophases and intensity categories
# ER Zylstra
# 13 May 2025

library(dplyr)
library(stringr)
library(lubridate)
library(rnpn)
library(jsonlite)

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

# Intensity categories labled as "peak" (categories 31-34) are confusing 
# because they can have numeric or qualitative responses.
filter(ia, grepl("peak", intensity_name))
filter(ia_orig, grepl("peak", category_name))

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

#########################
# One note on intensity categories

# Categories 13 and 41 (Leaf size, in %) seem to be identical
filter(ia, intensity_category_id %in% c(13, 41))
  # I think category 13 was replaced by 41
#########################

# Summarize information about each intensity category
ia_cats <- ia %>%
  group_by(intensity_category_id, intensity_name, intensity_type) %>%
  summarize(nlevels = n(), .groups = "keep") %>%
  data.frame()

# IMPORTANT NOTE: ##############################################################
# Want to match up phenophases and intensity categories. Originally I did this
# using status intensity data (downloaded all plant data in 2024!). However, I
# later learned about the Secondary Phenophase Details table available in the 
# API that lists all species-phenophase-abundance/intensity categories (they
# call this the SSPI data). Originally I checked that these two sources of data
# provided matching information, which they did. Removed code that downloaded
# massive status-intensity dataset and just kept the SPI data after that
################################################################################

# Obtain SSPI data (species-phenophase-intensity combos) ----------------------#

json_data <- read_json(
  "https://services.usanpn.org/npn_portal/phenophases/getSecondaryPhenophaseDetails.json",
  simplifyVector = TRUE
)

colnames(json_data)
  # species_specific_phenophase_id # unique identifier
  # phenophase_id             # 149 
  # species_id                # 1887
  # additional_definition 
  # active                    # 0/1
  # effective_datetime        # First date this was effective
  # deactivation_datetime     # Last date this was effective ("" = still active)
  # multimedia_uri 
  # multimedia_credit
  # abundance_category        # Two-digit ID or ""
  # extent_min                # Not sure what this is (1 or "")
  # extent_max                # Not sure what this is (1000000 or "")

spi <- json_data %>%
  select(-c(multimedia_uri, multimedia_credit, extent_min, extent_max)) %>%
  rename(sspi = species_specific_phenophase_id,
         startdate = effective_datetime,
         stopdate = deactivation_datetime) %>%
  mutate(startdate = as_date(startdate),
         stopdate = as_date(stopdate),
         abundance_category = ifelse(abundance_category == "", 
                                     NA, abundance_category)) %>%
  rename(intensity_category_id = abundance_category)
# Will leave the additional_definition column in there for now, though it's 
# probably unnecessary

# Write to file
# write.csv(spi, "npn-data/spps-phenophases-intensities.csv", row.names = FALSE)

# Identify unique combinations of phenophase and intensity category -----------#

# Remove species-phenophase combinations that don't have intensity categories or
# have intensity categories that were not used in 2024-2025
spi <- spi %>%
  select(-c(additional_definition, sspi)) %>%
  mutate(stopyear = year(stopdate)) %>%
  filter(!is.na(intensity_category_id)) %>%
  filter(stopyear >= 2024 | is.na(stopyear)) %>%
  select(-stopyear)

# Group SPI data
ph_int <- spi %>%
  group_by(phenophase_id, intensity_category_id) %>%
  summarize(n_spp = n_distinct(species_id), .groups = "keep") %>%
  mutate(intensity_category_id = as.integer(intensity_category_id)) %>%
  data.frame()

# Add information from other datasets
ph_int <- ph_int %>%
  left_join(ph, by = "phenophase_id") %>%
  left_join(ia_cats, by = "intensity_category_id") %>%
  # Remove any phenophases not associated with plants
  filter(!is.na(class_id)) %>%
  relocate(c(phenophase_name, class_id, class_name, category), 
           .after = "phenophase_id")  %>%
  relocate(c(intensity_name, intensity_type, nlevels),
           .after = "intensity_category_id") %>%
  rename(intensity_levels = nlevels) %>%
  # Sort
  arrange(class_id, phenophase_id)

# Write to file
# write.csv(ph_int, 
#           "npn-data/phenophases-intensities-2024.csv", 
#           row.names = FALSE)


# Write dataframe with intensity values/midpoints to file (just those categories
# used in 2024-2025)
ia_yr <- ia %>%
  filter(intensity_category_id %in% unique(ph_int$intensity_category_id))

# write.csv(ia_yr,
#           "npn-data/intensity-values-2024.csv",
#           row.names = FALSE)


#########################
# Few notes/questions on phenophase-intensity associations

select(ph_int, -c(intensity_levels, n_spp))

# One phenophase class is missing from these phenophase-intensity summaries: 
# 5 (Falling leaves or needles). I think this is because there are no intensity 
# categories associated with phenophases in this class. 

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
