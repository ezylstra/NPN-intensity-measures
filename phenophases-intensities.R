# Exploring relationships between plant phenophases and intensity categories
# ER Zylstra
# 27 Feb 2025

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

ia_orig <- npn_abundance_categories() 
  # Returns tibble with category ids, each with a dataframe that lists the value
  # ids and value names (which is what appears in intensity_value column in a
  # status-intensity dataset)

# Remove blank entries 
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

# Create new dataframe, that identifies whether values are a number/count, 
# percent, or if they're qualitative. If not qualitative, extract bounding 
# values to calculate an appoximate midpoint.
ia <- ia_orig %>%
  select(-value_description) %>%
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

# Download status-intensity data ----------------------------------------------#

# Want to match up phenophases and intensity categories, but I'm not sure how to
# do this without using a status-intensity (si) dataset.

# For now, will download recent si data, knowing that it won't include certain 
# phenophases and/or intensity categories that have been phased out (which might
# be fine)

# Commenting this out since I don't want to accidentally start a new download 
# or overwrite the existing file

# # Limit download to plants by finding relevant taxonomic orders first
#   all_spp <- npn_species()
#   plant_spp <- all_spp %>%
#     filter(kingdom == "Plantae") %>%
#     select(species_id, common_name, functional_type, class_id, class_name)
#   plant_classes <- unique(plant_spp$class_id)
# 
# # Download status-intensity data (this takes a long time! millions of records)
#   si_orig <- npn_download_status_data(
#     request_source = "erinz",
#     years = year(Sys.Date()) - 1,
#     class_ids = plant_classes,
#     climate_data = FALSE,
#     additional_fields = c("site_name",
#                           "species_functional_type")
#   )
# 
# # Remove unnecessary columns 
#   si_orig <- si_orig %>%
#     select(-c(update_datetime, site_name, elevation_in_meters, genus, species,
#               kingdom, day_of_year, abundance_value))
# 
# # Don't need all the observations. Just want a table that provides combinations
# # of species, phenophase, and intensity_category_id (when provided). 
#   si_combos <- si_orig %>%
#     filter(intensity_category_id != -9999) %>%
#     distinct(species_id, common_name, species_functional_type, phenophase_id,
#              phenophase_description, intensity_category_id, intensity_value)
# # Write this to file
#   write.csv(si_combos,
#             paste0("npn-data/si-plants-", year(Sys.Date()) - 1,
#                    "-spp-ph-int-combos.csv"),
#             row.names = FALSE)

# Load si data and combine with phenophase, intensity information -------------#

# Read in the species-phenophase-intensity category combination file
si_combos <- read.csv(paste0("npn-data/si-plants-", year(Sys.Date()) - 1, 
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

# Create dataframe with unique combinations of phenophase, intensity categories
# that were used in recent year of observations
ph_int <- si_combos %>%
  group_by(phenophase_id, phenophase_description, intensity_category_id) %>%
  summarize(n_spp = n_distinct(common_name), .groups = "keep") %>%
  left_join(select(ph, phenophase_id, category, class_id, class_name),
            by = "phenophase_id") %>%
  left_join(ia_cats, by = "intensity_category_id") %>%
  data.frame() %>%
  arrange(class_id, phenophase_id)
select(ph_int, -nlevels)

#########################
# Few notes/questions on phenophase-intensity associations

select(ph_int, -nlevels)

# No observations in two phenophase classes: 5 (Falling leaves or needles) and 
# 9 (End of flowering)

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
# None of the "peak" categories were used.
#########################
