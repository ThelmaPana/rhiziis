#--------------------------------------------------------------------------#
# Project: rhiziis
# Script purpose: Apply confidence threshold to predictions
# Date: 18/09/2022
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

library(tidyverse)
library(googlesheets4)
library(arrow)

gs4_auth(use_oob = TRUE)
2

overwrite_tree <- FALSE


## Read files ----
#--------------------------------------------------------------------------#

# Read export from ecotaxa
df <- read_delim(
  "data/ecotaxa_export_all_rhizaria_cc.tsv", 
  delim = "\t", 
  escape_double = FALSE, 
  trim_ws = TRUE
) %>% # select and rename columns
  select(
    orig_id = object_id,
    objid,
    classif_id_1 = object_annotation_category,
    cnn_score = object_cnn_score,
    classif_qual = object_annotation_status,
    area = object_area,
    eccentricity = object_eccentricity,
    esd = object_equivalent_diameter,
    orientation = object_orientation,
    img_name = acq_id
  )


# Read thresholds
thr <- read_csv("data_git/07.taxa_threshold_rhizaria.csv") %>% select(taxon = taxa, score_threshold)


## Rename taxa ----
#--------------------------------------------------------------------------#
# Overwrite taxonomic sheet if required
if (overwrite_tree){
  # Google sheet for taxonomy match
  ss <- "https://docs.google.com/spreadsheets/d/100Hmfx9BVTFtBzw4pUETj6m3I4KjpibGvuRvDEVRYC8/edit?usp=sharing"
  # Build taxonomic count 
  tc <- df %>% count(classif_id_1)
  # Start by erasing previous data (3 first columns) in spreadsheet
  range_flood(ss, sheet = "tc", range = "tc!A:B", reformat = FALSE)
  # Write new tree
  range_write(ss, data = tc) 
  # Open it in browser tab to make sure everything is ok
  gs4_browse(ss)
  # Read taxo match
  tc <- read_sheet(ss) %>% select(classif_id_1, new_taxon, group, use)
} else { # otherwise, read the CSV file containing taxonomy match
  tc <- read_csv("data_git/taxonomy_match.csv")
}

# Keep only selected taxa
df <- df %>% 
  left_join(tc) %>% 
  select(-classif_id_1) %>% 
  drop_na(new_taxon) %>% 
  rename(taxon = new_taxon) %>% 
  filter(taxon != "detritus") %>% 
  filter(use == 1) %>% 
  select(-use)


## Threshold objects ----
#--------------------------------------------------------------------------#
# Keep all validated objects
df_v <- df %>% filter(classif_qual == "validated") %>% select(img_name, orig_id, taxon, group, area:orientation)

# Filter predicted objects based on CNN score
df_p <- df %>% 
  filter(classif_qual == "predicted") %>% 
  left_join(thr) %>% 
  drop_na(score_threshold) %>% 
  filter(cnn_score >= score_threshold) %>% 
  select(img_name, orig_id, taxon, group, area:orientation)


# Put back together validated and predicted objects
df <- bind_rows(df_v, df_p)
taxa <- unique(df$taxon)


## Cleaning ----
#--------------------------------------------------------------------------#
# Order taxa and groups
df <- df %>% 
  mutate(
    taxon = factor(taxon, levels = rev(taxa)),
    group = factor(group, levels = c("Acantharea", "Phaeodaria", "Collodaria (solitary)", "Collodaria (colonial)"))
  )


## Save objects ----
#--------------------------------------------------------------------------#
write_parquet(df, sink = "data/08.thresholded_rhizaria.parquet")
