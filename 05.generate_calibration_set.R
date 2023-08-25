#--------------------------------------------------------------------------#
# Project: visiis_rhizaria
# Script purpose: Generate a calibration set for Rhizaria
# Date: 22/08/2022
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

library(tidyverse)
library(arrow)
library(lubridate)
library(reticulate)
library(googlesheets4)
gs4_auth(use_oob = TRUE)

source("lib/library_zoo.R")
source("lib/library_env.R")
source("lib/library_plot.R")

transect_plot <- "cross_current_02"

data_dir <- "data/distribution"
plot_dir_env <- "plot/env"
plot_dir_pl <- "plot/plankton"

parallel <- TRUE
n_cores <- 24

# Proportion of predictions to keep among best predictions
threshold <- 0.5

# Taxa regroupment
ss <- "https://docs.google.com/spreadsheets/d/100Hmfx9BVTFtBzw4pUETj6m3I4KjpibGvuRvDEVRYC8/edit?usp=sharing"


## Read files ----
#--------------------------------------------------------------------------#
# Read extracted objects and keep only validated
df <- read_parquet("data/all_rhizaria.parquet")

# Read predicted images
predicted_images <- read_parquet("data/isiis_images_env.parquet") %>% filter(str_detect(transect, "cross"))


## Select relevant taxa ----
#--------------------------------------------------------------------------#
# Build taxonomic count 
tc <- df %>% count(classif_id_1)

# Start by erasing previous data (3 first columns) in spreadsheet
range_flood(ss, sheet = "tc", range = "tc!A:B", reformat = FALSE)
# Write new tree
range_write(ss, data = tc) 
# Open it in browser tab to make sure everything is ok
gs4_browse(ss)

# Read taxo match
tc <- read_sheet(ss)

# Keep only selected taxa
df <- df %>% 
  left_join(tc) %>% 
  select(-classif_id_1) %>% 
  drop_na(new_taxon) %>% 
  select(orig_id, objid, path_to_img, classif_qual, cnn_score, img_name, n, taxon = new_taxon)


## List taxa to calibrate ----
#--------------------------------------------------------------------------#
taxa_calib <- tc %>% 
  drop_na(threshold_score) %>% 
  filter(threshold_score < 1) %>% 
  pull(new_taxon)
taxa_calib


## Generate calibration set ----
#--------------------------------------------------------------------------#
# For each taxa, randomly select n objects stratified by score decile
n_sample <- 200 # objects per class

# Initiate emtpy tibble
df_sample <- tibble()

for (my_taxon in taxa_calib){
  df_taxon <- df %>% 
    filter(taxon == my_taxon & classif_qual == "P") %>% 
    arrange(cnn_score) %>% 
    mutate(decile = cut(cnn_score, breaks = unique(quantile(cnn_score, probs = seq(0, 1, 0.1))), include.lowest = TRUE)) %>% 
    group_by(decile) %>% 
    slice_sample(n = n_sample / 10, replace = TRUE) %>% 
    ungroup()
  
  df_sample <- bind_rows(df_sample, df_taxon)
}

df_sample %>% count(taxon, decile)

df_sample <- df_sample %>% distinct(orig_id, .keep_all = TRUE)


## Save tsv table for ecotaxa ----
#--------------------------------------------------------------------------#
pd <- import("pandas")

df_sample <- df_sample %>% 
  mutate(
    object_annotation_status = "predicted",
    img_file_name = paste0("images/", img_name, "/", orig_id, ".png")
  ) %>% 
  select(
    object_id = orig_id,
    img_file_name,
    object_annotation_category = taxon,
    object_annotation_status,
    object_cnn_score = cnn_score
    )

# first_row as Dataframe row with appropriate headers
first_row = c('[t]', '[t]', '[t]', '[t]', '[f]')
first_row = t(pd$DataFrame(first_row))
colnames(first_row) <- colnames(df_sample)

# concat first_row and dataframe
df <- rbind(first_row, df_sample)

# save table
write_tsv(df, 'data/ecotaxa_rhizaria_calibration.tsv')

