#--------------------------------------------------------------------------#
# Project: visiis_copepoda
# Script purpose: Prepare a learning set from validated images to train a CNN classifier
# Date: 24/05/2022
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

library(tidyverse)
library(arrow)
library(yaml)


## Read files ----
#--------------------------------------------------------------------------#
# Read configuration file
cfg <- read_yaml("config.yaml")

# Read extracted objects and keep only validated
df <- read_parquet("data/all_rhizaria.parquet") %>% 
  filter(classif_qual == "V") %>% 
  select(-classif_qual)


## Ignore classes with too few objects ----
#--------------------------------------------------------------------------#
df <- df %>% 
  group_by(classif_id_1) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n > cfg$grouping$n_min) %>% 
  select(-n)


## Flag plankton classes ----
#--------------------------------------------------------------------------#
# List non plankton classes
non_plankton <- c("other<living", "detritus")
# Set plankton to TRUE for plankton classes
df <- df %>% mutate(plankton = ifelse(classif_id_1 %in% non_plankton, FALSE, TRUE))
# Quick check
df %>% select(classif_id_1, plankton) %>% unique()


## Split training and validation set ----
#--------------------------------------------------------------------------#
df <- df %>% 
  group_by(classif_id_1) %>% 
  # shuffle rows
  sample_frac(1) %>%
  mutate(
    # percent rank (which is random)
    r=(1:n()) / n(),
    # assign in set based in this
    set=case_when(
      r >= 1 - cfg$split_props$valid ~ "valid",
      TRUE ~ "train"
    )
  ) %>% 
  select(-r) %>%
  ungroup() %>% 
  relocate(set, .after = plankton)


## Fix path to img ----
#--------------------------------------------------------------------------#
df <- df %>% 
  mutate(path_to_img = str_c("data/images/", path_to_img))


## Save table ----
#--------------------------------------------------------------------------#
write_csv(df, file = file.path(cfg$base_dir, "rhizaria_learning_set.csv"))
