#--------------------------------------------------------------------------#
# Project: rhiziis
# Script purpose: Process annotations from keypoint tool annotation.
# Date: 13/10/2022
# Author: Thelma Panaiotis
#--------------------------------------------------------------------------#


library(tidyverse)
library(rjson)


## Read json files containing annotations ----
#--------------------------------------------------------------------------#
# Outputs of keypoint annotation tool are json files.
# List json files of annotations
json_files <- list.files(path = "data/images_taxo", pattern = ".json", full.names = TRUE, recursive = TRUE)

# Initiate empty tibble to store annotations
key_points <- tibble()

# Loop over json files and store points in a dataframe
for (file in json_files){
  orig_id <- str_split_fixed(file, "/", n = 5)[,5] %>% str_remove(".json")
  taxon <- str_split_fixed(file, "/", n = 5)[,3]
  res <- fromJSON(file = file) %>% 
    as.data.frame() %>% 
    mutate(
      orig_id = orig_id,
      taxon = taxon
      )
  
  key_points <- key_points %>% bind_rows(res)
}


## Rename columns and save data ----
#--------------------------------------------------------------------------#
key_points <- key_points %>% select(orig_id, taxon, x1 = tooltips.x, y1 = tooltips.y, x2 = tooltips.x.1, y2 = tooltips.y.1)
write_csv(key_points, file = "data/08.keypoint_annotations.csv")
