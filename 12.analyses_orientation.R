#--------------------------------------------------------------------------#
# Project: rhiziis
# Script purpose: Analyse annotations from keypoint tool annotation.
# Date: 13/10/2022
# Author: Thelma Panaiotis
#--------------------------------------------------------------------------#


library(tidyverse)
library(arrow)
library(ggpmisc)
library(ggpubr)

source("lib/library_env.R")


## List validated objects ----
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

predicted_images <- read_parquet("data/isiis_images_env.parquet") %>% 
  filter(str_detect(transect, "cross")) %>% 
  filter(transect != "cross_current_01")



## Match annotations with metadata ----
#--------------------------------------------------------------------------#
# Read annotations
key_points <- read_csv("data/keypoint_annotations.csv")

# List files of processed images
img_files <- list.files(path = "data/images_taxo", pattern = "_seg.png", full.names = TRUE, recursive = TRUE) %>% 
  str_remove_all("_seg") %>% 
  unique()

# Add images without annotations to keypoints
key_points <- tibble(orig_id = str_split_fixed(img_files, "/", n = 5)[,5] %>% str_remove(".png")) %>% 
  left_join(key_points)

# Keep only validated solitary black with vacuole validated in ecotaxa
vac_val <- df %>% filter(classif_id_1 %in% c("with-vacuole<solitaryblack", "with-large-vacuole", "Arthracanthida", "Aulosphaeridae") & classif_qual == "validated")
key_points <- key_points %>% filter(orig_id %in% vac_val$orig_id)


## Compute orientation of the line joining the two points ----
#--------------------------------------------------------------------------#
key_points <- key_points %>% 
  # Reverse y axis so it goes from bottom to top (in images it is top to bottom)
  mutate(
    y1 = -y1,
    y2 = -y2
  ) %>% 
  drop_na(x1) %>% # Drop images without annotation
  mutate(angle = atan2(y2 - y1, x2 - x1)) # Compute angle from arctan

# Correct orientation by pitch
key_points <- key_points %>% 
  rename(orig_id = orig_id) %>% 
  left_join(df %>% select(orig_id, img_name)) %>% 
  left_join(predicted_images %>% select(img_name, depth, dist, pitch)) %>% 
  mutate(
    pitch = deg2rad(pitch),
    angle = angle + pitch
  )


## Compute length between points ----
#--------------------------------------------------------------------------#
key_points <- key_points %>% 
  mutate(
    length = sqrt((x2 - x1)^2 + (y2 - y1)^2),
    length = length * 51 / 1000 # convert from px to mm
  )


## Angle of solitary vacuoles ----
#--------------------------------------------------------------------------#
key_points %>% 
  filter(taxon == "solitary_vacuoles") %>% 
  ggplot() +
  geom_histogram(aes(x = angle), fill = "white", color = "black", breaks = seq(0, 2*pi, by = pi/50), boundary = 0) +
  scale_x_continuous(
    limits = c(-pi, pi),
    breaks = c(-pi/2, 0, pi/2, pi),
    labels = c(expression(frac(-pi, 2)), "0", expression(frac(pi, 2)), expression(pi))
  ) +
  coord_polar(start = pi/2, direction = -1) +
  theme_minimal() +
  labs(x = "")
ggsave(file = "figures/collodaria_vacuole_angle.pdf", width = 12, height = 12, unit = "cm", dpi = 300)


## Length of vacuoles ----
#--------------------------------------------------------------------------#
key_points %>% 
  filter(taxon == "solitary_vacuoles") %>% 
  ggplot(aes(x = length, y = -depth)) +
  geom_point(size = 0.5) +
  stat_poly_line(orientation = "y") +
  stat_cor(label.sep='\n', aes(label = paste(..rr.label.., ..p.label.., sep = "~")), p.accuracy = 0.001, r.accuracy = 0.01, label.x = 3) +
  scale_y_continuous(limits = c(-132, 0), expand = c(0, 0)) +
  labs(x = "Vacuole length (mm)", y = "Depth (m)") +
  theme_minimal()
ggsave(file = "figures/collodaria_vacuole_size.pdf", width = 16, height = 12, unit = "cm", dpi = 300)

# Compute a linear model of vacuole length VS depth
mod <- lm(length ~ depth, data = key_points)
summary(mod)
summary(mod)$r.squared


## Distribution of phaeodium angle ----
#--------------------------------------------------------------------------#
key_points %>% 
  filter(taxon == "aulosphaeridae") %>% 
  ggplot() +
  geom_histogram(aes(x = angle), fill = "white", color = "black", breaks = seq(-pi, pi, by = pi/50), boundary = 0) +
  #geom_density(aes(x = angle)) +
  scale_x_continuous(
    limits = c(-pi, pi),
    breaks = c(-pi/2, 0, pi/2, pi),
    labels = c(expression(frac(-pi, 2)), "0", expression(frac(pi, 2)), expression(pi))
  ) +
  coord_polar(start = pi/2, direction = -1) +
  theme_minimal() +
  labs(x = "")
ggsave(file = "figures/phaeodium_orientation.pdf", width = 12, height = 12, unit = "cm", dpi = 300)


## Distribution of Arthracanthida angle ----
#--------------------------------------------------------------------------#
key_points %>% 
  filter(taxon == "Arthracanthida") %>% 
  mutate(angle = angle + pi) %>% # Set angle between 0 and pi
  ggplot() +
  geom_histogram(aes(x = angle), fill = "white", color = "black", breaks = seq(-pi, pi, by = pi/50), boundary = 0) +
  scale_x_continuous(
    limits = c(-pi, pi),
    breaks = c(-pi/2, 0, pi/2, pi),
    labels = c(expression(frac(-pi, 2)), "0", expression(frac(pi, 2)), expression(pi))
  ) +
  coord_polar(start = pi/2, direction = -1) +
  theme_minimal() +
  labs(x = "")
ggsave(file = "figures/arthracanthida_angle.pdf", width = 12, height = 12, unit = "cm", dpi = 300)
