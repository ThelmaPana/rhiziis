#--------------------------------------------------------------------------#
# Project: rhiziis
# Script purpose: Compute Rhizaria concentration and interpolate transects
# Date: 18/09/2022
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

library(tidyverse)
library(arrow)
#library(lubridate)
#library(cmocean)
#library(castr)
library(googlesheets4)

gs4_auth(use_oob = TRUE)
2

source("lib/library_zoo.R")
source("lib/library_env.R")
source("lib/library_plot.R")

# Settings for plots
transect_plot <- "cross_current_02" # Choose a transect for plots
test_plot <- FALSE # whether to plot intermediate stuff

plot_dir_env <- "plot/env"
plot_dir_pl <- "plot/plankton"

# Settings for interpolation
parallel <- TRUE
n_cores <- 24

# Taxa regroupment
ss <- "https://docs.google.com/spreadsheets/d/100Hmfx9BVTFtBzw4pUETj6m3I4KjpibGvuRvDEVRYC8/edit?usp=sharing"


## Read files ----
#--------------------------------------------------------------------------#
# Extracted objects
df <- read_parquet("data/08.thresholded_rhizaria.parquet")

# Metadata from images
predicted_images <- read_parquet("data/isiis_images_env.parquet") %>% 
  filter(str_detect(transect, "cross")) %>% 
  filter(transect != "cross_current_01")

# Make a list of taxa
taxa <- sort(unique(df$taxon))


## Join objects with metadata from images ----
#--------------------------------------------------------------------------#
df <- df %>% 
  left_join(predicted_images) %>% 
  select(transect, img_name, yo, yo_type, period, datetime, lat, lon, dist, depth, taxon, group, area:orientation, pitch)


## Additional corrections ----
#--------------------------------------------------------------------------#
# Keep only objects above 100 m
df <- df %>% filter(depth < 100)

# Convert esd from px to mm
# And correct esd for ISIIS speed
df <- df %>% 
  mutate(
    esd_mm = esd * 51 / 1000,
    esd_mm = esd_mm * sqrt(1.428)
  ) %>% 
  select(-esd)

# Correct for the pitch of ISIIS
df <- df %>% 
  mutate(
    pitch = deg2rad(pitch), # pitch from degree to radian
    orientation = orientation + pitch # Add pitch to orientation
  )

# Plot ISIIS pitch
pitch_lim <- max(abs(df$pitch))
df %>% 
  filter(transect == "cross_current_02") %>% 
  ggplot() +
  geom_point(aes(x = dist, y = -depth, color = pitch)) +
  scale_color_distiller(palette = "RdYlBu", limits = c(-pitch_lim, pitch_lim))

# Save corrected objects
write_parquet(df, sink = "data/09.corr_rhizaria.parquet")


## Concentration per image ----
#--------------------------------------------------------------------------#
df_img <- df %>% 
  group_by(img_name, taxon) %>% 
  summarise(abund = n()) %>% 
  ungroup()

# Generate all combination of image name and taxon
df_conc <- predicted_images %>% 
  filter(str_detect(transect, "cross")) %>% 
  crossing(taxon = taxa) %>% 
  left_join(df_img) %>% 
  mutate(
    conc = abund / image_vol,
    conc = ifelse(is.na(conc), 0, conc),
    conc = ifelse(keep, conc, NA)
    ) %>% 
  select(-c(abund, image_vol, seg_threshold, keep)) %>% 
  pivot_wider(names_from = taxon, values_from = conc)


## Concentration per meter bin ----
#--------------------------------------------------------------------------#
plankton_bin <- df_conc %>% 
  mutate(depth = round(depth)) %>% 
  select(transect, yo, yo_type, period, depth, everything()) %>% 
  group_by(transect, yo, yo_type, period, depth) %>% 
  summarise(across(.cols = c(datetime, lat:heading, all_of(taxa)), ~ mean(.x, na.rm = TRUE))) %>% 
  ungroup()

# Save
write_parquet(plankton_bin, sink = "data/09.plankton_bin_semantic.parquet")


## Coarse interpolation of plankton concentration ----
#--------------------------------------------------------------------------#
# Interpolate plankton concentrations and env data
# interpolation step in X axis: 1 km
# interpolation step in Y axis: 1 m
plankton_int_coarse <- coarse_interp(df = plankton_bin, vars = plankton_bin %>% select(temp:dens, all_of(taxa)) %>% colnames(), step_x = 1, step_y = 1, parallel = parallel, n_cores = n_cores)  

# Interpolate datetime from dist
plankton_int_coarse <- interp_datetime(df_bin = plankton_bin, df_int = plankton_int_coarse)

# Plot datetime VS dist
plankton_int_coarse %>% 
  ggplot() +
  geom_point(aes(x = dist, y = datetime)) +
  facet_wrap(~transect, scales = "free") +
  ggtitle("Datetime VS distance")


## Ignore interpolated data too far from original points ----
#--------------------------------------------------------------------------#
# NB: run this only once as it takes some time
if (!file.exists("data/09.bins_keep.parquet")){ # If file does not exist yet, run it!
  nn_x <- 1 # distance to which look for nearest neighbour in x (dist, km), default is 1
  nn_y <- 2 # distance to which look for nearest neighbour in y (depth, m), default is 2
  max_dist <- 2 # distance threshold above which interpolated data is deleted
  
  # Parallel cleaning
  if (parallel){
    intl <- mclapply(unique(plankton_int_coarse$transect), function(id) {
      
      # filter interpolated data for transect
      df_int <- plankton_int_coarse %>% filter(transect == id)
      
      # filter binned plankton data for transect
      df_bin  <- plankton_bin %>% 
        filter(transect == id) %>% 
        # and keep only rows where plankton is present
        mutate(plankton = rowSums(across(taxa))) %>% 
        drop_na(plankton)
      
      # delete unwanted values
      df_int <- clean_zoo_interp(df_int, df_bin, vars=taxa, nn_x=nn_x, nn_y=nn_y, max_dist=max_dist)
      
      return(df_int)
      
    }, mc.cores=n_cores) 
    # this returns a list, recombine it into a tibble
    plankton_int_coarse <- do.call(bind_rows, intl)
    
  } else{ # not parallel cleaning
    
    # Initiate empty tibble for interpolation results
    plankton_int_coarse_cl <- tibble()
    
    # Loop over transects
    for(id in unique(plankton_int_coarse$transect)) {
      
      # filter interpolated data for transect
      df_int <- plankton_int_coarse %>% filter(transect == id)
      
      # filter binned plankton data for transect
      df_bin  <- plankton_bin %>% filter(transect == id)
      
      # delete unwanted values
      df_int <- clean_zoo_interp(df_int, df_bin, vars=taxa, nn_x=nn_x, nn_y=nn_y, max_dist=max_dist)
      
      # Append transect rows to other transects
      plankton_int_coarse_cl <- bind_rows(plankton_int_coarse_cl, df_int)
      plankton_int_coarse <- plankton_int_coarse_cl
      del(plankton_int_coarse_cl)
    }
  }
  
  # Save output
  bins_keep <- plankton_int_coarse %>% 
    mutate(to_keep = ifelse(is.na(Acantharea), FALSE, TRUE)) %>% 
    select(transect, dist, depth, to_keep)
  write_parquet(bins_keep, sink = "data/09.bins_keep.parquet")

} else{ # Otherwise, just read the file
  bins_keep <- read_parquet(file = "data/09.bins_keep.parquet")
}

# Set plankton conc to NA in bins not to be kept
plankton_int_coarse <- plankton_int_coarse %>% 
  left_join(bins_keep) %>% 
  mutate(across(all_of(taxa), ~ ifelse(to_keep, ., NA))) %>% 
  select(-to_keep)

if (test_plot){
  # Plot Aulacanthidae
  plankton_int_coarse %>% 
    #filter(transect == transect_plot) %>% 
    ggplot() +
    geom_tile(aes(x = dist, y = -depth, fill = Aulacanthidae)) +
    scale_fill_viridis_c(trans = "log1p", na.value = NA) +
    theme_minimal() +
    facet_wrap(~transect, scales = "free_x")
  
  # Plot salinity
  plankton_int_coarse %>% 
    filter(transect == transect_plot) %>%
    #filter(transect == "along_current_04") %>% 
    ggplot() +
    geom_tile(aes(x = dist, y = -depth, fill = sal)) +
    scale_fill_cmocean(name = "haline", na.value = NA) +
    labs(fill = "Salinity<br>(ppt)") +
    theme_light() +
    theme(legend.title = element_markdown())
}

## Plot all coarse interpolation
plot_interp_conc(plankton_int_coarse, taxa, output_file = file.path(plot_dir_pl, paste0("09.plankton_conc_int_coarse_semantic.pdf")))
plot_interp_env(plankton_int_coarse, output_file = file.path(plot_dir_env, "09.env_int_coarse.pdf"))

# Save interpolated data
write_parquet(plankton_int_coarse, sink = "data/09.plankton_int_coarse_semantic.parquet")


## Fine interpolate plankton distribution ----
#--------------------------------------------------------------------------#
## Fine interpolation
step_x <- 0.2 # fine interpolation step in X axis: 0.2 km or 200 m
step_y <- 0.5 # fine interpolation step in Y axis: 0.5 m
theta <- 0.5 # bandwidth or scale parameter

plankton_int_fine <- fine_interp(
  df = plankton_int_coarse, 
  vars = plankton_bin %>% select(temp:dens, all_of(taxa)) %>% colnames(), 
  step_x = 0.2, # fine interpolation step in X axis: 0.2 km or 200 m
  step_y =  0.5, # fine interpolation step in Y axis: 0.5 m
  theta = 0.5, # bandwidth or scale parameter
  parallel = TRUE, 
  n_cores = n_cores
)

# Interpolate datetime from dist
plankton_int_fine <- interp_datetime(df_bin = plankton_bin, df_int = plankton_int_fine)

if (test_plot){
  # Plot datetime VS dist
  plankton_int_fine %>% 
    group_by(transect) %>% 
    sample_n(100) %>% 
    ungroup() %>% 
    ggplot() +
    geom_point(aes(x = dist, y = datetime)) +
    facet_wrap(~transect, scales = "free") +
    ggtitle("Datetime VS distance")
  
  plankton_int_fine %>% 
    ggplot() +
    geom_raster(aes(x = dist, y = -depth, fill = Acantharea)) +
    facet_wrap(~transect, scales = "free") +
    scale_fill_viridis_c(trans = "log1p", na.value = NA) 
  
  plankton_int_fine %>% 
    filter(transect == transect_plot) %>% 
    pivot_longer(Acantharea:Collodaria_solitaryglobule, names_to = "taxon", values_to = "conc") %>% 
    ggplot() +
    geom_tile(aes(x = dist, y = -depth, fill = conc)) +
    geom_contour(aes(x = dist, y = -depth, z = fluo), breaks = c(0.05, 0.1, 0.15), size = 0.5, colour = "gray") +
    facet_wrap(~taxon, scales = "free") +
    scale_fill_viridis_c(trans = "log1p", na.value = NA) 
  
  plankton_int_fine %>% 
    filter(transect == transect_plot) %>% 
    ggplot() +
    geom_tile(aes(x = dist, y = -depth, fill = fluo)) +
    scale_fill_cmocean(name = "algae") +
    geom_contour(aes(x = dist, y = -depth, z = fluo), breaks = c(0.05, 0.1, 0.15), size = 0.5, colour = "gray") 
  }

# Plot all coarse interpolation
plot_interp_conc(plankton_int_fine, taxa, output_file = file.path(plot_dir_pl, paste0("09.plankton_conc_int_fine_semantic.pdf")))
plot_interp_env(plankton_int_fine, output_file = file.path(plot_dir_env, "09.env_int_fine.pdf"))

# Save interpolated data
write_parquet(plankton_int_fine, sink = "data/09.plankton_int_fine_semantic.parquet")

