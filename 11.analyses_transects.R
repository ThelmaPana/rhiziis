#--------------------------------------------------------------------------#
# Project: rhiziis
# Script purpose: Perform analyses
# Date: 18/09/2022
# Author: Thelma Panaiotis
#--------------------------------------------------------------------------#

library(tidyverse)
library(arrow)
library(cmocean)
library(castr)
library(egg)
library(patchwork)
library(ggridges)

transect_plot <- "cross_current_02"

# Function to generate limits for plots
lims <- function(x){c(min(x, na.rm = TRUE), max(x, na.rm = TRUE))}


## Read data ----
#--------------------------------------------------------------------------#
# Extracted objects
df <- read_parquet("data/09.corr_rhizaria.parquet")
length(unique(df$taxon))


# Interpolated data (fine-scale)
plankton_int_fine <- read_parquet("data/09.plankton_int_fine_semantic.parquet")


## Dataset composition ----
#--------------------------------------------------------------------------#
df_comp <- df %>%
  count(group, taxon) %>% 
  arrange(desc(n)) %>% 
  mutate(taxon = factor(taxon, levels = taxon))

df_comp %>% 
  ggplot() +
  geom_col(aes(x = taxon, y = n, fill = group)) + 
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

df_comp


## ESD distribution by group ----
#--------------------------------------------------------------------------#
# Plot ESD per group as well as number of objects per group
df %>% 
  group_by(taxon) %>% 
  mutate(
    n = n(),
    median = median(esd_mm)
  ) %>% 
  ungroup() %>% 
  arrange(median) %>% 
  mutate(taxon = factor(taxon, levels = rev(unique(taxon)))) %>% 
  #add_count(taxon) %>% 
  ggplot(aes(y = taxon, x = esd_mm)) +
  geom_violin(aes(fill = n), scale = "width", alpha = 0.5, linewidth = 0.2) +
  geom_boxplot(width = 0.2, outlier.colour = NA, linewidth = 0.2) +
  scale_x_continuous(trans = "log1p", breaks = c(0.1, 0.5, 1, 5, 10, 20, 40), limits = c(0.1, 40), expand = c(0,0)) +
  #annotation_logticks(sides = "b", scaled = FALSE) +
  scale_fill_viridis_c(trans = "log1p", breaks = c(100, 1000, 10000, 100000), limits = c(100, max(df$n))) +
  #scale_fill_hp(option = "slytherin", trans = "log1p", breaks = c(1000, 10000, 100000)) +
  facet_grid(group ~ ., scales = "free", space = "free") +
  theme_minimal() +
  labs(x = "ESD (mm)", y = "Taxonomic group") +
  theme(legend.key.height = unit(1.5, 'cm'))
ggsave(file = "figures/esd_distribution.pdf", width = 20, height = 12, unit = "cm", dpi = 300)


## Position of DCM ----
#--------------------------------------------------------------------------#
# Compute the position of the DCM as the depth of maximum fluorescence (distance-wise)
dcm <- plankton_int_fine %>% 
  filter(transect != "merged") %>% 
  filter(depth > 40) %>% 
  group_by(transect, dist) %>% 
  summarise(
    n = sum(!is.na(fluo)),
    dcm = maxd(fluo, depth = depth)
  ) %>% 
  filter(n > 100) %>% 
  ungroup() %>% 
  select(-n)

# Smooth the DCM with a moving average
dcm_sm <- dcm %>% 
  group_by(transect) %>% 
  mutate(
    dcm = despike(dcm, k = 5),
    dcm = smooth(dcm, k = 5, n = 2)
  ) %>% 
  ungroup()

# Plot the DCM
dcm_sm %>% 
  ggplot() +
  geom_path(aes(x = dist, y = -dcm, group = transect)) +
  facet_wrap(~transect, scales = "free_x")

# On top of chla
plankton_int_fine %>% 
  filter(transect != "merged") %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = fluo)) +
  geom_path(aes(x = dist, y = -dcm, group = transect), data = dcm_sm) +
  scale_fill_cmocean(name = "algae", na.value = NA) +
  theme_minimal() +
  facet_wrap(~transect, ncol = 2)

# On top of oxygen
plankton_int_fine %>% 
  filter(transect != "merged") %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = oxy)) +
  geom_path(aes(x = dist, y = -dcm, group = transect), data = dcm_sm) +
  scale_fill_distiller(palette = "Blues", direction = 1, na.value = NA) +
  theme_minimal() +
  facet_wrap(~transect, ncol = 2)



## Vertical distance to of organisms DCM ----
#--------------------------------------------------------------------------#
# For each organism, compute its distance to the DCM
plankton_dcm <- df %>% 
  mutate(dist = roundp(dist, precision = 0.2)) %>%
  left_join(dcm) %>% 
  drop_na(dcm) %>% 
  mutate(
    dcm_dist = dcm - depth
  ) 

# Plot distance to DCM for all taxa, day VS night
plankton_dcm %>% 
  mutate(
    dn = ifelse(as.integer(str_split_fixed(transect, "_", n = 3)[,3]) %in% c(2, 4, 7), "night", "day"),
    depth_d = ifelse(dn == "day", dcm_dist, NA),
    depth_n = ifelse(dn == "night", dcm_dist, NA),
    taxon = str_replace_all(taxon, "_", "\n")
  ) %>% 
  ggplot(aes(fill = dn, group = dn)) +
  geom_density(aes(y = depth_d, x = -after_stat(density)), outline.type = "both", alpha = 0.5) +
  geom_density(aes(y = depth_n, x = after_stat(density)), outline.type = "both", alpha = 0.5) +
  geom_hline(yintercept = 0, alpha = 0.2, linewidth = 2, color = "darkgreen") +
  scale_fill_manual(values = c("white", "gray")) +
  labs(x = "Taxon", y = "Distance to DCM (m)", fill = "", title = "Distance to DCM") +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(), 
    axis.text.x = element_blank(), 
    panel.spacing.x = unit(0, "null"), 
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x.bottom = element_text(angle = 90, vjust = 1)
  ) +
  facet_wrap(~taxon, scales = "free_x", nrow = 1, strip.position = "bottom") 
# No differences between day and night -> merge together

# Collodaria
plankton_dcm %>% 
  mutate(
    dn = ifelse(as.integer(str_split_fixed(transect, "_", n = 3)[,3]) %in% c(2, 4, 7), "night", "day"),
    taxon = str_replace_all(taxon, "_", "\n")
  ) %>% 
  filter(str_detect(taxon, "Collodaria")) %>% 
  ggplot() +
  geom_violin(aes(x = taxon, y = dcm_dist), fill = "gray95", scale = "width") +
  geom_hline(yintercept = 0, alpha = 0.2, linewidth = 2, color = "darkgreen") +
  labs(x = "Taxon", y = "Distance to DCM (m)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
  ) +
  facet_grid(.~group, scales = "free", space = "free") 
ggsave(file = "figures/collodaria_dcm_dist.pdf", width = 16, height = 12, unit = "cm", dpi = 300)


# Acantharia
plankton_dcm %>% 
  mutate(
    dn = ifelse(as.integer(str_split_fixed(transect, "_", n = 3)[,3]) %in% c(2, 4, 7), "night", "day"),
    taxon = str_replace_all(taxon, "_", "\n")
  ) %>% 
  filter(group == "Acantharea") %>% 
  ggplot() +
  geom_violin(aes(x = taxon, y = dcm_dist), fill = "gray95", scale = "width") +
  geom_hline(yintercept = 0, alpha = 0.2, linewidth = 2, color = "darkgreen") +
  labs(x = "Taxon", y = "Distance to DCM (m)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
  ) +
  facet_grid(.~group, scales = "free", space = "free") 
ggsave(file = "figures/acantharea_dcm_dist.pdf", width = 8, height = 12, unit = "cm", dpi = 300)


# Phaeodaria
plankton_dcm %>% 
  mutate(
    dn = ifelse(as.integer(str_split_fixed(transect, "_", n = 3)[,3]) %in% c(2, 4, 7), "night", "day"),
    taxon = str_replace_all(taxon, "_", "\n")
  ) %>% 
  filter(group == "Phaeodaria") %>% 
  ggplot() +
  geom_violin(aes(x = taxon, y = dcm_dist), linewidth = "gray95", scale = "width") +
  geom_hline(yintercept = 0, alpha = 0.2, size = 2, color = "darkgreen") +
  labs(x = "Taxon", y = "Distance to DCM (m)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
  ) +
  facet_grid(.~group, scales = "free", space = "free") 
ggsave(file = "figures/phaeodaria_dcm_dist.pdf", width = 8, height = 12, unit = "cm", dpi = 300)


## Plot interpolated environment ----
#--------------------------------------------------------------------------#
# Temp
p1 <- plankton_int_fine %>% 
  filter(transect == transect_plot) %>% 
  filter(dist < 59) %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = temp)) +
  scale_fill_cmocean(name = "thermal", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + 
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = "Temperature\n(°C)", tag = "A") +
  theme_minimal() +
  theme(text = element_text(size = 8)) + 
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 4))

# Salinity
p2 <- plankton_int_fine %>% 
  filter(transect == transect_plot) %>% 
  filter(dist < 59) %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = sal)) +
  geom_contour(aes(x = dist, y = -depth, z = sal), breaks = c(38.2, 38.3), colour = "white", linewidth = 0.3) +
  scale_fill_cmocean(name = "haline", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + 
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = "Salinity", tag = "B") +
  theme_minimal() +
  theme(text = element_text(size = 8)) + 
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 4))

# Dens
p3 <- plankton_int_fine %>% 
  filter(transect == transect_plot) %>% 
  filter(dist < 59) %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = dens)) +
  scale_fill_cmocean(name = "dense", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + 
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = "Density\n(kg m⁻³)", tag = "C") +
  theme_minimal() +
  theme(text = element_text(size = 8)) + 
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 4))

# Fluo
p4 <- plankton_int_fine %>% 
  filter(transect == transect_plot) %>% 
  filter(dist < 59) %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = fluo)) +
  geom_path(aes(x = dist, y = -dcm), linewidth = 0.3, color = "white", data = dcm %>% filter(transect == transect_plot & dist < 59)) +
  scale_fill_cmocean(name = "algae", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + 
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = "Fluorescence\n(V)", tag = "D") +
  theme_minimal() +
  theme(text = element_text(size = 8)) + 
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 4))

# Oxy
p5 <- plankton_int_fine %>% 
  filter(transect == transect_plot) %>% 
  filter(dist < 59) %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = oxy)) +
  scale_fill_distiller(palette = "Blues", direction = 1, na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + 
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = "Oxygen\n(µmol kg⁻¹)", tag = "E") +
  theme_minimal() +
  theme(text = element_text(size = 8)) + 
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 4))


p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 2)
ggsave(file = "figures/interpolated_env.pdf", width = 178, height = 150, unit = "mm", dpi = 300)


## Aulacanthidae and environment ----
#--------------------------------------------------------------------------#

## Aulacanthidae, flat Aulacanthidae and oxygen
# With oxygen
p1 <- plankton_int_fine %>% 
  filter(transect == transect_plot) %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = Aulacanthidae)) +
  geom_contour(aes(x = dist, y = -depth, z = oxy), breaks = c(210, 230, 250), colour = "#3d87c0", linewidth = 0.3) +
  geom_path(aes(x = dist, y = -dcm), colour = "#b2df8a", linewidth = 0.3, data = dcm_sm %>% filter(transect == "cross_current_02")) +
  scale_fill_gradient(high = "white", low = "black", trans = "log1p", na.value = NA, breaks = c(0, 10, 40)) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = "Aulac.\n(m⁻³)", tag = "A") +
  theme_minimal() +
  theme(text = element_text(size = 8)) +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 4))


p2 <- plankton_int_fine %>% 
  filter(transect == transect_plot) %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = Aulacanthidae_flat)) +
  geom_contour(aes(x = dist, y = -depth, z = oxy), breaks = c(210, 230, 250), colour = "#3d87c0", linewidth = 0.3) +
  geom_path(aes(x = dist, y = -dcm), colour = "#b2df8a", linewidth = 0.3, data = dcm_sm %>% filter(transect == "cross_current_02")) +
  scale_fill_gradient(high = "white", low = "black", trans = "log1p", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = "Flat Aulac.\n(m⁻³)", tag = "B") +
  theme_minimal() +
  theme(text = element_text(size = 8)) +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 4))


p3 <- plankton_int_fine %>% 
  filter(transect == transect_plot) %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = oxy)) +
  geom_contour(aes(x = dist, y = -depth, z = oxy), breaks = c(210, 230, 250), colour = "#eff3ff", size = 0.3) +
  geom_path(aes(x = dist, y = -dcm), colour = "#b2df8a", linewidth = 0.3, data = dcm_sm %>% filter(transect == "cross_current_02")) +
  scale_fill_distiller(palette = "Blues", direction = 1, na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = "Oxygen\n(µmol kg⁻¹)", tag = "C") +
  theme_minimal() +
  theme(text = element_text(size = 8)) +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 4))

p <- p1 + p2 + p3  + plot_layout(ncol = 2)
p
ggsave(plot = p, file = "figures/aulacanthidae_transect.png", width = 178, height = 100, unit = "mm", dpi = 300)


## Aulacanthidae, downwelling waters and DCM ----
#--------------------------------------------------------------------------#
# Smoothed DCM with the same precision as plankton data
dcm_sm <- dcm %>% 
  group_by(transect) %>% 
  mutate(
    dcm = despike(dcm, k = 5),
    dcm = smooth(dcm, k = 5, n = 2),
    dcm = roundp(dcm, precision = 0.5) # round to 0.5 to match with data
  ) %>% 
  ungroup()

# Generate short names for transects
transect_names <- tribble(
  ~transect, ~short_name,
  "cross_current_02", "2",
  "cross_current_03", "3",
  "cross_current_04", "4",
  "cross_current_05", "5",
  "cross_current_06", "6",
  "cross_current_07", "7",
)

my_taxon <- "Aulacanthidae"

# Get Aulacanthidae and oxygen along the DCM
dcm_conc <- dcm_sm %>% 
  rename(depth = dcm) %>% 
  left_join(plankton_int_fine %>% select(transect:depth, oxy, all_of(my_taxon))) %>% 
  pivot_longer(all_of(my_taxon), names_to = "taxon", values_to = "conc") %>% 
  # Normalise values between 0 and 1 for each transect
  group_by(transect) %>% 
  mutate(
    oxy = (oxy - min(oxy, na.rm = TRUE)) / (max(oxy, na.rm = TRUE) - min(oxy, na.rm = TRUE)),
    conc = (conc - min(conc, na.rm = TRUE)) / (max(conc, na.rm = TRUE) - min(conc, na.rm = TRUE))
  ) %>% 
  ungroup()

# Get distance of beginning and end of the transects to plot the DCM
line_dist <- dcm_conc %>% 
  drop_na(oxy) %>% 
  group_by(transect) %>% 
  summarise(
    min_dist = min(dist),
    max_dist = max(dist),
  ) %>% 
  left_join(transect_names)

# Concentrations (Aulacanthidae and oxygen) along the DCM
dcm_conc %>% 
  left_join(transect_names) %>% 
  select(-transect) %>% 
  select(transect = short_name, everything()) %>% 
  mutate(
    transect = factor(transect, levels = rev(unique(transect)))
  ) %>% 
  mutate(conc = -conc) %>% 
  pivot_longer(c(oxy, conc)) %>% 
  mutate(group = str_c(transect, name, sep = "_")) %>% 
  ggplot() +
  geom_ridgeline(aes(x = dist, y = transect, group = group, height = value, color = name, fill = name), min_height = -1, scale = 0.5, alpha = 0.5) +
  scale_color_manual(
    breaks = c("oxy", "conc"),
    values = c("oxy" = "#3d87c0", "conc" = "gray50"),
    labels = c("Oxygen", "Aulac.")
  ) +
  scale_fill_manual(
    breaks = c("oxy", "conc"),
    values = c("oxy" = "#3d87c0", "conc" = "gray90"),
    labels = c("Oxygen", "Aulac.")
  ) + 
  geom_segment(aes(y = short_name, yend = short_name, x = min_dist, xend = max_dist), color = "#b2df8a", data = line_dist) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Transect", fill = "Relative \nconcentration", colour = "Relative \nconcentration") +
  theme_minimal() +
  theme(text = element_text(size = 8))
ggsave(file = "figures/aulacanthidae_dcm.pdf", width = 160, height = 80, unit = "mm", dpi = 300)

# Transect
plankton_int_fine %>% 
  filter(transect == "cross_current_02") %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = Aulacanthidae)) +
  geom_contour(aes(x = dist, y = -depth, z = oxy), breaks = c(210, 230, 250), colour = "#3d87c0", linewidth = 0.3) +
  geom_path(aes(x = dist, y = -dcm), colour = "#b2df8a", linewidth = 0.3, data = dcm_sm %>% filter(transect == "cross_current_02")) +
  scale_fill_gradient(high = "white", low = "black", trans = "log1p", na.value = NA, breaks = c(0, 10, 40)) +
  #scale_fill_viridis_c(trans = "log1p", na.value = NA, breaks = c(0, 10, 40), option = "G") +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = expression(paste("Aulac. (", m^{-3}, ")"))) +
  theme_minimal() +
  theme(text = element_text(size = 8)) + 
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 4))
ggsave(file = "figures/aulacanthidae_conc.pdf", width = 90, height = 45, unit = "mm", dpi = 300)



## Organisms not affected by recirculation ----
#--------------------------------------------------------------------------#
# Maps
p1 <- plankton_int_fine %>% 
  filter(transect == transect_plot) %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = Arthracanthida)) +
  geom_contour(aes(x = dist, y = -depth, z = oxy), breaks = c(210, 230, 250), colour = "#3d87c0", linewidth = 0.5) +
  geom_path(aes(x = dist, y = -dcm), colour = "#b2df8a", linewidth = 0.5, data = dcm_sm %>% filter(transect == "cross_current_02")) +
  scale_fill_gradient(high = "white", low = "black", trans = "log1p", na.value = NA, limits = lims(plankton_int_fine$Arthracanthida), breaks = c(0, 5, 20)) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = expression(paste("Conc. (", m^{-3}, ")")), tag = "A", title = "Arthracanthida") +
  theme_minimal() +
  theme(text = element_text(size = 8)) +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 4))

p2 <- plankton_int_fine %>% 
  filter(transect == transect_plot) %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = Collodaria_solitaryglobule)) +
  geom_contour(aes(x = dist, y = -depth, z = oxy), breaks = c(210, 230, 250), colour = "#3d87c0", linewidth = 0.5) +
  geom_path(aes(x = dist, y = -dcm), colour = "#b2df8a", linewidth = 0.5, data = dcm_sm %>% filter(transect == "cross_current_02")) +
  scale_fill_gradient(high = "white", low = "black", trans = "log1p", na.value = NA, limits = lims(plankton_int_fine$Arthracanthida), breaks = c(0, 5, 20)) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = expression(paste("Conc. (", m^{-3}, ")")), tag = "B", title = "Collodaria solitary") +
  theme_minimal() +
  theme(text = element_text(size = 8), axis.title.y = element_blank(), axis.text.y = element_blank()) +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 4))

p <- p1 + p2 + plot_layout(guides = "collect")

ggsave(plot = p, file = "figures/not_affected_sinking_map.pdf", width = 164, height = 60, unit = "mm", dpi = 300)

# Along DCM
my_taxon <- "Arthracanthida"
dcm_conc <- dcm_sm %>% 
  rename(depth = dcm) %>% 
  left_join(plankton_int_fine %>% select(transect:depth, oxy, all_of(my_taxon))) %>% 
  pivot_longer(all_of(my_taxon), names_to = "taxon", values_to = "conc") %>% 
  # Normalise values between 0 and 1 for each transect
  group_by(transect) %>% 
  mutate(
    oxy = (oxy - min(oxy, na.rm = TRUE)) / (max(oxy, na.rm = TRUE) - min(oxy, na.rm = TRUE)),
    conc = (conc - min(conc, na.rm = TRUE)) / (max(conc, na.rm = TRUE) - min(conc, na.rm = TRUE))
  ) %>% 
  ungroup()

line_dist <- dcm_conc %>% 
  drop_na(oxy) %>% 
  group_by(transect) %>% 
  summarise(
    min_dist = min(dist),
    max_dist = max(dist),
  ) %>% 
  left_join(transect_names)


dcm_conc %>% 
  left_join(transect_names) %>% 
  select(-transect) %>% 
  select(transect = short_name, everything()) %>% 
  mutate(
    transect = factor(transect, levels = rev(unique(transect)))
  ) %>% 
  #filter(isopyc == "iso_28.7") %>% 
  mutate(conc = -conc) %>% 
  pivot_longer(c(oxy, conc)) %>% 
  mutate(group = str_c(transect, name, sep = "_")) %>% 
  ggplot() +
  geom_ridgeline(aes(x = dist, y = transect, group = group, height = value, color = name, fill = name), min_height = -1, scale = 0.5, alpha = 0.5) +
  scale_color_manual(
    breaks = c("oxy", "conc"),
    values = c("oxy" = "#3d87c0", "conc" = "gray50"),
    labels = c("Oxygen", "Arthra.")
  ) +
  #guides(colour = "none") +
  scale_fill_manual(
    breaks = c("oxy", "conc"),
    values = c("oxy" = "#3d87c0", "conc" = "gray90"),
    labels = c("Oxygen", "Arthra.")
  ) + 
  geom_segment(aes(y = short_name, yend = short_name, x = min_dist, xend = max_dist), color = "#b2df8a", data = line_dist) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Transect", fill = "Relative \nconcentration", colour = "Relative \nconcentration") +
  theme_minimal() +
  theme(text = element_text(size = 8))
ggsave(file = "figures/arthracanthida_dcm.pdf", width = 140, height = 60, unit = "mm", dpi = 300)


my_taxon <- "Collodaria_solitaryglobule"
dcm_conc <- dcm_sm %>% 
  rename(depth = dcm) %>% 
  left_join(plankton_int_fine %>% select(transect:depth, oxy, all_of(my_taxon))) %>% 
  pivot_longer(all_of(my_taxon), names_to = "taxon", values_to = "conc") %>% 
  # Normalise values between 0 and 1 for each transect
  group_by(transect) %>% 
  mutate(
    oxy = (oxy - min(oxy, na.rm = TRUE)) / (max(oxy, na.rm = TRUE) - min(oxy, na.rm = TRUE)),
    conc = (conc - min(conc, na.rm = TRUE)) / (max(conc, na.rm = TRUE) - min(conc, na.rm = TRUE))
  ) %>% 
  ungroup()

line_dist <- dcm_conc %>% 
  drop_na(oxy) %>% 
  group_by(transect) %>% 
  summarise(
    min_dist = min(dist),
    max_dist = max(dist),
  ) %>% 
  left_join(transect_names)

dcm_conc %>% 
  left_join(transect_names) %>% 
  select(-transect) %>% 
  select(transect = short_name, everything()) %>% 
  mutate(
    transect = factor(transect, levels = rev(unique(transect)))
  ) %>% 
  #filter(isopyc == "iso_28.7") %>% 
  mutate(conc = -conc) %>% 
  pivot_longer(c(oxy, conc)) %>% 
  mutate(group = str_c(transect, name, sep = "_")) %>% 
  ggplot() +
  geom_ridgeline(aes(x = dist, y = transect, group = group, height = value, color = name, fill = name), min_height = -1, scale = 0.5, alpha = 0.5) +
  scale_color_manual(
    breaks = c("oxy", "conc"),
    values = c("oxy" = "#3d87c0", "conc" = "gray50"),
    labels = c("Oxygen", "Coll. sol.")
  ) +
  #guides(colour = "none") +
  scale_fill_manual(
    breaks = c("oxy", "conc"),
    values = c("oxy" = "#3d87c0", "conc" = "gray90"),
    labels = c("Oxygen", "Coll. sol.")
  ) + 
  geom_segment(aes(y = short_name, yend = short_name, x = min_dist, xend = max_dist), color = "#b2df8a", data = line_dist) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Transect", fill = "Relative \nconcentration", colour = "Relative \nconcentration") +
  theme_minimal() +
  theme(text = element_text(size = 8))
ggsave(file = "figures/collodaria_solitary_dcm.pdf", width = 140, height = 60, unit = "mm", dpi = 300)


## Additional illustration figure ----
#--------------------------------------------------------------------------#
# Relative conc between 0 and 1 for all groups
plankton_scale <- plankton_int_fine %>% 
  filter(transect == transect_plot) %>% 
  pivot_longer(Acantharea:Collodaria_colonial_budding, names_to = "taxon", values_to = "conc") %>% 
  group_by(taxon) %>%
  mutate(conc = (conc / max(conc, na.rm = TRUE))) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "taxon", values_from = "conc")

p1 <- plankton_scale %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = Arthracanthida)) +
  scale_fill_gradient(high = "white", low = "black", trans = "log1p", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = expression(paste("Conc. (", m^{-3}, ")")), title = "Arthracanthida") +
  theme_minimal() +
  theme(text = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank())

p2 <- plankton_scale %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = Acantharea_small)) +
  scale_fill_gradient(high = "white", low = "black", trans = "log1p", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = expression(paste("Conc. (", m^{-3}, ")")), title = "small Acantharia") +
  theme_minimal() +
  theme(text = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())

p3 <- plankton_scale %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = Aulacanthidae)) +
  scale_fill_gradient(high = "white", low = "black", trans = "log1p", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = expression(paste("Conc. (", m^{-3}, ")")), title = "Aulacanthidae") +
  theme_minimal() +
  theme(text = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())

p4 <- plankton_scale %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = oxy)) +
  scale_fill_distiller(palette = "Blues", direction = 1, na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = bquote("µmol"~kg^-1), title = "Oxygen") +
  theme_minimal() +
  theme(text = element_text(size = 8), axis.text.y = element_blank(), axis.title.y = element_blank())

p5 <- plankton_scale %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = Collodaria_solitaryglobule)) +
  scale_fill_gradient(high = "white", low = "black", trans = "log1p", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = expression(paste("Conc. (", m^{-3}, ")")), title = "Coll. sol.") +
  theme_minimal() +
  theme(text = element_text(size = 8))

p6 <- plankton_scale %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = `Collodaria_solitaryblack_with-vacuole`)) +
  scale_fill_gradient(high = "white", low = "black", trans = "log1p", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = expression(paste("Conc. (", m^{-3}, ")")), title = "Coll. sol. vac.") +
  theme_minimal() +
  theme(text = element_text(size = 8), axis.text.y = element_blank(), axis.title.y = element_blank())

p7 <- plankton_scale %>% 
  ggplot() +
  geom_raster(aes(x = dist, y = -depth, fill = Collodaria_colonial)) +
  scale_fill_gradient(high = "white", low = "black", trans = "log1p", na.value = NA) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  labs(x = "Distance from shore (km)", y = "Depth (m)", fill = expression(paste("Conc. (", m^{-3}, ")")), title = "Coll. colonial") +
  theme_minimal() +
  theme(text = element_text(size = 8), axis.text.y = element_blank(), axis.title.y = element_blank())

p <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + guide_area() + plot_layout(nrow = 2, guides = 'collect')
ggsave(plot = p, file = "figures/illust_dist.pdf", width = 180, height = 100, unit = "mm", dpi = 300)

