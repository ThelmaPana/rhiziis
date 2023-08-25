#--------------------------------------------------------------------------#
# Project: VISIIS
# Script purpose: Regroup functions used for environemental data processing of visufront
# Date: 02/12/2020
# Author: Thelma Panaiotis
#--------------------------------------------------------------------------#

library(cmocean)
library(egg)
library(ggtext)
library(grid)
library(gridExtra)
library(qpdf)
library(scales)


## Plot points for temperature, salinity, density, fluo, oxy and irrad for each transect and save as pdf ----
#--------------------------------------------------------------------------#
plot_raw_env <- function(df, output_file, step=1){
  #' Make a combined pdf of temp, sal, density, fluo, oxy and irrad variables for each transect.
  #'
  #' The pdf has one page per transect, 6 plots per page, using the cmocean color map. 
  #' @param df Dataframe containing ISIIS data to plot
  #' @param output_file name of output file
  #' @param step step for plotted point (plot one point every <step> point)
  
  # Create output_dir
  output_dir <- tools::file_path_sans_ext(output_file)
  dir.create(output_dir)
  
  # Empty output_dir
  files <- list.files(output_dir, full.names = TRUE)
  file.remove(files)
  
  # Make list of transects
  transects <- df %>% pull(transect) %>% unique()
  
  # Loop over transects
  for (my_transect in transects){
    
    # Nice string for plot title
    my_transect_title <- str_replace_all(my_transect, "_", " ") %>% str_to_title()
    
    # Filter data for this transect and subsample rows
    d <- df %>% 
      filter(transect == my_transect) %>% 
      filter(row_number() %% step == 1)
    
    # Make plots
    p1 <- d %>% 
      ggplot() +
      geom_point(aes(x = dist, y = -depth, color = temp)) +
      scale_color_cmocean(name = "thermal") +
      labs(title = "Temperature") +
      labs(x = "Dist (km)", y = "Depth (m)", color = "Temperature<br>(°C)", title = "Temperature") +
      scale_x_continuous(expand = c(0,0))  + scale_y_continuous(expand = c(0,0)) +
      theme_minimal() +
      theme(legend.title = element_markdown())
    
    p2 <- d %>% 
      ggplot() +
      geom_point(aes(x = dist, y = -depth, color = sal)) +
      scale_color_cmocean(name = "haline") +
      labs(title = "Salinity") +
      labs(x = "Dist (km)", y = "Depth (m)", color = "Salinity<br>(ppt)", title = "Salinity") +
      scale_x_continuous(expand = c(0,0))  + scale_y_continuous(expand = c(0,0)) +
      theme_minimal() +
      theme(legend.title = element_markdown())
    
    p3 <- d %>% 
      ggplot() +
      geom_point(aes(x = dist, y = -depth, color = dens)) +
      scale_color_cmocean(name = "dense") +
      labs(title = "Density") +
      labs(x = "Dist (km)", y = "Depth (m)", color = "Potential Density<br>Anomaly (kg.m<sup>-3</sup> )", title = "Density") +
      scale_x_continuous(expand = c(0,0))  + scale_y_continuous(expand = c(0,0)) +
      theme_minimal() +
      theme(legend.title = element_markdown())
    
    p4 <- d %>% 
      ggplot() +
      geom_point(aes(x = dist, y = -depth, color = fluo)) +
      scale_color_cmocean(name = "algae") +
      labs(title = "Fluorescence") +
      labs(x = "Dist (km)", y = "Depth (m)", color = "Fluorescence<br>(V)", title = "Fluorescence") +
      scale_x_continuous(expand = c(0,0))  + scale_y_continuous(expand = c(0,0)) +
      theme_minimal() +
      theme(legend.title = element_markdown())
    
    p5 <- d %>% 
      ggplot() +
      geom_point(aes(x = dist, y = -depth, color = oxy)) +
      #scale_color_cmocean(name = "oxy") +
      scale_color_distiller(palette = "Blues", na.value=NA, direction = 1) +
      labs(title = "Oxygen") +
      labs(x = "Dist (km)", y = "Depth (m)", color = "Oxygen<br>(µmol.kg<sup>-1</sup> )", title = "Oxygen") +
      scale_x_continuous(expand = c(0,0))  + scale_y_continuous(expand = c(0,0)) +
      theme_minimal() +
      theme(legend.title = element_markdown())
    
    p6 <- d %>% 
      ggplot() +
      geom_point(aes(x = dist, y = -depth, color = irrad)) +
      scale_color_cmocean(name = "solar") +
      labs(title = "Irradiance") +
      labs(x = "Dist (km)", y = "Depth (m)", color = "Irradiance<br>(UE.cm<sup>-2</sup> )", title = "Irradiance") +
      scale_x_continuous(expand = c(0,0))  + scale_y_continuous(expand = c(0,0)) +
      theme_minimal() +
      theme(legend.title = element_markdown())
    
    # Arrange them on one page and save pdf with name of transect
    g <- ggarrange(p1, p2, p3, p4, p5, p6, ncol=2, top = textGrob(my_transect_title, gp=gpar(fontsize=15,font=8)), draw=FALSE)
    ggsave(g, file=str_c(output_dir, "/", my_transect, ".pdf"), width=15, height=15)
  }
  
  # List of newly created plots 
  my_plots <- list.files(output_dir, full.names = TRUE)
  
  # Combine all pages in one pdf
  pdf_combine(input = my_plots, output = str_c(output_dir, "/all_transects.pdf"))
  
  # Delete single page plots
  file.remove(my_plots)
  
  # Move multi transect file out of temporary dir
  file.rename(from = str_c(output_dir, "/all_transects.pdf"), to = output_file)
  
  # Delete temporary dir
  unlink(output_dir, recursive = TRUE)
}


## Plot interpolated temperature, salinity, density, fluo and oxy for each transect and save as pdf ----
#--------------------------------------------------------------------------#
plot_interp_env <- function(df, output_file){
  #' Make a combined pdf of plots of interpolated temp, sal, dens, fluo and oxy variables for each transect.
  #'
  #' The pdf has one page per transect, 5 plots per page, using the cmocean color map. 
  #' @param df Dataframe containing interpolated ISIIS data to plot
  #' @param output_file name of output file
  
  
  # Create output_dir
  output_dir <- tools::file_path_sans_ext(output_file)
  dir.create(output_dir)
  
  # Empty output_dir
  files <- list.files(output_dir, full.names = TRUE)
  file.remove(files)
  
  # Make list of transects
  transects <- df %>% pull(transect) %>% unique()
  
  # Loop over transects
  for (my_transect in transects){
    
    # Nice string for plot title
    my_transect_title <- str_replace_all(my_transect, "_", " ") %>% str_to_title()
    
    # Filter data for this transect
    d <- df %>% filter(transect == my_transect)
    
    # Make plots
    p1 <- d %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = temp)) +
      scale_fill_cmocean(name = "thermal", na.value=NA) +
      labs(x = "Dist (km)", y = "Depth (m)", fill = "Temperature<br>(°C)", title = "Temperature") +
      theme_light() +
      theme(legend.title = element_markdown())
    
    p2 <- d %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = sal)) +
      scale_fill_cmocean(name = "haline", na.value=NA) +
      labs(x = "Dist (km)", y = "Depth (m)", fill = "Salinity<br>(ppt)", title = "Salinity") +
      geom_contour(aes(x = dist, y = -depth, z = sal), breaks=c(38.2,38.3), color="white") +
      theme_light() +
      theme(legend.title = element_markdown())
    
    p3 <- d %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = dens)) +
      scale_fill_cmocean(name = "dense", na.value=NA) +
      labs(x = "Dist (km)", y = "Depth (m)", fill = "Potential Density<br>Anomaly (kg.m<sup>-3</sup> )", title = "Density") +
      theme_light() +
      theme(legend.title = element_markdown())
    
    p4 <- d %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = fluo)) +
      scale_fill_cmocean(name = "algae", na.value=NA) +
      labs(x = "Dist (km)", y = "Depth (m)", fill = "Fluorescence<br>(V)", title = "Fluorescence") +
      theme_light() +
      theme(legend.title = element_markdown())
    
    p5 <- d %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = oxy)) +
      scale_fill_distiller(palette = "Blues", na.value=NA, direction = 1) +
      #scale_fill_distiller(palette = "Blues", na.value=NA, direction = 1) +
      labs(x = "Dist (km)", y = "Depth (m)", fill = "Oxygen<br>(µmol.kg<sup>-1</sup> )", title = "Oxygen") +
      theme_light() +
      theme(legend.title = element_markdown())
    
    p6 <- d %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = irrad)) +
      scale_fill_cmocean(name = "solar", na.value=NA) +
      labs(x = "Dist (km)", y = "Depth (m)", fill = "Irradiance<br>(UE.cm<sup>-2</sup> )", title = "Irradiance") +
      theme_light() +
      theme(legend.title = element_markdown())
    
    # Arrange them on one page and save pdf with name of transect
    g <- ggarrange(p1, p2, p3, p4, p5, p6, ncol=2, top = textGrob(my_transect_title, gp=gpar(fontsize=15,font=8)), draw=FALSE)
    ggsave(g, file=str_c(output_dir, "/", my_transect, ".pdf"), width=15, height=15)
  }
  
  # List of newly created plots 
  my_plots <- list.files(output_dir, full.names = TRUE)
  
  # Combine all pages in one pdf
  pdf_combine(input = my_plots, output = str_c(output_dir, "/all_transects.pdf"))
  
  # Delete single page plots
  file.remove(my_plots)
  
  # Move multi transect file out of temporary dir
  file.rename(from = str_c(output_dir, "/all_transects.pdf"), to = output_file)
  
  # Delete temporary dir
  unlink(output_dir, recursive = TRUE)
}


## Plot interpolated temperature, salinity, density, fluo and oxy for each transect alongside depth_wise anomaly and save as pdf ----
#--------------------------------------------------------------------------#
plot_interp_env_anom <- function(df, output_file){
  #' Make a combined pdf of plots of interpolated temp, sal, dens, fluo and oxy variables and their depth-wise anomaly for each transect.
  #'
  #' The pdf has one page per transect, 10 plots per page, using the cmocean color map for variables and a diverging colormap for anomalies.
  #' @param df Dataframe containing interpolated ISIIS data to plot
  #' @param output_file name of output file
  
  
  # Create output_dir
  output_dir <- tools::file_path_sans_ext(output_file)
  dir.create(output_dir)
  
  # Empty output_dir
  files <- list.files(output_dir, full.names = TRUE)
  file.remove(files)
  
  # Make list of transects
  transects <- df %>% pull(transect) %>% unique()
  
  # Loop over transects
  for (my_transect in transects){
    
    # Nice string for plot title
    my_transect_title <- str_replace_all(my_transect, "_", " ") %>% str_to_title()
    
    # Filter data for this transect
    d <- df %>% filter(transect == my_transect)
    
    # Make plots
    p1 <- d %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = temp)) +
      scale_fill_cmocean(name = "thermal", na.value = NA) +
      labs(title = "Temperature<br>(°C)") +
      theme_light() +
      theme(legend.title = element_markdown())
    
    p2 <- d %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = temp_anom)) +
      scale_fill_gradient2(high = muted("red"), low = muted("blue"), na.value=NA) +
      labs(title = "Temperature anomaly") +
      theme_light() +
      theme(legend.title = element_markdown())
    
    p3 <- d %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = sal)) +
      scale_fill_cmocean(name = "haline", na.value=NA) +
      labs(title = "Salinity<br>(ppt)") +
      theme_light() +
      theme(legend.title = element_markdown())
    
    p4 <- d %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = sal_anom)) +
      scale_fill_gradient2(high = muted("red"), low = muted("blue"), na.value=NA) +
      labs(title = "Salinity anomaly") +
      theme_light() +
      theme(legend.title = element_markdown())
    
    p5 <- d %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = dens)) +
      scale_fill_cmocean(name = "dense", na.value=NA) +
      labs(title = "Potential Density<br>Anomaly (kg.m<sup>-3</sup> )") +
      theme_light() +
      theme(legend.title = element_markdown())
    
    p6 <- d %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = dens_anom)) +
      scale_fill_gradient2(high = muted("red"), low = muted("blue"), na.value=NA) +
      labs(title = "Density anomaly") +
      theme_light() +
      theme(legend.title = element_markdown())
    
    p7 <- d %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = fluo)) +
      scale_fill_cmocean(name = "algae", na.value=NA) +
      labs(title = "Fluorescence<br>(V)") +
      theme_light() +
      theme(legend.title = element_markdown())
    
    p8<- d %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = fluo_anom)) +
      scale_fill_gradient2(high = muted("red"), low = muted("blue"), na.value=NA) +
      labs(title = "Fluorescence anomaly") +
      theme_light() +
      theme(legend.title = element_markdown())
    
    p9 <- d %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = oxy)) +
      scale_fill_distiller(palette = "Blues", na.value=NA, direction = 1) +
      labs(title = "Oxygen<br>(µmol.kg<sup>-1</sup> )") +
      theme_light() +
      theme(legend.title = element_markdown())
    
    p10 <- d %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = oxy_anom)) +
      scale_fill_gradient2(high = muted("red"), low = muted("blue"), na.value=NA) +
      labs(title = "Oxygen anomaly") +
      theme_light() +
      theme(legend.title = element_markdown())
    
    # Arrange them on one page and save pdf with name of transect
    g <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, ncol=2, top = textGrob(my_transect_title, gp=gpar(fontsize=15,font=8)))
    ggsave(g, file=str_c(output_dir, "/", my_transect, ".pdf"), width=15, height=30)
  }
  
  # List of newly created plots 
  my_plots <- list.files(output_dir, full.names = TRUE)
  
  # Combine all pages in one pdf
  pdf_combine(input = my_plots, output = str_c(output_dir, "/all_transects.pdf"))
  
  # Delete single page plots
  file.remove(my_plots)
  
  # Move multi transect file out of temporary dir
  file.rename(from = str_c(output_dir, "/all_transects.pdf"), to = output_file)
  
  # Delete temporary dir
  unlink(output_dir, recursive = TRUE)
}


## Plot interpolated oxy data after lag correction for 6 lag values between 0 and 5 for each transect ----
#--------------------------------------------------------------------------#
plot_oxy_lags <- function(df, output_dir){
  #' Make a combined pdf of plots of oxy data after lag correction for 6 lag values between 0 and 5 for each transect
  #'
  #' The pdf has one page per transect, 6 plots per page, using the cmocean color map. 
  #' @param df Dataframe containing ISIIS oxy data to plot
  #' @param output_dir directory to save pdf of plots
  

  # Create output_dir
  dir.create(output_dir)
  
  # Empty output_dir
  files <- list.files(output_dir, full.names = TRUE)
  file.remove(files)
  
  # Make list of transects
  transects <- df %>% pull(transect) %>% unique()
  
  # Loop over transects
  for (my_transect in transects){
    
    # Nice string for plot title
    my_transect_title <- str_replace_all(my_transect, "_", " ") %>% str_to_title()
    
    # Filter data for this transect
    d <- df %>% filter(transect == my_transect)
    
    # Make plots
    p1 <- d %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = lag_0)) +
      scale_fill_distiller(palette = "Blues", na.value=NA, direction = 1) +
      labs(title = str_c(my_transect, " - lag_0")) +
      theme_light()
    
    p2 <- d %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = lag_1)) +
      scale_fill_distiller(palette = "Blues", na.value=NA, direction = 1) +
      labs(title = str_c(my_transect, " - lag_1")) +
      theme_light()
    
    p3 <- d %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = lag_2)) +
      scale_fill_distiller(palette = "Blues", na.value=NA, direction = 1) +
      labs(title = str_c(my_transect, " - lag_2")) +
      theme_light()
    
    p4 <- d %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = lag_3)) +
      scale_fill_distiller(palette = "Blues", na.value=NA, direction = 1) +
      labs(title = str_c(my_transect, " - lag_3")) +
      theme_light()
    
    p5 <- d %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = lag_4)) +
      scale_fill_distiller(palette = "Blues", na.value=NA, direction = 1) +
      labs(title = str_c(my_transect, " - lag_4")) +
      theme_light()
    
    p6 <- d %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = lag_5)) +
      scale_fill_distiller(palette = "Blues", na.value=NA, direction = 1) +
      labs(title = str_c(my_transect, " - lag_5")) +
      theme_light()
    
    # Arrange them on one page and save pdf with name of transect
    g <- arrangeGrob(p1, p2, p3, p4, p5, p6, ncol=2)
    ggsave(g, file=str_c(output_dir, "/", my_transect, ".pdf"), width=15, height=15)
  }
  # List of newly created plots 
  my_plots <- list.files(output_dir, full.names = TRUE)
  
  # Combine all pages in one pdf
  pdf_combine(input = my_plots, output = str_c(output_dir, "/all_transects.pdf"))
  
  # Delete single page plots
  file.remove(my_plots)
}


## Plot concentration per image (or depth bin) and per taxon and save plot ----
#--------------------------------------------------------------------------#
plot_raw_conc <- function(df, taxa, output_file){
  #' Plot values of concentration per image or depth bin per taxon along a transect and save plot.
  #'
  #' Values are plotted as -depth = f(dist) with concentration in color. Plots are saved in a multipage pdf. 
  #' @param df dataframe with concentration data to plot
  #' @param taxa list of taxa to plot
  #' @param output_file name of saved file 
  
  
  # Create output_dir
  output_dir <- tools::file_path_sans_ext(output_file)
  dir.create(output_dir)
  
  # Empty output_dir in case it already exists
  files <- list.files(output_dir, full.names = TRUE)
  file.remove(files)
  
  for (my_transect in unique(df$transect)){
    # Nice string for plot title
    my_transect_title <- str_replace_all(my_transect, "_", " ") %>% str_to_title()
    
    # Get data for transect
    d <- df %>% 
      filter(transect == my_transect) %>% 
      select(dist, depth, all_of(taxa)) %>% 
      gather(all_of(taxa), key = taxon, value = conc)
    
    # Plot 
    p <- d %>% 
      ggplot() +
      geom_point(aes(x = dist, y = -depth, color = conc)) +
      scale_color_viridis_c(trans = "log1p") +
      facet_wrap(~taxon, ncol = 4) +
      labs(x = "Distance (km)", y = "Depth (m)", color = "Concentration<br>(ind.m<sup>-3</sup> )") +
      theme_minimal() +
      theme(legend.title = element_markdown(), text = element_text(size = 14)) +
      ggtitle(my_transect_title)
    
    # Save
    ggsave(p, file=str_c(output_dir, "/", my_transect, ".pdf"), width=15, height=15)
    
  }    
  # List of newly created plots 
  my_plots <- list.files(output_dir, full.names = TRUE)
  
  # Combine all pages in one pdf
  pdf_combine(input = my_plots, output = str_c(output_dir, "/all_transects.pdf"))
  
  # Delete single page plots
  file.remove(my_plots)
  
  # Move multi transect file out of temporary dir
  file.rename(from = str_c(output_dir, "/all_transects.pdf"), to = output_file)
  
  # Delete temporary dir
  unlink(output_dir, recursive = TRUE)
}


## Plot interpolated values of concentration  per taxon and save plot ----
#--------------------------------------------------------------------------#
plot_interp_conc <- function(df, taxa, output_file){
  #' Plot interpolated values of concentration per taxon along a transect and save plot.
  #'
  #' Values are plotted as -depth = f(dist) with concentration in color. Plots are saved in a multipage pdf. 
  #' @param df dataframe with interpolated concentration data to plot
  #' @param taxa list of taxa to plot 
  #' @param output_file name of saved file 
  
  
  # Create output_dir
  output_dir <- tools::file_path_sans_ext(output_file)
  dir.create(output_dir)
  
  # Empty output_dir in case it already exists
  files <- list.files(output_dir, full.names = TRUE)
  file.remove(files)
  
  for (my_transect in unique(df$transect)){
    # Nice string for plot title
    my_transect_title <- str_replace_all(my_transect, "_", " ") %>% str_to_title()
    
    # Get data for this transect
    d <- df %>% 
      filter(transect == my_transect) %>% 
      select(dist, depth, fluo, all_of(taxa)) %>% 
      gather(all_of(taxa), key = taxon, value = conc)
    
    # Plot 
    p <- d %>% 
      ggplot() +
      geom_raster(aes(x = dist, y = -depth, fill = conc)) +
      scale_fill_viridis_c(trans = "log1p", na.value = NA) +
      geom_contour(aes(x = dist, y = -depth, z = fluo), breaks = c(0.05, 0.1), size = 0.3, color = "white", alpha = 0.8) +
      facet_wrap(~taxon, ncol = 3) +
      labs(x = "Distance (km)", y = "Depth (m)", fill = "Concentration<br>(ind.m<sup>-3</sup> )") +
      theme_minimal() +
      theme(legend.title = element_markdown(), text = element_text(size = 14)) +
      ggtitle(my_transect_title)
    
    # Save
    ggsave(p, file=str_c(output_dir, "/", my_transect, ".pdf"), width=15, height=15)
  }
  # List of newly created plots 
  my_plots <- list.files(output_dir, full.names = TRUE)
  
  # Combine all pages in one pdf
  pdf_combine(input = my_plots, output = str_c(output_dir, "/all_transects.pdf"))
  
  # Delete single page plots
  file.remove(my_plots)
  
  # Move multi transect file out of temporary dir
  file.rename(from = str_c(output_dir, "/all_transects.pdf"), to = output_file)
  
  # Delete temporary dir
  unlink(output_dir, recursive = TRUE)
  
}


## Split violin plot ----
#--------------------------------------------------------------------------#
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "width", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}



## Large mixed color palette ----
#--------------------------------------------------------------------------#
yo_cols <- c(
  "#a4db9a",
  "#aa3ad5",
  "#72e64b",
  "#704ee4",
  "#d7e93b",
  "#4228aa",
  "#73bd2c",
  "#dc50d0",
  "#93e373",
  "#532488",
  "#53e392",
  "#d43e99",
  "#47a242",
  "#b26dd5",
  "#a2b935",
  "#6472e1",
  "#e0be3b",
  "#3f195a",
  "#d2ea7e",
  "#92378c",
  "#7f9a3d",
  "#4f5099",
  "#e78e2e",
  "#5c8ed1",
  "#de5624",
  "#6ae1cb",
  "#db3438",
  "#6bc2d9",
  "#d83b6b",
  "#429d77",
  "#a93931",
  "#b7ddd1",
  "#391625",
  "#ddc572",
  "#2c264e",
  "#98892f",
  "#a385ca",
  "#3e6a28",
  "#e296d5",
  "#1a2825",
  "#dbd0a7",
  "#823260",
  "#889d68",
  "#cc6e8d",
  "#32543c",
  "#e27e67",
  "#4d8b9f",
  "#b87b30",
  "#305067",
  "#cb9b74",
  "#796f90",
  "#894e26",
  "#b9b9de",
  "#722627",
  "#7d9585",
  "#65404e",
  "#ddafb8",
  "#4b3b24",
  "#996f6a",
  "#6a6331"
)

