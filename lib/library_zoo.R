#--------------------------------------------------------------------------#
# Project: VISIIS
# Script purpose: Functions relavive to zooplankton data
# Date: 27/11/2020
# Author: Thelma Panaiotis
#--------------------------------------------------------------------------#


## Function to read TSV tables with object prediction ----
#--------------------------------------------------------------------------#
read_preds <- function(file, ignore_classes) {
  #' Read a TSV file of ISIIS particles predictions
  #'
  #' Returns the content of predicted and living organisms after a sample_id clean
  #' @param file TSV file with predictions to read
  #' @param ignore_classes list of classes to ignore (eg "detritus")
  
  library(lubridate)
  
  # Read TSV file
  pred <- read_tsv(file, skip=2, 
                   col_names = read_tsv(file) %>% colnames(),
                   col_types = cols(object_date = col_character(), object_time = col_character(), sample_id = col_character(), acq_id = col_character()))

  # Clean data
  pred <-  pred %>% 
    select(sample_id, acq_id, date=object_date, time=object_time, lat=object_lat, lon=object_lon,
           depth=object_depth_min, taxon=object_annotation_category, score=object_cnn_score) %>%  # keep relevant columns
    mutate(datetime = as_datetime(str_c(date, time), format = "%Y%m%d%H%M%S", tz = "Europe/Paris")) %>% # compute datetime with local time zone
    select(-c(date, time)) %>% 
    filter(!(taxon %in% ignore_classes)) #%>% # ignore non relevant taxa
    # and correct sample_id
    #mutate(
    #  transect = str_split_fixed(sample_id, '_', 2)[,1], # transect type
    #  yo = gsub("^.*_", "", sample_id), # extract yo
    #  yo = ifelse((yo == "prelag1") | (yo == "yonan"), NA, yo), # NA for yo of pre lag 1 and yonan
    #  yo = as.numeric(str_sub(yo, 3)), # convert yo number to numeric
    #) %>% 
    #fill(yo) %>% # fill missing value from top to down. This is ok because data is arranged following datetime. 
    #mutate(
    #  sample_id = paste(transect, str_pad(yo, 2, pad = "0"), sep = "_yo"), # is transect is not "prelag1", compute sample name from transect and yo, else keep "prelag1"
    #  cast_type = ifelse(yo %% 2 == 0, 'up', 'down'),
    #  yo = factor(yo) # convert yo to factor
    #) 
  
  return (pred)
}


## Interpolate (coarse) multiple variables across multiple transects ----
#--------------------------------------------------------------------------#

coarse_interp = function(df, vars, step_x, step_y, parallel = TRUE, n_cores = 12) {
  #' Interpolate variables across transects.
  #' 
  #' Compute the linear interpolation of one or multiple variables on one or multiple transects.
  #' with onto a specified output grid. 
  #' @param df Dataframe containing ISIIS data to interpolate, must contain a "transect" column
  #' @param vars list of variables to interpolate
  #' @param step_x dimension of output grid in x direction (dist)
  #' @param step_y dimension of output grid in y direction (depth)
  #' @param parallel whether to perform parallel computation
  #' @param n_cores number of cores to use for parallel computation
  
  library(akima)
  library(broom)
  library(tidyverse)
  library(parallel)
  
  # Parallel interpolation
  if (parallel){
    intl <- mclapply(unique(df$transect), function(id) {
      
      # Initiate empty tibble for this transect
      transect_int <- tibble()
      
      # filter data for transect
      df_t <- df %>% filter(transect == id)
      
      # Generate output grid for transect
      xo <- seq(floor(min(df_t$dist)),  ceiling(max(df_t$dist)), by = step_x)
      yo <- seq(floor(min(df_t$depth)), ceiling(max(df_t$depth)), by = step_y)
      
      # Loop over variables to interpolate
      for (variable in vars){
        # Compute variable interpolation
        var_int <- coarse_interp_1v(df_t, variable, xo=xo, yo=yo)
        
        # Join to table with other variables
        if (length(transect_int) == 0) { # If first interpolated variable on this transect
          transect_int <- var_int # Replace transect table by newly computed interpolation
        } else { # Else perform a left join with previously interpolated variables
          transect_int <- left_join(transect_int, var_int, by = c("transect", "dist", "depth"))
        } 
      }
      
      # Reorder columns
      transect_int <- transect_int %>% select(transect, dist, depth, everything())
      
      return(transect_int)
      
    }, mc.cores=n_cores) 
    # this returns a list, recombine it into a tibble
    df_int <- do.call(bind_rows, intl)
    
  } else{ # non parallel interpolation
    
    # Initiate empty tibble for interpolation results
    df_int <- tibble()
    
    # Loop over transects
    for(id in unique(df$transect)) {
      
      # Initiate empty tibble for this transect
      transect_int <- tibble()
      
      # filter data for transect
      df_t <- df %>% filter(transect == id)
      
      # Generate output grid for transect
      xo <- seq(floor(min(df_t$dist)),  ceiling(max(df_t$dist)), by = step_x)
      yo <- seq(floor(min(df_t$depth)), ceiling(max(df_t$depth)), by = step_y)
      
      # Loop over variables to interpolate
      for (variable in vars){
        # Compute variable interpolation
        var_int <- coarse_interp_1v(df_t, variable, xo=xo, yo=yo)
        
        # Join to table with other variables
        if (length(transect_int) == 0) { # If first interpolated variable on this transect
          transect_int <- var_int # Replace transect table by newly computed interpolation
        } else { # Else perform a left join with previously insterpolated variables
          transect_int <- left_join(transect_int, var_int)
        } 
      }
      
      # Append transect rows to other transects
      df_int <- bind_rows(df_int, transect_int)
    }
  }
  
  # Reorder columns
  df_int <- df_int %>% select(transect, dist, depth, everything())
  
  return(df_int)
}


## Interpolate (coarse) a variable across a transect ----
#--------------------------------------------------------------------------#
coarse_interp_1v = function(df, variable, xo, yo) {
  #' Interpolate a variable across a transect.
  #' 
  #' Compute the linear interpolation of a variable on a given transect 
  #' with onto a specified output grid. 
  #' @param df Dataframe containing ISIIS data to interpolate
  #' @param variable name of variable to interpolate
  #' @param xo vector of x-coordinate of output grid (dist)
  #' @param yo vector of y-coordinate of output grid (depth)
  
  library(akima)
  library(broom)
  library(tidyverse)

  
  # Get transect name
  transect <- df %>% pull(transect) %>% unique()
  
  # Keep relevant columns and drop NA values
  df <- df %>% 
    select(dist, depth, all_of(variable)) %>% 
    drop_na(all_of(variable))
  
  # Perform the interpolation 
  df_int <- akima::interp(
    x=df$dist, 
    y=df$depth, 
    z=df[variable] %>% pull(), 
    linear=TRUE, 
    duplicate = "mean",
    xo = xo,
    yo = yo,
  ) %>% 
    tidy() %>% # convert list to dataframe
    rename(dist = x, depth = y, value = z) %>% # rename columns
    # change name of column containing newly interpolated variable with mutate and spread
    mutate(
      variable = variable,
      transect = transect,
    ) %>% 
    spread(variable, value)
  
  return(df_int)
}


## Interpolate datetime from dist ----
#--------------------------------------------------------------------------#
# Interpolate datetime from dist in 1D
interp_datetime = function(df_bin, df_int) {
  #' Interpolate datetime from dist on multiple transects.
  #' 
  #' Compute the linear interpolation (with extrapolation) of datetime fram dist.
  #' @param df_bin Dataframe of binned data with datetime and dist
  #' @param df_int Dataframe of interpolated data on which to interpolate datetime
  
  library(castr)
  library(tidyverse)
  
  # initiate output tibble
  df_out <- tibble()
  
  # loop over transects
  for (tr in sort(unique(df_bin$transect))){
    # filter data for this transect
    d_bin <- df_bin %>% filter(transect == tr)
    d_int <- df_int %>% filter(transect == tr)
    
    # interpolate datetime for this transect
    d_int <- d_int %>% 
      mutate(
        datetime = interpolate(x = d_bin$dist, y = d_bin$datetime, xout = d_int$dist),
        datetime = as_datetime(datetime, tz = "Europe/Paris"),
        .after = depth
      )
    
    df_out <- bind_rows(df_out, d_int)
  }
  
  return(df_out)
}


## Delete interpolated values too far from original points ----
#--------------------------------------------------------------------------#
clean_zoo_interp = function(df_int, df_bin, vars, nn_x=1, nn_y=2, max_dist=2){
  #' Delete interpolated plankton data too far from original points.
  #' 
  #' For each interpolated point (x, y), compute distance to closest original point among 
  #' close points in the range (x ± nn_x, y ± nn_y). If distance is NA (no point in the 
  #' explored range) or above max_dist, set interpolated value associated to this point to NA.
  #' @param df_int dataframe of interpolated data
  #' @param df_bin dataframe of non interpolated data
  #' @param vars list of variables to perform cleaning on
  #' @param nn_x distance to which look for nearest neighbour in x (dist, km), default is 1
  #' @param nn_y distance to which look for nearest neighbour in y (depth, m), default is 2
  #' @param max_dist distance threshold above which interpolated data is deleted
  
  
  # Container column for distance to closest point
  df_int$dist_ref <-  NA
  
  # Loop over interpolated points
  for(i in 1:nrow(df_int)){
    x1 <-  df_int$dist[i]
    y1 <-  df_int$depth[i]
    
    # Keep only close points from non interpolated data
    df_ref <- df_bin %>% 
      filter(between(dist, x1 - nn_x, x1 + nn_x) & between(depth, y1 - nn_y, y1 + nn_y))
    
    if (nrow(df_ref) > 0){ # compute distance to all points in df_ref
      
      # loop over points of df_ref
      for(j in 1:nrow(df_ref)){
        x2 <- df_ref$dist[j]
        y2 <- df_ref$depth[j]
        current_dist <- euc_dist(x1, y1, x2, y2) # compute distance
        
        if(j == 1){ # if this is the first point in df_ref, set minimal distance to current distance
          min_dist <- current_dist
        }
        
        if(current_dist < min_dist){ # if current distance is lower than minimal distance, update min dist
          min_dist <- current_dist
        }
      }
      
      # fill values in the containers
      df_int$dist_ref[i] <-  min_dist
    }
  }
  
  # Delete interpolated plankton conc values too far from original points
  df_int <- df_int %>% 
    # keep only rows where distance was computed and is lowen than max_dist
    mutate(to_keep = ifelse(dist_ref > max_dist | is.na(dist_ref), FALSE, TRUE)) %>% 
    # replace interpolated plankton values by NAs when required
    mutate(across(all_of(vars), ~ ifelse(to_keep, ., NA))) %>% 
    select(-c(dist_ref, to_keep))
  
  return(df_int)
}


## Interpolate (coarse) multiple variables across multiple transects ----
#--------------------------------------------------------------------------#
fine_interp = function(df, vars, step_x, step_y, theta, parallel = TRUE, n_cores = 12) {
  #' Interpolate variables across transects.
  #' 
  #' Compute the linear interpolation of one or multiple variables on one or multiple transects.
  #' with onto a specified output grid. 
  #' @param df Dataframe containing ISIIS data to interpolate, must contain a "transect" column
  #' @param vars list of variables to interpolate
  #' @param step_x dimension of output grid in x direction (dist)
  #' @param step_y dimension of output grid in y direction (depth)
  #' @param parallel whether to perform parallel computation
  #' @param n_cores number of cores to use for parallel computation
  
  library(akima)
  library(broom)
  library(tidyverse)
  library(parallel)
  
  # Parallel interpolation
  if (parallel){
    intl <- mclapply(unique(df$transect), function(id) {
      
      # Initiate empty tibble for this transect
      transect_int <- tibble()
      
      # filter data for transect
      df_t <- df %>% filter(transect == id)
      
      # Loop over variables to interpolate
      for (variable in vars){
        # Compute variable interpolation
        var_int <- fine_interp_1v(df_t, variable, plankton = taxa, step_x_fine = step_x, step_y_fine = step_y, theta = theta)
        
        # Join to table with other variables
        if (length(transect_int) == 0) { # If first interpolated variable on this transect
          transect_int <- var_int # Replace transect table by newly computed interpolation
        } else { # Else perform a left join with previously interpolated variables
          transect_int <- left_join(transect_int, var_int, by = c("transect", "dist", "depth"))
        } 
      }
      
      # Reorder columns
      transect_int <- transect_int %>% select(transect, dist, depth, everything())
      
      return(transect_int)
      
    }, mc.cores=n_cores) 
    # this returns a list, recombine it into a tibble
    df_int <- do.call(bind_rows, intl)
    
  } else { # non parallel interpolation
    
    # Initiate empty tibble for interpolation results
    df_int <- tibble()
    
    # Loop over transects
    for(id in unique(df$transect)) {
      
      # Initiate empty tibble for this transect
      transect_int <- tibble()
      
      # filter data for transect
      df_t <- df %>% filter(transect == id)
      
      # Loop over variables to interpolate
      for (variable in vars){
        # Compute variable interpolation
        var_int <- fine_interp_1v(df_t, variable, plankton = taxa, step_x_fine = step_x, step_y_fine = step_y, theta = theta)
        
        # Join to table with other variables
        if (length(transect_int) == 0) { # If first interpolated variable on this transect
          transect_int <- var_int # Replace transect table by newly computed interpolation
        } else { # Else perform a left join with previously insterpolated variables
          transect_int <- left_join(transect_int, var_int, by = c("dist", "depth", "transect"))
        } 
      }
      
      # Append transect rows to other transects
      df_int <- bind_rows(df_int, transect_int)
    }
  }
  
  # Reorder columns
  df_int <- df_int %>% select(transect, dist, depth, everything())
  
  return(df_int)
}


## Interpolate (fine) a variable across a transect ----
#--------------------------------------------------------------------------#
fine_interp_1v = function(df_int, variable, plankton, step_x_fine, step_y_fine, theta = 0.5) {
  #' Interpolate a variable across a transect.
  #'
  #' Compute the linear interpolation of a variable on a given transect
  #' with onto a specified output grid.
  #' @param df_int dataframe of coarse interpolated data
  #' @param variable name of variale to interpolate 
  #' @param plankton list of plankton taxa (to ignore data around the thermocline)
  #' @param step_x_fine dimension of fine output grid in x direction (dist)
  #' @param step_y_fine dimension of fine output grid in y direction (depth)
  #' @param theta bandwidth or scale parameter
  
  
  library(akima)
  library(broom)
  library(tidyverse)
  
  # Get transect name
  transect <- df_int %>% pull(transect) %>% unique()
  
  ## Keep relevant columns and drop NA values
  df_int <- df_int %>%
    select(dist, depth, all_of(variable)) %>% 
    arrange(dist, depth)
  
  # Reformat data as akima::interp output
  x <- unique(df_int$dist)
  y <- unique(df_int$depth)
  z <- matrix(df_int[[variable]], nrow = length(x), byrow = TRUE)
  list_int <- list(x=x, y=y, z=z)
  
  # Perform fine interpolation
  df_int_fine <- list_int %>% 
    fields::interp.surface.grid(grid.list = list(
      x = seq(floor(min(x)), ceiling(max(x)), by = step_x_fine),
      y = seq(floor(min(y)), ceiling(max(y)), by = step_y_fine)
    )) %>%
    fields::image.smooth(theta = theta) %>%
    broom::tidy() %>%
    dplyr::rename(dist = x, depth = y, value = z) %>% # rename columns
    # change name of column containing newly interpolated variable with mutate and spread
    mutate(
      variable = variable,
      transect = transect
    ) %>%
    spread(variable, value)
  
  # Delete data in holes for plankton
  if (variable %in% plankton){
    
    df_int <- df_int %>% 
      rename(dist_round = dist, depth_round = depth, present = variable)
    
    # computed rounded distance and depth for fine interpolation grid
    df_int_fine <- df_int_fine %>% 
      mutate(
        dist_round = round(dist),
        depth_round = round(depth)
      ) %>% 
      # match with coarse grid
      left_join(df_int, by = c("dist_round", "depth_round")) %>% 
      mutate_at(variable, ~ ifelse(is.na(present), NA, .)) %>% 
      select(transect, dist, depth, all_of(variable))
    
  }
  
  return(df_int_fine)
}

## Compute euclidean distance between two points in 2D ----
#--------------------------------------------------------------------------#
euc_dist <-  function (x1, y1, x2, y2){
  dist = sqrt(((x2 - x1)^2) + ((y2-y1)^2))
  return(dist)
}


## Standardize a vector to mean = 0 and sd = 1 ----
#--------------------------------------------------------------------------#
scale2 <- function(x, na.rm = T){(x - mean(x, na.rm = na.rm)) / sd(x, na.rm)}


## Standardize a vector for values between 0 and 1 ----
#--------------------------------------------------------------------------#
scale3 <- function(x, na.rm = T){x / max(x, na.rm)}


## Renumber output clusters of mvpart ----
#--------------------------------------------------------------------------#
renumber_cl <- function(gr) { # Renumber output clusters of mvpart
  gr2 <- rep(0, length(gr))
  for (i in 1:length(gr)){
    gr2[i] = which(unique(gr) == gr[i])
  }
  return(gr2)  
} 

