# Generate sampling zones

library(terra)
library(tidyr)
library(dplyr)
library(magrittr)
library(tibble)

# Load samplekmeans
# library(devtools)
# 
# install_github("anbm-dk/samplekmeans", force = TRUE)

library(samplekmeans)

# First try 2025-03-12

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/UDKIK_data/")

dir_sample_input <- dir_dat %>%
  paste0(., "/Sampling_Input/")
dir_sample_output <- dir_dat %>%
  paste0(., "/Sampling_zones_pts/")

tmpfolder <- paste0(dir_dat, "/Temp/") %T>% dir.create()
terraOptions(
  # memfrac = 0.02,
  tempdir = tmpfolder
  )

# Load shapefile with administrative area (may change)

peat_admin_area <- dir_sample_input %>%
  paste0(., "/p68_fields_exp_CVR_1ha.shp") %>%
  vect()

# Load new peat probability map (using 2025 combination)
# 
# peat_prob <- dir_sample_input %>%
#   paste0(., "/peat_probability_combined_sel.tif") %>%
#   rast()

# Reduce the file size

outname_prob <- "peat_probability_2025"
# 
# math(
#   peat_prob,
#   "round",
#   digits = 1,
#   filename = paste0(dir_sample_input, outname_prob, ".tif"),
#   overwrite = TRUE,
#   gdal = "TILED=YES",
#   datatype = "FLT4S",
#   NAflag = -9,
#   names = outname_prob
#   )

peat_prob_2025 <- paste0(dir_sample_input, outname_prob, ".tif") %>%
  rast()

rast("C:/Users/au542768/UDKIK/Kulstof2022_usikkerheder_2025/Layers_delivered/peat_probability_2025.tif")

# Smoothen probability raster

myfilter1 <- round(
  focalMat(peat_prob_2025, c(10, 20), "Gauss"),
  3
)

# prob_sum <- terra::focal(
#   peat_prob_2025,
#   w = myfilter1,
#   fun = "sum",
#   na.policy = "all",
#   na.rm = TRUE
# )
# 
# prob_ext_sum <-  terra::focal(
#   is.finite(peat_prob_2025),
#   w = myfilter1,
#   fun = "sum",
#   na.policy = "all",
#   na.rm = TRUE
# )
# 
# math(
#   prob_sum/prob_ext_sum,
#   "round",
#   digits = 1,
#   filename = paste0(dir_sample_input, outname_prob, "_smooth.tif"),
#   overwrite = TRUE,
#   gdal = "TILED=YES",
#   datatype = "FLT4S",
#   NAflag = -9,
#   names = outname_prob
#   )

# peat_prob_2025_smooth2 <- resample(
#   peat_prob_2025_smooth,
#   dem,
#   method = "near",
#   filename = paste0(dir_sample_input, outname_prob, "_smooth2.tif"),
# )

peat_prob_2025_smooth <- paste0(
  dir_sample_input, 
  outname_prob, 
  "_smooth.tif"
) %>%
  rast()

# Setup

dem <- paste0(
  dir_sample_input, 
  "/dhm2015_terraen_10m.tif"
) %>%
  rast()

# vdchn <- paste0(
#   dir_sample_input, 
#   "/vdtochn.tif"
# ) %>%
#   rast()

ext(peat_prob_2025)

# peat_prob_2025_res <- resample(
#   peat_prob_2025,
#   dem,
#   method = "near",
#   filename = paste0(dir_sample_input, outname_prob, "_resample.tif"),
# )
  
peat_prob_2025_res <- paste0(
  dir_sample_input,
  outname_prob, "_resample.tif"
  ) %>%
  rast()

sample_layers <- c(
  # peat_prob_2025_smooth,
  peat_prob_2025_res,
  dem
  # ,
  # vdchn
)

# General statistics for the sample layers

# sample_lyrs_pts <- sample_layers %>%
#   crop(
#     peat_admin_area,
#     mask = TRUE,
#     touches = FALSE 
#   ) %>%
#   spatSample(
#     size = 1000,  # Increase this number later
#     xy = TRUE,
#     na.rm = TRUE,
#     exhaustive = TRUE
#   )
# 
# iqr_sample_lyrs <- apply(
#   sample_lyrs_pts,
#   2,
#   function(x) {
#     quantile(x, 0.75, na.rm = TRUE) - quantile(x, 0.25, na.rm = TRUE)
#   }
# ) %>%
#   unlist()

# IQR loop
iqrs <- matrix(numeric(), nrow = nrow(peat_admin_area), ncol = 4)

for (i in 1:nrow(peat_admin_area)) {
  input_i <- terra::crop(
    sample_layers,
    peat_admin_area[i],
    mask = TRUE,
    touches = FALSE  # 2025-07-16
  )
  
  iqr_i <- input_i %>%
    as.data.frame(xy = TRUE) %>%
    apply(
      2,
      function(x) {
        quantile(x, 0.75, na.rm = TRUE) - quantile(x, 0.25, na.rm = TRUE)
      }
    ) %>%
    unlist()
  
  iqrs[i, ] <- iqr_i
}

iqrs_mean <- iqrs %>% apply(2, mean)
# [1] 189.0578042 179.0017571   8.5085100   0.5386088

# iqrs %>% apply(2, function(x) { weighted.mean(x, peat_admin_area$Shape_Area) } )

# Prepare sampling loop
  
n_samples_area <- peat_admin_area$Shape_Area %>%
  divide_by(10^4) %>%
  round(digits = 0)

set.seed(1)
ID_random <- sample.int(sum(n_samples_area), sum(n_samples_area))

start_ID_area <- c(0, cumsum(n_samples_area)) + 1

dir_zones_small <- dir_sample_output %>%
  paste0(., "/zones_small/") %T>%
  dir.create()

outpts <- list()
  
# Sampling loop
  
  
for (i in 1:nrow(peat_admin_area)) {
# for (i in 1:3) {
  
  input_i <- terra::crop(
    sample_layers,
    peat_admin_area[i],
    mask = TRUE,
    touches = FALSE  # 2025-07-16
  )
  
  weights_i <- is.finite(sum(input_i)) %>%
    focal(
      w = myfilter1,
      # w = 3,
      fun = "sum",
      na.policy = "all",
      na.rm = TRUE
      ) %>%
    mask(sum(input_i)) %>%
    raise_to_power(2)
  
  # Define input layers
  
  # threshold_i <- c(
  #   7,
  #   0.0687*sqrt(n_samples_area[i]) + 0.2747
  # )
  # 
  # iqr_i <- global(
  #   input_i,
  #   function(x) {
  #     quantile(x, 0.75, na.rm = TRUE) - quantile(x, 0.25, na.rm = TRUE)
  #   }
  # ) %>%
  #   unlist()
  # 
  # iqrs[i, ] <- iqr_i
  
  # if (sum(iqr_i > threshold_i) == 0) {
  #   xy_only_i <- TRUE
  # } else {
  #   xy_only_i <- FALSE
  #   if (sum(iqr_i > threshold_i) == 1) {
  #     input_i <- subset(input_i, iqr_i > threshold_i)
  #   }
  # }
  
  xy_only_i <- FALSE

  layers_scalers_i <-   unname(iqrs[i, ] / iqrs_mean)

  xy_weights_i <- mean(layers_scalers_i[1:2])
  layer_weights_i <- layers_scalers_i[3:4]
  
  # Run clustering
  
  samples_i <- sample_kmeans(
    input = input_i,
    weights = weights_i,
    clusters = n_samples_area[i],
    use_xy = TRUE,
    only_xy = xy_only_i,
    initializer = "kmeans++",
    pca = !xy_only_i,
    tol_pca = 0.05,
    xy_weight = xy_weights_i,
    layer_weights = layer_weights_i
  )
  
  IDs_i <- ID_random[(start_ID_area[i]):(start_ID_area[i + 1] - 1)]
  
  outpts_i <- samples_i$points %>%
    mutate(ID = IDs_i[ID])
  
  classify(
    samples_i$clusters,
    rcl = c(1:length(IDs_i), IDs_i) %>% matrix(ncol = 2),
    filename = paste0(dir_zones_small, "zones_", i, ".tif"),
    overwrite = TRUE,
    gdal = "TILED=YES",
    datatype = "INT4S",
    NAflag = -9,
    names = "zones"
  )
  
  outpts[[i]] <- outpts_i
}

outpts %<>% bind_rows()

outzones <- dir_zones_small %>%
  list.files(full.names = TRUE) %>%
  sprc() %>%
  merge()

writeRaster(
  outzones,
  filename = paste0(dir_sample_output, "/zones_all.tif"),
  overwrite = TRUE,
  gdal = "TILED=YES",
  datatype = "INT4S",
  NAflag = -9,
  names = outname_prob
)

outpts_vec <- vect(outpts, geom = c("x", "y"), crs = crs(outzones))

writeVector(
  outpts_vec,
  filename = paste0(dir_sample_output, "/samplepoints_all.shp"),
  overwrite = TRUE 
)

# Debug

plot(samples_i$clusters)
points(samples_i$points, pch = 21, bg = "white")

samples_debug <- sample_kmeans(
  input = input_i,
  weights = weights_i,
  clusters = clusters_try,
  use_xy = TRUE,
  only_xy = xy_only_i,
  initializer = "kmeans++",
  # ,
  # pca = TRUE,
  # tol_pca = 0.05
)

plot(samples_debug$clusters)
points(samples_debug$points)

# END

