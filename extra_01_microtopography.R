# Microtography

# Notes
# I am using 1 m resolution instead of the original 0.4 resolution, as the the 
# highest resolution seems to contain linear artifacts. Maybe related to 
# overlaps between flights?
# Also, resampling the DEM to 1 m increases the speed of computation.

library(terra)
library(magrittr)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/UDKIK_data/")

tmpfolder <- paste0(dir_dat, "/Temp/") %T>% dir.create()
terraOptions(
  # memfrac = 0.02,
  tempdir = tmpfolder
)

# Folder for output tiles

dir_microdem_tiles <- dir_dat %>%
  paste0(., "/microdem_tiles/") %T>%
  dir.create()

dir_tiles_nmins <- dir_microdem_tiles %>%
  paste0(., "/nmins/") %T>%
  dir.create()

dir_tiles_demmad <- dir_microdem_tiles %>%
  paste0(., "/demmad/") %T>%
  dir.create()

dir_tiles_aspsd <- dir_microdem_tiles %>%
  paste0(., "/aspsd/") %T>%
  dir.create()

dir_tiles_flowsd <- dir_microdem_tiles %>%
  paste0(., "/flowsd/") %T>%
  dir.create()

dir_tiles_slopeaspsd <- dir_microdem_tiles %>%
  paste0(., "/slopeaspsd/") %T>%
  dir.create()

dir_microdem_merged <- dir_dat %>%
  paste0(., "/microdem_merged/") %T>%
  dir.create()

# Focal matrix

myfocalmat <- matrix(c(1,1,1,1,NA,1,1,1,1), nrow = 3)

r_na <- rast(ncols=180, nrows=180, xmin=0)
mygaussmat <- focalMat(r_na, c(1,1), "Gauss")

mygaussmat2 <- focalMat(r_na, c(1,2), "Gauss")
mygaussmat5 <- focalMat(r_na, c(2.5,5), "Gauss")

# Local (pseudo) flow accumulation

localflow <- function(x) {
  x_center <- x[ceiling(length(x)/2)]
  out <- sum((x > x_center), na.rm = TRUE) - sum((x < x_center), na.rm = TRUE)
  return(out)
}

# C++ implementation of the above

library(Rcpp)

cppFunction( 
  "#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector localflowCpp(NumericVector x, std::size_t ni, std::size_t nw) {
  NumericVector out(ni);
  if (nw == 0) return out;            // nothing to do

  std::size_t start = 0;
  const std::size_t center_offset = (nw - 1) / 2; // window size is odd in focal

  for (std::size_t i = 0; i < ni; ++i) {
    std::size_t end = start + nw;

    double center = x[start + center_offset];
    int greater_count = 0;
    int smaller_count = 0;

    if (!NumericVector::is_na(center)) {
      for (std::size_t j = start; j < end; ++j) {
        double v = x[j];
        if (NumericVector::is_na(v)) continue;
        if (v > center)      ++greater_count;
        else if (v < center) ++smaller_count;
      }
    }
    // if center is NA, both counts stay 0 — same as R's na.rm=TRUE behavior
    out[i] = static_cast<double>(greater_count - smaller_count);
    start = end;
  }
  return out;
}"
)

# Function to calculate the angle based on the chord

cppFunction( 
  '#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double chordAngleCpp(NumericVector x) {
  if (x.size() < 6) {
    stop("Input vector must have at least 6 elements.");
  }

  double x1 = x[0], x2 = x[1], x3 = x[2],
         x4 = x[3], x5 = x[4], x6 = x[5];

  // Propagate NA like R would
  if (NumericVector::is_na(x1) || NumericVector::is_na(x2) ||
      NumericVector::is_na(x3) || NumericVector::is_na(x4) ||
      NumericVector::is_na(x5) || NumericVector::is_na(x6)) {
    return NA_REAL;
  }

  // crd <- sqrt((x[1]*x[3]-x[4]*x[6])^2 + (x[2]*x[3]-x[5]*x[6])^2)
  double t1 = x1 * x3 - x4 * x6;
  double t2 = x2 * x3 - x5 * x6;
  double crd = std::sqrt(t1 * t1 + t2 * t2);

  // out <- asin(crd/2) * 2  (in radians)
  double arg = crd / 2.0;
  // Clamp to avoid tiny numerical excursions beyond [-1, 1]
  if (arg > 1.0) arg = 1.0;
  if (arg < -1.0) arg = -1.0;

  return 2.0 * std::asin(arg);
}'
)

# Function to calculate differences in aspect

differ_aspect <- function(x) {
  library(circhelp)
  
  x <- demtile
  
  asp <- terrain(x, "aspect", unit = "radians")

  dem_focal <- focal(x, mygaussmat, fun = "mean", na.rm = TRUE)
  
  asp_focal <- dem_focal %>% terrain("aspect", unit = "radians")
  
  asp_diff2 <- angle_diff_rad(asp, asp_focal)
  
  asp_diff_agg <- asp_diff2 %>%
    focal(mygaussmat, fun = "sd", na.rm = TRUE) %>%
    aggregate(fact = 2, fun = "mean") %>%
    focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
    aggregate(fact = 2, fun = "mean") %>%
    focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
    aggregate(fact = 2, fun = "mean") %>%
    focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
    resample(r10, "bilinear") %>%
    round(3)
  
  return(asp_diff_agg)
}

# Path for DEMs

dem_zips <- list.files(
  "O:/AUIT_Geodata/Denmark/Digital_elevation_models/Lidar/DHM_2023_working",
  pattern = ".zip",
  full.names = TRUE
)

for (i in 1:length(dem_zips)) {
  
  # i <- 300
  
  unlink(paste0(tmpfolder, "*"))   # Delete everything in the temp folder
  
  unzip(
    dem_zips[i],
    exdir = tmpfolder
  )
  
  rasters <- list.files(
    tmpfolder,
    pattern = ".tif",
    full.names = TRUE
  )

  rasterlist_mins <- list()
  rasterlist_dem_mad <- list()
  rasterlist_aspsd <- list()
  rasterlist_flowsd <- list()
  rasterlist_slopeaspsd <- list()
  
  for (j in 1:length(rasters)) {
    # for (j in 1) {
    
    # j <- 1
    
    demtile0 <- rast(rasters[j])
    
    r2 <- rast(ext(demtile0), resolution = 1)
    
    r10 <- rast(ext(demtile0), resolution = 10)
    
    # demtile <- resample(demtile, r2, "average")
    
    demtile <- demtile0 %>%
      focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
      aggregate(fact = 2, fun = "mean")
    
    # Local depressions
    
    focalmin <- focal(demtile, myfocalmat, fun = "min", na.rm = TRUE)
    
    lessthanmin <- demtile < focalmin
    
    n_mins <- lessthanmin %>%
      focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
      aggregate(fact = 2, fun = "sum") %>%
      focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
      aggregate(fact = 2, fun = "sum") %>%
      focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
      aggregate(fact = 2, fun = "sum") %>%
      focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
      resample(r10, "bilinear") %>%
      multiply_by(100/6.4^2) %>%
      round(1)
    
    # n_mins <- aggregate(lessthanmin, fact = 10, fun = "sum")
    
    rasterlist_mins[[j]] <- n_mins
    
    # Roughness
    
    focalmeans <- focal(demtile, mygaussmat, fun = "mean", na.rm = TRUE)
    
    demdiffs <- abs(focalmeans - demtile)
    
    # dem_mad <- aggregate(demdiffs, fact = 10, fun = "median", na.rm = TRUE) %>%
    #   round(3)
    
    dem_mad <- demdiffs %>%
      focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
      aggregate(fact = 2, fun = "median") %>%
      focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
      aggregate(fact = 2, fun = "median") %>%
      focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
      aggregate(fact = 2, fun = "median") %>%
      focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
      resample(r10, "bilinear") %>%
      round(3)
    
    rasterlist_dem_mad[[j]] <- dem_mad
    
    # Differences in aspect
    
    tileasp_sd_agg <- differ_aspect(demtile)
    
    rasterlist_aspsd[[j]] <- tileasp_sd_agg
    
    # Tendency towards local flow accumulation
    
    tileflow <- focalCpp(
      demtile,
      w = 3,
      fun = localflowCpp,
      fillvalue = NA_real_
      )
    
    tileflow_sd <- tileflow %>%
      raise_to_power(2) %>%
      focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
      aggregate(fact = 2, fun = "mean") %>%
      focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
      aggregate(fact = 2, fun = "mean") %>%
      focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
      aggregate(fact = 2, fun = "mean") %>%
      focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
      resample(r10, "bilinear") %>%
      sqrt() %>%
      round(2)
    
    # Note: It is also possible to do channels/ridges and the mean value.
    # However, this will have to wait.
    
    # tileflow_sd <- aggregate(
    #   tileflow^2,
    #   fact = 10,
    #   fun = "sd",
    #   na.rm = TRUE
    # ) %>%
    #   sqrt() %>%
    #   round(2)
    
    rasterlist_flowsd[[j]] <- tileflow_sd
    
    # Standard deviation of slope aspect (giving less weight to flat areas)
    
    tileslope_sin <- terrain(demtile, "slope", unit="radians") %>% sin()

    tileasp <- terrain(demtile, "aspect", unit="radians")
    tileaspcos <- tileasp %>% cos()
    tileaspsin <- tileasp %>% sin()
    
    tileasp_focal <- terrain(focalmeans, "aspect", unit="radians")
    tileaspsin_focal <- tileasp_focal %>% sin()
    tileaspcos_focal <- tileasp_focal %>% cos()
    
    tileslope_smooth <- terrain(focalmeans, "slope", unit="radians") 
    tileslope_sin_smooth <- sin(tileslope_smooth)

    mystack3 <- c(
      tileaspcos, tileaspsin, tileslope_sin,
      tileaspcos_focal, tileaspsin_focal, tileslope_sin_smooth
    )

    angletest5 <- app(mystack3, chordAngleCpp)
    
    tileaslopeasp_sd_agg <- angletest5 %>%
      raise_to_power(2) %>%
      focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
      aggregate(fact = 2, fun = "mean") %>%
      focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
      aggregate(fact = 2, fun = "mean") %>%
      focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
      aggregate(fact = 2, fun = "mean") %>%
      focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
      resample(r10, "bilinear") %>%
      sqrt() %>%
      round(3)
    
    rasterlist_slopeaspsd[[j]] <- tileaslopeasp_sd_agg
    
    print(paste0("i = ", i, "; j = ", j))
  }
  
  # Aggregate 
  
  alltiles_nmins <- rasterlist_mins %>%
    sprc() %>%
    merge(
      filename = paste0(
        dir_tiles_nmins,
        "nmins_",
        sprintf("%03d", i),
        ".tif"
      ),
      names = "nmins",
      overwrite = TRUE
    )
  
  alltiles_demmad <- rasterlist_dem_mad %>%
    sprc() %>%
    merge(
      filename = paste0(
        dir_tiles_demmad,
        "demmad_",
        sprintf("%03d", i),
        ".tif"
      ),
      names = "nmins",
      overwrite = TRUE
    )
  
  alltiles_aspsd <- rasterlist_aspsd %>%
    sprc() %>%
    merge(
      filename = paste0(
        dir_tiles_aspsd,
        "aspsd_",
        sprintf("%03d", i),
        ".tif"
      ),
      names = "aspsd",
      overwrite = TRUE
    )
  
  alltiles_flowsd <- rasterlist_flowsd %>%
    sprc() %>%
    merge(
      filename = paste0(
        dir_tiles_flowsd,
        "flowsd_",
        sprintf("%03d", i),
        ".tif"
      ),
      names = "flowsd",
      overwrite = TRUE
    )
  
  alltiles_slopeaspsd <- rasterlist_slopeaspsd %>%
    sprc() %>%
    merge(
      filename = paste0(
        dir_tiles_slopeaspsd,
        "slopeaspsd_",
        sprintf("%03d", i),
        ".tif"
      ),
      names = "slopeaspsd",
      overwrite = TRUE
    )
}

# Mask and Write results to files

dem <- dir_dat %>%
  paste0(., "/Sampling_Input/dhm2015_terraen_10m.tif") %>%
  rast()

micro_nmins <- dir_tiles_nmins %>%
  list.files(".tif", full.names = TRUE) %>%
  sprc() %>%
  merge() %>%
  mask(
    mask = dem,
    filename = paste0(
      dir_microdem_merged,
      "micro_nmins.tif"
    ),
    names = "micro_nmins",
    overwrite = TRUE
  )

micro_demmad <- dir_tiles_demmad %>%
  list.files(".tif", full.names = TRUE) %>%
  sprc() %>%
  merge() %>%
  mask(
    mask = dem,
    filename = paste0(
      dir_microdem_merged,
      "micro_demmad.tif"
    ),
    names = "micro_demmad",
    overwrite = TRUE
  )
    

micro_aspsd <- dir_tiles_aspsd %>%
  list.files(".tif", full.names = TRUE) %>%
  sprc() %>%
  merge() %>%
  mask(
    mask = dem,
    filename = paste0(
      dir_microdem_merged,
      "micro_aspsd.tif"
    ),
    names = "micro_aspsd",
    overwrite = TRUE
  )


micro_flowsd <- dir_tiles_flowsd %>%
  list.files(".tif", full.names = TRUE) %>%
  sprc() %>%
  merge() %>%
  mask(
    mask = dem,
    filename = paste0(
      dir_microdem_merged,
      "micro_flowsd.tif"
    ),
    names = "micro_flowsd",
    overwrite = TRUE
  )

micro_slopeaspsd <- dir_tiles_slopeaspsd %>%
  list.files(".tif", full.names = TRUE) %>%
  sprc() %>%
  merge() %>%
  mask(
    mask = dem,
    filename = paste0(
      dir_microdem_merged,
      "micro_slopeaspsd.tif"
    ),
    names = "micro_slopeaspsd",
    overwrite = TRUE
  )

# Inspect results from first loop

# plot(demtile)
# plot(n_mins)
# plot(dem_mad)
# plot(tileasp_sd_agg)
# plot(tileflow_sd)

# Plot tiles for second loop

# alltiles_dem <- rasters %>%
#   sprc() %>%
#   merge() %>%
#   aggregate(fact = 25, fun = "mean")
# plot(alltiles_dem)
plot(alltiles_nmins)   # 1 decimal
plot(alltiles_demmad)   # 3 decimals
plot(alltiles_aspsd)  # 3 decimals
plot(alltiles_flowsd)  # 2 decimals
plot(alltiles_slopeaspsd) # 3 decimals

# Plot merged results

plot(micro_nmins)
plot(micro_demmad)
plot(micro_aspsd)
plot(micro_flowsd)
plot(micro_flowsd)

# Old stuff

# Test asp diff with slope
# plot(demtile)
# 
# tileslope <- terrain(demtile, "slope", unit="radians") 
# 
# tileslope_cotan <- tan(pi/2 - tileslope)
# 
# tileslope_sin <- sin(tileslope)
# 
# plot(tileslope_cotan)
# 
# tileasp <- terrain(demtile, "aspect", unit="radians")
# 
# tileaspcos <- tileasp %>% cos()
# 
# plot(tileaspcos)
# 
# tileaspsin <- tileasp %>% sin()
# 
# plot(tileaspsin)
# 
# tileasp_focal <- terrain(focalmeans, "aspect", unit="radians")
# 
# tileaspsin_focal <- tileasp_focal %>% sin()
# 
# plot(tileaspsin_focal)
# 
# tileaspcos_focal <- tileasp_focal %>% cos()
# 
# plot(tileaspcos_focal)
# 
# # tileslope_smooth <- focal(tileslope, mygaussmat, fun = "mean", na.rm = TRUE)
# 
# tileslope_smooth <- terrain(focalmeans, "slope", unit="radians") 
# 
# # tileslope_cotan_smooth <- focal(tileslope_cotan, mygaussmat, fun = "mean", na.rm = TRUE)
# 
# tileslope_cotan_smooth <- tan(pi/2 - tileslope_smooth)
# 
# tileslope_sin_smooth <- sin(tileslope_smooth)
# 
# plot(tileaspsin - tileaspsin_focal)
# plot(tileaspcos - tileaspcos_focal)
# 
# plot(tileaspsin*tileslope_sin_smooth)
# 
# plot(tileaspsin*tileslope_sin - tileaspsin_focal*tileslope_sin_smooth)
# 
# plot(tileslope_cotan_smooth)
# 
# mystack3 <- c(
#   tileaspcos, tileaspsin, tileslope_sin,
#   tileaspcos_focal, tileaspsin_focal, tileslope_sin_smooth
# )
# 
# angletest4 <- app(mystack3, function(x) {
#   a <- c(x[1], x[2], x[3])
#   b <- c(x[4], x[5], x[6])
#   dot.prod <- a%*%b 
#   norm.a <- norm(a,type = "2")
#   norm.b <- norm(b,type = "2")
#   theta <- acos(dot.prod / (norm.a * norm.b))
#   out <- as.numeric(theta)
#   return(out)
# })
# 
# plot(angletest4)
# 
# angletest5 <- app(mystack3, function(x) {
#   crd <- sqrt((x[1]*x[3]-x[4]*x[6])^2 + (x[2]*x[3]-x[5]*x[6])^2)
#   out <- asin(crd/2)*2
#   return(out)
# }
# )
# 
# plot(angletest5)
# 
# aggregate(angletest5, fact = 10, na.rm = TRUE) %>% plot()
# 
# angletest5 <- app(mystack3, function(x) {
#   z <- (x[3] + x[6])/2
#   a <- c(x[1], x[2], z)
#   b <- c(x[4], x[5], z)
#   dot.prod <- a%*%b 
#   norm.a <- norm(a,type = "2")
#   norm.b <- norm(b,type = "2")
#   theta <- acos(dot.prod / (norm.a * norm.b))
#   out <- as.numeric(theta)
#   return(out)
# })
# 
# angletest5 <- app(mystack3, chordAngleCpp)
# 
# plot(angletest5)
# 
# cppFunction( 
#   '#include <Rcpp.h>
# #include <cmath>
# using namespace Rcpp;
# 
# // [[Rcpp::export]]
# double chordAngleCpp(NumericVector x) {
#   if (x.size() < 6) {
#     stop("Input vector must have at least 6 elements.");
#   }
# 
#   double x1 = x[0], x2 = x[1], x3 = x[2],
#          x4 = x[3], x5 = x[4], x6 = x[5];
# 
#   // Propagate NA like R would
#   if (NumericVector::is_na(x1) || NumericVector::is_na(x2) ||
#       NumericVector::is_na(x3) || NumericVector::is_na(x4) ||
#       NumericVector::is_na(x5) || NumericVector::is_na(x6)) {
#     return NA_REAL;
#   }
# 
#   // crd <- sqrt((x[1]*x[3]-x[4]*x[6])^2 + (x[2]*x[3]-x[5]*x[6])^2)
#   double t1 = x1 * x3 - x4 * x6;
#   double t2 = x2 * x3 - x5 * x6;
#   double crd = std::sqrt(t1 * t1 + t2 * t2);
# 
#   // out <- asin(crd/2) * 2  (in radians)
#   double arg = crd / 2.0;
#   // Clamp to avoid tiny numerical excursions beyond [-1, 1]
#   if (arg > 1.0) arg = 1.0;
#   if (arg < -1.0) arg = -1.0;
# 
#   return 2.0 * std::asin(arg);
# }'
# )
# 
# angletest5 <- app(mystack3, chordAngleCpp)
# 
# plot(angletest5)
# 
# aggregate(angletest5, fact = 10, na.rm = TRUE) %>% plot()
# 
# angletest5 %>%
#   raise_to_power(2) %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r10, "bilinear") %>%
#   sqrt() %>%
#   round(3) %>% plot()
# 
# angle <- function(x,y){
#   dot.prod <- x%*%y 
#   norm.x <- norm(x,type="2")
#   norm.y <- norm(y,type="2")
#   theta <- acos(dot.prod / (norm.x * norm.y))
#   as.numeric(theta)
# }


# Tests of different approaches

# plot(terrain(demtile, "slope"))
# 
# terrain(demtile, "slope") %>%
#   aggregate( fact = 25, fun = "mean", na.rm = TRUE) %>%
#   plot()
# 
# plot(terrain(aggregate(demtile, fact = 25, fun = "mean", na.rm = TRUE), "slope"))
# 
# terrain(demtile, "TRIrmsd") %>%
#   aggregate(fact = 25, fun = "median", na.rm = TRUE) %>%
#   plot()
# 
# terrain(demtile, "slope") %>%
#   aggregate(fact = 25, fun = "median", na.rm = TRUE) %>%
#   plot()
# 
# 
# 
# tileasp_diff_agg <- differ_aspect(demtile)
# 
# plot(tileasp_diff_agg)
# 
# detrend_asp_diff <- differ_aspect(demdiffs)
# 
# plot(detrend_asp_diff)
# 
# # Aspect
# 
# tileslope <- terrain(demtile, "slope", unit = "radians") %>%
#   sin()
# 
# tileasp <- terrain(demtile, "aspect", unit = "radians")
# 
# plot(tileasp)
# 
# # sine
# 
# tileaspsin <- tileasp %>% sin()
# 
# plot(tileaspsin)
# 
# tileaspsin_focal <- focal(tileaspsin, mygaussmat, fun = "mean", na.rm = TRUE)
# 
# plot(tileaspsin_focal)
# 
# tileaspsin_diff <- tileaspsin - tileaspsin_focal
# 
# plot(tileaspsin_diff)
# 
# tileaspsin_agg <- aggregate(tileaspsin_diff^2, fact = 25, fun = "median", na.rm = TRUE) %>%
#   sqrt()
# 
# plot(tileaspcos_agg)
# 
# # cosine
# 
# tileaspcos <- tileasp %>% cos()
# 
# plot(tileaspcos)
# 
# tileaspcos_focal <- focal(tileaspcos, mygaussmat, fun = "mean", na.rm = TRUE)
# 
# plot(tileaspcos_focal)
# 
# tileaspcos_diff <- tileaspcos - tileaspcos_focal
# 
# plot(tileaspcos_diff)
# 
# tileaspcos_agg <- aggregate(tileaspcos_diff^2, fact = 25, fun = "median", na.rm = TRUE) %>%
#   sqrt()
# 
# plot(tileaspcos_agg)
# 
# # Combine aggregated sine and cosine
# 
# tileaspagg_diff <- atan2(tileaspsin_agg, tileaspcos_agg)
# 
# plot(tileaspagg_diff)
# 
# # Combine raw sine and cosine
# 
# library(circhelp)
# 
# tileasp2 <- atan2(tileaspsin, tileaspcos)
# 
# plot(tileasp2)
# 
# tileasp_focal <- atan2(tileaspsin_focal, tileaspcos_focal)
# 
# plot(tileasp_focal)
# 
# tileasp_diff2 <- angle_diff_rad(tileasp2, tileasp_focal) %>% abs()
# 
# plot(tileasp_diff2)
# 
# aggregate(
#   tileasp_diff2^2, fact = 25, fun = "mean", na.rm = TRUE) %>%
#   sqrt() %>%
#   plot()
# 
# tileasp_diff <- atan2(tileaspsin_diff, tileaspcos_diff)
# 
# plot(tileasp)
# 
# plot(tileasp_diff)
# 
# plot(abs(tileasp_diff))
# 
# plot(abs(tileasp_diff) > pi*0.5)
# 
# aggregate(
#   abs(tileasp_diff) > pi*0.5, fact = 25, fun = "sum", na.rm = TRUE) %>%
#   plot()
# 
# aggregate(
#   abs(tileasp_diff), fact = 25, fun = "mean", na.rm = TRUE) %>%
#   plot()
# 
# 
# tileasp_diff_agg <- aggregate(
#   tileasp_diff^2, fact = 25, fun = "mean", na.rm = TRUE) %>%
#   sqrt()
# 
# plot(tileasp_diff_agg)  # This one has promise
# 
# tileasp_diff_sd <- aggregate(
#   tileasp_diff^2, fact = 25, fun = "sd", na.rm = TRUE) %>%
#   sqrt()
# 
# plot(tileasp_diff_sd)
# 
# # Slope
# 
# plot(tileslope)
# 
# focalslope <- focal(tileslope, mygaussmat, fun = "mean", na.rm = TRUE)
# 
# slopediff <- tileslope - focalslope
# 
# plot(slopediff)
# 
# slopediff_agg <- aggregate(
#   slopediff^2, fact = 25, fun = "median", na.rm = TRUE) %>%
#   sqrt()
# 
# plot(slopediff_agg)
# 
# # Different focal filters
# 
# # Laplacian filter
# 
# dem_lap <- focal(demtile, matrix(c(0,1,0,1,-4,1,0,1,0), nrow=3), fun = "sum")
# 
# plot(dem_lap)
# 
# dem_lap_agg <- aggregate(
#   dem_lap, fact = 25, fun = "mean", na.rm = TRUE)
# 
# plot(dem_lap_agg)
# 
# # Sobel filter
# 
# fx <- matrix(c(-1,-2,-1,0,0,0,1,2,1), nrow=3)
# 
# fy <- matrix(c(1,0,-1,2,0,-2,1,0,-1), nrow=3)
# 
# dem_sob <- focal(demtile, fx, fun = "sum")
# 
# plot(dem_sob)
# 
# # Asp also using slope
# 
# plot(tileasp_diff*tileslope)
# 
# tileaspslope_diff_agg <- aggregate(tileasp_diff*tileslope, fact = 25, fun = "median", na.rm = TRUE)
# 
# plot(tileaspslope_diff_agg)
# 
# tileaspsin_slope <- tileasp %>% sin() %>% multiply_by(tileslope)
# 
# tileaspsin_slope_focal <- focal(tileaspsin_slope, mygaussmat, fun = "mean", na.rm = TRUE)
# 
# plot(tileaspsin_slope - tileaspsin_slope_focal)
# 
# aggregate((tileaspsin_slope - tileaspsin_slope_focal)^2, fact = 25, fun = "median", na.rm = TRUE) %>% sqrt() %>% plot()
# 
# terrain(demtile, "aspect", unit = "radians") %>%
#   cos() %>%
#   aggregate(fact = 25, fun = "mean", na.rm = TRUE) %>%
#   plot()
# 
# demtile %>%
#   aggregate(fact = 25, fun = "mean", na.rm = TRUE) %>%
#   terrain("aspect", unit = "radians") %>%
#   cos() %>%
#   plot()
# 
# # using smooth dem as baseline
# 
# smoothdem_aspsin <- focal(demtile, mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   terrain("aspect", unit = "radians") %>%
#   sin()
#   
# plot(smoothdem_aspsin)
# 
# plot(tileaspsin - smoothdem_aspsin)
# 
# aggregate((tileaspsin - smoothdem_aspsin)^2, fact = 25, fun = "median", na.rm = TRUE) %>% sqrt() %>% plot()
# 
# # Local flow accumulation
# 
# localflow <- function(x) {
#   xmat <- matrix(x, nrow = sqrt(length(x))) %>% rast()
#   flow <- terrain(xmat, "flowdir") %>% flowAccumulation()
#   out <- flow %>% values() %>% unlist() %>% magrittr::extract(5)
#   return(out)
# }
# 
# f <- system.file("ex/elev.tif", package="terra")
# r <- rast(f)
# 
# plot(r)
# 
# r_focalized <- focal(r, w = 3, localflow)
# 
# plot(r_focalized)
# 
# mynumbers <- c(1:9)
# 
# localflow(mynumbers)
# 
# localflow_R <- function(x) {
#   x_center <- x[ceiling(length(x)/2)]
#   out <- sum((x > x_center), na.rm = TRUE) - sum((x < x_center), na.rm = TRUE)
#   return(out)
# }
# 
# r_focalized <- focal(demtile, w = 3, localflow)
# 
# plot(r_focalized)
# 
# aggregate(r_focalized^2, fact = 25, fun = "mean", na.rm = TRUE) %>% sqrt() %>% plot()
# 
# aggregate(abs(r_focalized), fact = 25, fun = "mean", na.rm = TRUE) %>% plot()
# 
# aggregate(r_focalized, fact = 25, fun = "sd", na.rm = TRUE) %>% plot()
# 
# raw <- focalCpp(demtile, w=3, fun=localflowCpp, fillvalue=NA_real_)

# END