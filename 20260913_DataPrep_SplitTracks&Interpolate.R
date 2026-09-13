# ---------- #
### 202609013: Hidden Markov model prep code
### split triangulated fish position data into tracks of desired duration
### interpolate series of positions into temporally regular tracks
# ---------- #

# ---------- #
### LOAD PACKAGES
# ---------- #
library(sf) # (Pebesma et al. 2026) 
library(sp) # (Pebesna et al. 2026)
library(tidyverse) # (Wickham 2023)
library(lubridate) # (Spinu et al. 2026)
library(ggplot2) # (Wickham et al. 2026)
library(crawl) # (Johnson et al. 2026)
library(dplyr) # (Wickham et al. 2026)

# ---------- # 
### LOAD DATA
# ---------- #
dat <- read.csv("...habitats_SO368b.28093_NEPaP.csv")


# ---------- #
### SPLIT TRACKS
# ---------- #

# Michelot's split tracks function calls for specific columns / headers
dat$ID <- dat$FullId
dat$time <- dat$DateTime
dat$lon <- dat$Longitude
dat$lat <- dat$Latitude

dat <- dat[, c("ID", "time", "lon", "lat")]

# load split tracks function
source("...splittracks_Michelot22.R")
# split tracks
data_split <- split_at_gap(data = dat, max_gap = 0.05*60, shortest_track = 0.25*60)

# view split tracks
ggplot(data_split, aes(lon, lat, col = ID)) + 
  geom_point(size = 0.5) + geom_path() +
  coord_equal()

# number of tracks
tracks1 <- length(unique(data_split$ID))


# ---------- #
### INTERPOLATE TRACKS
# ---------- #

# --------------------------------------------------
# 1 UNIT CONVERSIONS
# --------------------------------------------------

# ensure timestamp is R-compatible
data_split$time <- as.POSIXct(
  data_split$time,
  format = "%Y/%m/%d %H:%M:%S",
  tz = "UTC"
)

# convert lat and long into UTM
dat_sf <- st_as_sf(
  data_split,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

dat_utm <- st_transform(
  dat_sf,
  crs = 32616
)

xy <- st_coordinates(dat_utm)

dat_utm$x <- xy[, 1]
dat_utm$y <- xy[, 2]

# --------------------------------------------------
# 2 FIT ONE TRACK FUNCTION
# --------------------------------------------------

fit_one_track <- function(d, theta_values, attempts = 20) {
  
  for (i in seq_along(theta_values)) {
    
    theta <- theta_values[[i]]
    
    message(
      "    Attempt ", i,
      " | theta = ", paste(theta, collapse = ", ")
    )
    
    fit <- tryCatch(
      crwMLE(
        data = d,
        mov.model = ~ 1,
        err.model = NULL,
        coord = c("x", "y"),
        Time.name = "time",
        time.scale = "seconds",
        theta = theta,
        fixPar = c(NA, NA),
        attempts = attempts
      ),
      error = function(e) e
    )
    
    # If fitting succeeded, return the fit AND
    # information about which attempt worked
    if (!inherits(fit, "error")) {
      
      return(
        list(
          fit = fit,
          attempt = i,
          theta = theta
        )
      )
    }
  }
  
  # None of the starting values worked
  return(NULL)
}

# --------------------------------------------------
# 3 STARTING VALUES
# --------------------------------------------------

theta_values <- list(
  c(4, 0),    # original
  c(2, 0),
  c(6, 0),
  c(4, 1),
  c(4, -1),
  c(2, 1),
  c(6, -1)
)

# --------------------------------------------------
# 4 OBJECTS TO STORE RESULTS
# --------------------------------------------------

fits <- list()

failed <- character(0)

fit_log <- data.frame(
  ID = character(),
  attempt = integer(),
  theta1 = numeric(),
  theta2 = numeric(),
  stringsAsFactors = FALSE
)

# --------------------------------------------------
# 5 FIRST PASS: fit every track using theta = c(4, 0)
# --------------------------------------------------

for (id in unique(dat_utm$ID)) {
  
  message("\nFitting ID: ", id)
  
  d <- dat_utm %>%
    filter(ID == id) %>%
    arrange(time) %>%
    filter(
      !is.na(time),
      !is.na(x),
      !is.na(y)
    )
  
  result <- fit_one_track(
    d = d,
    theta_values = list(c(4, 0)),
    attempts = 20
  )
  
  if (is.null(result)) {
    
    failed <- c(failed, id)
    
    message("  FAILED")
    
  } else {
    
    fits[[id]] <- result$fit
    
    fit_log <- rbind(
      fit_log,
      data.frame(
        ID = id,
        attempt = result$attempt,
        theta1 = result$theta[1],
        theta2 = result$theta[2]
      )
    )
    
    message("  SUCCESS")
  }
}

failed

# --------------------------------------------------
# 6 SECOND PASS: retry ONLY failed tracks
# --------------------------------------------------

if (length(failed) > 0) {
  
  message("\n\nRETRYING FAILED TRACKS")
  
  for (id in failed) {
    
    message("\nRetrying ID: ", id)
    
    d <- dat_utm %>%
      filter(ID == id) %>%
      arrange(time) %>%
      filter(
        !is.na(time),
        !is.na(x),
        !is.na(y)
      )
    
    # Don't try c(4,0) again.
    # Start with the alternative theta values.
    retry_theta <- theta_values[-1]
    
    result <- fit_one_track(
      d = d,
      theta_values = retry_theta,
      attempts = 20
    )
    
    if (is.null(result)) {
      
      message("  STILL FAILED")
      
    } else {
      
      # Add successful retry to the existing fits
      fits[[id]] <- result$fit
      
      fit_log <- rbind(
        fit_log,
        data.frame(
          ID = id,
          attempt = result$attempt + 1,
          theta1 = result$theta[1],
          theta2 = result$theta[2]
        )
      )
      
      message(
        "  SUCCESS ON RETRY ",
        result$attempt + 1
      )
    }
  }
}

# --------------------------------------------------
# 7 Identify tracks that still failed
# --------------------------------------------------

failed_final <- setdiff(
  unique(dat_utm$ID),
  names(fits)
)

failed_final


# --------------------------------------------------
# 8 Look at fitting history
# --------------------------------------------------

fit_log

# --------------------------------------------------
# 9 CREATE TRACK PREDICTIONS 
# --------------------------------------------------
predictions <- lapply(names(fits), function(id) {
  
  pred <- crwPredict(
    fits[[id]],
    predTime = "180 secs"
  )
  
  # Remove crawl's special class BEFORE using dplyr
  class(pred) <- "data.frame"
  
  # Keep only predicted locations
  pred <- pred[pred$locType == "p", , drop = FALSE]
  
  # Make a clean data frame containing the smoothed locations
  pred <- data.frame(
    ID = id,
    ID_old = unique(
      dat_utm$ID_old[dat_utm$ID == id]
    ),
    time = pred$time,
    x = pred$mu.x,
    y = pred$mu.y
  )
  
  pred
})

names(predictions) <- names(fits)

smooth_tracks <- dplyr::bind_rows(predictions)

head(smooth_tracks)

# -------------------------------------------------
# 10 What failed and why?
# -------------------------------------------------
failed_final

# -------------------------------------------------
# 11 Add long and lat back for mapping
# -------------------------------------------------
smooth_sf <- st_as_sf(
  smooth_tracks,
  coords = c("x", "y"),
  crs = 32616
)

smooth_sf <- st_transform(
  smooth_sf,
  crs = 4326
)

coords <- st_coordinates(smooth_sf)

smooth_tracks$lon <- coords[, "X"]
smooth_tracks$lat <- coords[, "Y"]

head(smooth_tracks)