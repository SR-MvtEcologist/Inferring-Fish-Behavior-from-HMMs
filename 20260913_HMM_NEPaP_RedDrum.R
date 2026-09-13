##### --------------- #####
### 20260913: HMM For Assigning Behavior States to Triangulated Red Drum Positions at NEPaP ###
##### --------------- #####

## ------------------------- ##
# LOAD PACKAGES #
## ------------------------- ##
library(sf) # (Pebesma et al. 2026) 
library(sp) # (Pebesma et al. 2026) 
library(ggplot2) # (Wickham 2023)
library(momentuHMM) # (McClintock & Michelot 2025)
library(tidyverse) # (Wickham 2023)
library(lubridate) # (Spinu et al. 2026)
library(circular) # (Lund et al. 2025)


## ------------------------- ##
# LOAD DATA #
## ------------------------- ##
read.csv("...RedDrumHMMData.csv")


## ------------------------- ##
# PREP DATA #
## ------------------------- ##

# this step, unavoidably, takes a minute to run
data_hmm <- prepData(data = pos_hmm)


## ------------------------- ##
# Determine a sensible place in parameter space to start the model optimizer #
## ------------------------- ##

### Step Lengths ###
step <- data_hmm$step

summary(step)
mean(step, na.rm = TRUE)
sd(step, na.rm = TRUE)
quantile(step, c(.01, .05, .10, .25, .50, .75, .90, .95, .99),
         na.rm = TRUE)

hist(step,
     breaks = 100,
     probability = TRUE,
     main = "Step length",
     xlab = "Step length")
# note: a gamma distribution is "business as usual" for fish step length data

sum(data_hmm$step == 0, na.rm = TRUE)
mean(data_hmm$step == 0, na.rm = TRUE)
# no zeroes, no need for zero inflation parameter

# k-means on log step length to get initial state parameters
# Note Par order: mean1, mean2, (mean3), sd1, sd2, (sd3)
# 2-state
keep <- is.finite(step) & step >= 0
step_valid <- step[keep]
z <- log1p(step_valid)

km2 <- kmeans(z, centers = 2, nstart = 50)
g2 <- km2$cluster

mu2 <- tapply(step_valid, g2, mean)
sd2 <- tapply(step_valid, g2, sd)

mu2
sd2

ord2 <- order(mu2)
mu2 <- mu2[ord2]
sd2 <- sd2[ord2]

Par0.step.2 <- c(mu2, sd2)
Par0.step.2

# 3-state
km3 <- kmeans(z, centers = 3, nstart = 50)
g3 <- km3$cluster

mu3 <- tapply(step_valid, g3, mean)
sd3 <- tapply(step_valid, g3, sd)

mu3
sd3

ord3 <- order(mu3)
mu3 <- mu3[ord3]
sd3 <- sd3[ord3]

Par0.step.3 <- c(mu3, sd3)
Par0.step.3


### Turn Angles ###
range(data_hmm$angle, na.rm = TRUE) # should be -3.14 to 3.14, not -180 to 180
theta_valid <- data_hmm$angle[keep]
circ_stats <- function(x) {
  x <- x[is.finite(x)]
  
  z <- mean(exp(1i * x))
  R <- abs(z)
  mu <- Arg(z)
  
  c(mean = mu, R = R)
}

kappa_from_R <- function(R) {
  if (R < 0.53) {
    2 * R + R^3 + (5 * R^5) / 6
  } else if (R < 0.85) {
    -0.4 + 1.39 * R + 0.43 / (1 - R)
  } else {
    1 / (R^3 - 4 * R^2 + 3 * R)
  }
}

# 2-state
stats_angle_2 <- lapply(1:2, function(i) circ_stats(theta_valid[g2 == i]))
angle_mu_2 <- sapply(stats_angle_2, function(x) x["mean"])
angle_R_2  <- sapply(stats_angle_2, function(x) x["R"])

angle_mu_2 <- angle_mu_2[ord2]
angle_R_2  <- angle_R_2[ord2]

angle_kappa_2 <- sapply(angle_R_2, kappa_from_R)

Par0.angle.vm.2 <- c(angle_mu_2, angle_kappa_2)
Par0.angle.wrpcauchy.2 <- c(angle_mu_2, angle_R_2)

# 3-state
stats_angle_3 <- lapply(1:3, function(i) circ_stats(theta_valid[g3 == i]))
angle_mu_3 <- sapply(stats_angle_3, function(x) x["mean"])
angle_R_3  <- sapply(stats_angle_3, function(x) x["R"])

angle_mu_3 <- angle_mu_3[ord3]
angle_R_3  <- angle_R_3[ord3]

angle_kappa_3 <- sapply(angle_R_3, kappa_from_R)

Par0.angle.vm.3 <- c(angle_mu_3, angle_kappa_3)
Par0.angle.wrpcauchy.3 <- c(angle_mu_3, angle_R_3)


## ------------------------- ##
# 2-state with von Mises #
## ------------------------- ##

HMM2_vm <- fitHMM(
  data = data_hmm,
  nbStates = 2,
  dist = list(step = "gamma", angle = "vm"),
  Par0 = list(
    step = Par0.step.2,
    angle = Par0.angle.vm.2
  ),
  estAngleMean = list(angle = TRUE),
  retryFits = 10
)

HMM2_vm

states_m2_vm <- viterbi(HMM2_vm)
data_hmm$state_m2_vm <- states_m2_vm

## ------------------------- ##
# 2-state with wrapped Cauchy #
## ------------------------- ##

HMM2_wc <- fitHMM(
  data = data_hmm,
  nbStates = 2,
  dist = list(step = "gamma", angle = "wrpcauchy"),
  Par0 = list(
    step = Par0.step.2,
    angle = Par0.angle.wrpcauchy.2
  ),
  estAngleMean = list(angle = TRUE),
  retryFits = 10
)

HMM2_wc

states_m2_wc <- viterbi(HMM2_wc)
data_hmm$state_m2_wc <- states_m2_wc

## ------------------------- ##
# 3-state with von Mises #
## ------------------------- ##

HMM3_vm <- fitHMM(
  data = data_hmm,
  nbStates = 3,
  dist = list(step = "gamma", angle = "vm"),
  Par0 = list(
    step = Par0.step.3,
    angle = Par0.angle.vm.3
  ),
  estAngleMean = list(angle = TRUE),
  retryFits = 10
)

HMM3_vm

states_m3_vm <- viterbi(HMM3_vm)
data_hmm$state_m3_vm <- states_m3_vm

## ------------------------- ##
# 3-state with wrapped Cauchy #
## ------------------------- ##

HMM3_wc <- fitHMM(
  data = data_hmm,
  nbStates = 3,
  dist = list(step = "gamma", angle = "wrpcauchy"),
  Par0 = list(
    step = Par0.step.3,
    angle = Par0.angle.wrpcauchy.3
  ),
  estAngleMean = list(angle = TRUE),
  retryFits = 10
)

HMM3_wc

states_m3_wc <- viterbi(HMM3_wc)
data_hmm$state_m3_wc <- states_m3_wc

## ------------------------- ##
# model visualizations #
## ------------------------- ##

# distribution of step lengths by state
par(mfrow = c(1, 3))

for(i in 1:3) {
  hist(
    data_hmm$step[data_hmm$state_m3_wc == i],
    breaks = 50,
    main = paste("State", i),
    xlab = "Step length",
    probability = TRUE
  )
}

# distribution of angles by state
par(mfrow = c(1, 3))

for(i in 1:3) {
  hist(
    data_hmm$angle[data_hmm$state_m3_wc == i],
    breaks = 50,
    main = paste("State", i),
    xlab = "Angle",
    probability = TRUE
  )
}
# good variation among behaviors, fat tails on each = wrapped Cauchy suites the data well

# do not let plot() go ham and print every single color coded track
# print the first two diagnostic plots + 20 random plots to check
# OR, run this code last, because you will have to restart R before making any additional plots
# signed, someone whose R session just crashed
plot(HMM3_wc, plotTracks = FALSE)
 