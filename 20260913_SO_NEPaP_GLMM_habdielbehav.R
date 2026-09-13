##### --------------- #####
### 20260913: GLMM
### predictors: habitat, diel period
### response: fish behavior
### random effect: fish ID
##### --------------- #####

### ---------------- ###
# LOAD PACKAGES #
### ---------------- ###
library(lubridate) # (Spinu et al. 2026)
library(dplyr) # (Wickham et al. 2026)
library(mclogit) # (Elff 2025)
library(tidyr) # (Wickham et al. 2025)
library(ggplot2) # (Wickham 2023)
library(patchwork) # (Pedersen 2025)
library(cowplot) # (Wilke 2025)
library(grid) # ()
library(gtable) # (Wickham et al. 2024)


### ---------------- ###
# LOAD DATA #
### ---------------- ###
dat <- read.csv("...habitats_NEPaP_SO_HMM3wc_diel.csv")

# replace any blank habitats with open water
dat$habitat[dat$habitat == ""] <- "open"

# combine WADs into one habitat type
dat$habitat <- ifelse(grepl("^WAD", dat$habitat), "WAUs", dat$habitat)


### ---------------- ###
# BUILD MODEL #
### ---------------- ###

# make sure everything is being read in as factors
dat$state_m3_w <- factor(dat$state_m3_w)
dat$diel <- factor(dat$diel)
dat$habitat <- factor(dat$habitat)
dat$ID_old <- factor(dat$ID_old)
dat$ID <- factor(dat$ID)

# insufficient data to support a 4-level diel period term, so group dawn with day, and dusk with night
dat_mod <- dat
dat_mod$diel3 <- ifelse(
  dat_mod$diel %in% c("dawn", "day"),
  "day",
  "night"
)

dat_mod$diel3 <- factor(dat_mod$diel3)

# model
levels(dat_mod$diel3)
levels(dat_mod$habitat)
levels(dat_mod$state_m3_w)

dat_mod$diel3 <- relevel(dat_mod$diel3, ref = "day")
dat_mod$habitat <- relevel(dat_mod$habitat, ref = "open")
dat_mod$state_m3_w <- relevel(dat_mod$state_m3_w, ref = "1")

# only 97 em_veg positions, 12 receiver positions, and 1 WAD position
# result = multiple interaction categories with 0 positions
# result = complete model failure
# as those habitats clearly don't mean much to red drum, we remove them
dat_mod <- dat_mod %>%
  filter(!habitat %in% c("em_veg", "receiver", "WAUs"))

dat_mod <- droplevels(dat_mod)

dat_mod$diel3 <- relevel(dat_mod$diel3, ref = "day")
dat_mod$habitat <- relevel(dat_mod$habitat, ref = "open")
dat_mod$state_m3_w <- relevel(dat_mod$state_m3_w, ref = "1")

mod_daynight <- mblogit(
  state_m3_w ~ diel3 * habitat,
  random = ~1 | ID_old,
  data = dat_mod
)

### ---------------- ###
# MODEL CHECKS #
### ---------------- ###

# observed vs fitted probabilities
pred_fixed <- predict(
  mod_daynight,
  type = "response",
  conditional = FALSE
)

pred_random <- predict(
  mod_daynight,
  type = "response",
  conditional = TRUE
)

head(pred_fixed)
head(pred_random)

pred_fixed <- as.data.frame(pred_fixed)

dat_diag <- dat_mod

dat_diag$pred1 <- pred_fixed$`1`
dat_diag$pred2 <- pred_fixed$`2`
dat_diag$pred3 <- pred_fixed$`3`

check <- dat_diag %>%
  group_by(diel3, habitat) %>%
  summarise(
    n = n(),
    
    observed_1 = mean(state_m3_w == "1"),
    predicted_1 = mean(pred1),
    
    observed_2 = mean(state_m3_w == "2"),
    predicted_2 = mean(pred2),
    
    observed_3 = mean(state_m3_w == "3"),
    predicted_3 = mean(pred3),
    
    .groups = "drop"
  )

check

check_long <- check %>%
  pivot_longer(
    cols = c(
      observed_1, predicted_1,
      observed_2, predicted_2,
      observed_3, predicted_3
    ),
    names_to = c("type", "state"),
    names_pattern = "(observed|predicted)_(1|2|3)",
    values_to = "probability"
  ) %>%
  pivot_wider(
    names_from = type,
    values_from = probability
  )

ggplot(check_long,
       aes(x = predicted, y = observed, color = state)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2,
    color = "gray50"
  ) +
  geom_point(size = 3) +
  coord_equal() +
  scale_color_manual(
    values = c(
      "1" = "#1b9e77",
      "2" = "#d95f02",
      "3" = "#7570b3"
    )
  ) +
  labs(
    x = "Population-level predicted probability",
    y = "Observed proportion",
    color = "State"
  ) +
  theme_bw()

### ---------------- ###
# ARE ALL TERMS SIGNIFICANT #
### ---------------- ###

# mblogit is not supported by drop1, so we go old school

# additive model
mod_daynight_add <- mblogit(
  state_m3_w ~ diel3 + habitat,
  random = ~1 | ID_old,
  data = dat_mod
)

# comparison
anova(mod_daynight_add, mod_daynight, test = "Chisq")
AIC(mod_daynight_add, mod_daynight)
AIC(mod_daynight) - AIC(mod_daynight_add)

### ---------------- ###
# MODEL PREDICTIONS #
### ---------------- ###

# mclogit is not compatible with emmeans, so we have to do this by hand
# do I have any idea what most of the below code is doing? not really
# does it appear to be producing a valid emmeans-esque table? yes

mblogit_predictions <- function(model, data, conf.level = 0.95) {
  
  # Levels of predictors and response
  diel_levels <- levels(data$diel3)
  habitat_levels <- levels(data$habitat)
  state_levels <- levels(data$state_m3_w)
  
  # All combinations for which predictions are wanted
  newdat <- expand.grid(
    diel3 = diel_levels,
    habitat = habitat_levels,
    ID_old = data$ID_old[1]
  )
  
  newdat$diel3 <- factor(
    newdat$diel3,
    levels = levels(data$diel3)
  )
  
  newdat$habitat <- factor(
    newdat$habitat,
    levels = levels(data$habitat)
  )
  
  newdat$ID_old <- factor(
    newdat$ID_old,
    levels = levels(data$ID_old)
  )
  
  # Population-level predictions
  pred <- predict(
    model,
    newdata = newdat,
    type = "response",
    se.fit = TRUE,
    conditional = FALSE
  )
  
  # Predictions
  fit_long <- as.data.frame(pred$fit) %>%
    mutate(
      diel3 = newdat$diel3,
      habitat = newdat$habitat
    ) %>%
    pivot_longer(
      cols = all_of(state_levels),
      names_to = "state_m3_w",
      values_to = "estimate"
    )
  
  # Standard errors
  se_long <- as.data.frame(pred$se.fit) %>%
    mutate(
      diel3 = newdat$diel3,
      habitat = newdat$habitat
    ) %>%
    pivot_longer(
      cols = all_of(state_levels),
      names_to = "state_m3_w",
      values_to = "SE"
    )
  
  # Combine
  out <- left_join(
    fit_long,
    se_long,
    by = c("diel3", "habitat", "state_m3_w")
  )
  
  # CI
  z <- qnorm(1 - (1 - conf.level) / 2)
  
  out %>%
    mutate(
      eta = qlogis(estimate),
      SE_eta = SE / (estimate * (1 - estimate)),
      lower.CL = plogis(eta - z * SE_eta),
      upper.CL = plogis(eta + z * SE_eta)
    ) %>%
    select(
      diel3,
      habitat,
      state_m3_w,
      estimate,
      SE,
      lower.CL,
      upper.CL
    )
}

# final predictions
pred_daynight <- mblogit_predictions(
  mod_daynight,
  dat_mod
)

pred_daynight

# value check
pred_daynight %>%
  group_by(diel3, habitat) %>%
  summarise(
    total = sum(estimate)
  )


### ---------------- ###
# FIGURE - COMMON THEME
### ---------------- ###

figure_theme <- theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      size = 11,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 11,
      colour = "black"
    ),
    
    axis.title = element_text(
      size = 14,
      face = "bold"
    ),
    
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    
    legend.position = "none",
    
    plot.margin = margin(
      t = 10,
      r = 15,
      b = 10,
      l = 10
    )
  )

### ---------------- ###
# FIGURE - BEHAVIOR BY HABITAT
### ---------------- ###

##### 1 - labels and ordering
habitat_labels <- c(
  "open"     = "Open water",
  "oyster"   = "Eroded oyster reef",
  "SAV"      = "SAV"
)

habitat_order <- c(
  "Open water",
  "Eroded oyster reef",
  "SAV"
)

###### 2 - population level predictions
pred_plot <- pred_daynight %>%
  mutate(
    habitat = recode(habitat, !!!habitat_labels),
    
    diel = recode(
      diel3,
      "day" = "Day",
      "night" = "Night"
    ),
    
    behavior = recode(
      state_m3_w,
      "1" = "Resting",
      "2" = "Searching",
      "3" = "Travelling"
    )
  ) %>%
  mutate(
    habitat = factor(
      habitat,
      levels = habitat_order
    ),
    
    diel = factor(
      diel,
      levels = c("Day", "Night")
    ),
    
    behavior = factor(
      behavior,
      levels = c(
        "Resting",
        "Searching",
        "Travelling"
      )
    )
  )


##### 3 - colors
diel_colors <- c(
  "Day" = "#E69F00",
  "Night" = "#163A5F"
)


##### 4 - figure
figure_habitat <- pred_plot %>%
  ggplot(
    aes(
      x = behavior,
      y = estimate,
      colour = diel
    )
  ) +
  
  geom_errorbar(
    aes(
      ymin = lower.CL,
      ymax = upper.CL
    ),
    position = position_dodge(width = 0.35),
    width = 0.10,
    linewidth = 1.3
  ) +
  
  geom_point(
    position = position_dodge(width = 0.35),
    size = 4.5
  ) +
  
  scale_colour_manual(
    values = diel_colors
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1),
    expand = c(0, 0)
  ) +
  
  facet_wrap(
    ~ habitat,
    ncol = 2
  ) +
  
  labs(
    x = "Behavior",
    y = "Probability of behavior",
    colour = "Diel period"
  ) +
  
  figure_theme +
  
  theme(
    legend.position = "none",
    strip.text = element_text(
      size = 14,
      face = "bold",
      colour = "black"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


figure_habitat

