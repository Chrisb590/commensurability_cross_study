library(lme4)
library(lmerTest)
library(MuMIn)
library(here)
library(readr)
library(dplyr)
library(stringr)

################################################################################
#
# Data
#
################################################################################

meta <- read_csv(here("general_network_information.csv"), show_col_types = FALSE) %>%
  mutate(
    ID = trimws(ID),
    TYPE = trimws(TYPE) %>%
      tolower() %>%
      str_replace_all(" ", "-") %>%
      str_to_title()
  )

# Since H2', NODF, weighted modularity is already calculated in Patefield 
# randomization, I can just read them in
indices <- read_csv(here("all_results_for_patefield_randomization.csv"), show_col_types = FALSE) %>%
  transmute(
    ID = str_trim(ID),
    H2 = as.numeric(H2),
    weighted_NODF = as.numeric(weighted_NODF),
    weighted_modularity = as.numeric(DIRTMod)  
  )

results <- meta %>%
  left_join(indices,by="ID")

results$Publication <- factor(results$Publication)

# Only pollination networks from publications that produced multiple networks
results_no_one_pollination <- results %>%
  filter(Publication != "One network per publication" & TYPE == "Pollination")

# Only seed-dispersal networks from publications that produced multiple networks
results_no_one_seed <- results %>%
  filter(Publication != "One network per publication" & TYPE == "Seed-Dispersal")

################################################################################
#
# Table 2: linear model vs linear mixed model
#
################################################################################

# Pollination

# Specialization linear mixed model
mod_H2 <- lmer(H2 ~ Latitude + (1 | Publication), data = results_no_one_pollination)
summary(mod_H2)
r.squaredGLMM(mod_H2)

# Specialization linear model
mod_H2 <- lm(H2 ~ Latitude, data = results_no_one_pollination)
summary(mod_H2)

# Nestedness linear mixed model
mod_wNODF <- lmer(weighted_NODF ~ Latitude  + (1 | Publication), data = results_no_one_pollination)
summary(mod_wNODF)
r.squaredGLMM(mod_wNODF)

# Nestedness linear model
mod_wNODF <- lm(weighted_NODF ~ Latitude , data = results_no_one_pollination)
summary(mod_wNODF)

# Modularity linear mixed model
mod_mod <- lmer(weighted_modularity ~ Latitude  + (1 | Publication), data = results_no_one_pollination)
summary(mod_mod)
r.squaredGLMM(mod_mod)

# Modularity linear model
mod_mod <- lm(weighted_modularity ~ Latitude, data = results_no_one_pollination)
summary(mod_mod)

# Seed-disperseal

# Specialization linear mixed model
mod_H2 <- lmer(H2 ~ Latitude + (1 | Publication), data = results_no_one_seed)
summary(mod_H2)
r.squaredGLMM(mod_H2)

# Specialization linear model
mod_H2 <- lm(H2 ~ Latitude, data = results_no_one_seed)
summary(mod_H2)

# Nestedness linear mixed model
mod_wNODF <- lmer(weighted_NODF ~ Latitude  + (1 | Publication), data = results_no_one_seed)
summary(mod_wNODF)
r.squaredGLMM(mod_wNODF)

# Nestedness linear model
mod_wNODF <- lm(weighted_NODF ~ Latitude , data = results_no_one_seed)
summary(mod_wNODF)

# Modularity linear mixed model
mod_mod <- lmer(weighted_modularity ~ Latitude  + (1 | Publication), data = results_no_one_seed)
summary(mod_mod)
r.squaredGLMM(mod_mod)

# Modularity linear model
mod_mod <- lm(weighted_modularity ~ Latitude, data = results_no_one_seed)
summary(mod_mod)

################################################################################
#
# Average variance in latitude per publication
#
################################################################################

results_no_one <- results %>%
  filter(Publication != "One network per publication")

lat_summary <- results_no_one %>%
  group_by(Publication) %>%
  summarise(
    n_networks = n(),
    lat_range = max(Latitude) - min(Latitude),
    .groups = "drop"
  )

mean(lat_summary$lat_range)
