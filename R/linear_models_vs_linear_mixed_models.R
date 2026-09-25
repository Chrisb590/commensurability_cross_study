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

# Directory where all figures and tables are written (created if missing)
output_dir <- here("output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

meta <- read_csv(here("data", "general_network_information.csv"), show_col_types = FALSE) %>%
  mutate(
    ID = trimws(ID),
    TYPE = trimws(TYPE) %>%
      tolower() %>%
      str_replace_all(" ", "-") %>%
      str_to_title()
  )

# Since H2', NODF, weighted modularity is already calculated in Patefield 
# randomization, I can just read them in
indices <- read_csv(here("data", "all_results_for_patefield_randomization.csv"), show_col_types = FALSE) %>%
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

# All model summaries are written to a single text file
sink(file.path(output_dir, "table_2_linear_models_vs_linear_mixed_models.txt"))

# Pollination
cat("\n==== Pollination ====\n")

# Specialization linear mixed model
cat("\n--- Specialization linear mixed model ---\n")
mod_H2 <- lmer(H2 ~ Latitude + (1 | Publication), data = results_no_one_pollination)
print(summary(mod_H2))
print(r.squaredGLMM(mod_H2))

# Specialization linear model
cat("\n--- Specialization linear model ---\n")
mod_H2 <- lm(H2 ~ Latitude, data = results_no_one_pollination)
print(summary(mod_H2))

# Nestedness linear mixed model
cat("\n--- Nestedness linear mixed model ---\n")
mod_wNODF <- lmer(weighted_NODF ~ Latitude  + (1 | Publication), data = results_no_one_pollination)
print(summary(mod_wNODF))
print(r.squaredGLMM(mod_wNODF))

# Nestedness linear model
cat("\n--- Nestedness linear model ---\n")
mod_wNODF <- lm(weighted_NODF ~ Latitude , data = results_no_one_pollination)
print(summary(mod_wNODF))

# Modularity linear mixed model
cat("\n--- Modularity linear mixed model ---\n")
mod_mod <- lmer(weighted_modularity ~ Latitude  + (1 | Publication), data = results_no_one_pollination)
print(summary(mod_mod))
print(r.squaredGLMM(mod_mod))

# Modularity linear model
cat("\n--- Modularity linear model ---\n")
mod_mod <- lm(weighted_modularity ~ Latitude, data = results_no_one_pollination)
print(summary(mod_mod))

# Seed-disperseal
cat("\n==== Seed-disperseal ====\n")

# Specialization linear mixed model
cat("\n--- Specialization linear mixed model ---\n")
mod_H2 <- lmer(H2 ~ Latitude + (1 | Publication), data = results_no_one_seed)
print(summary(mod_H2))
print(r.squaredGLMM(mod_H2))

# Specialization linear model
cat("\n--- Specialization linear model ---\n")
mod_H2 <- lm(H2 ~ Latitude, data = results_no_one_seed)
print(summary(mod_H2))

# Nestedness linear mixed model
cat("\n--- Nestedness linear mixed model ---\n")
mod_wNODF <- lmer(weighted_NODF ~ Latitude  + (1 | Publication), data = results_no_one_seed)
print(summary(mod_wNODF))
print(r.squaredGLMM(mod_wNODF))

# Nestedness linear model
cat("\n--- Nestedness linear model ---\n")
mod_wNODF <- lm(weighted_NODF ~ Latitude , data = results_no_one_seed)
print(summary(mod_wNODF))

# Modularity linear mixed model
cat("\n--- Modularity linear mixed model ---\n")
mod_mod <- lmer(weighted_modularity ~ Latitude  + (1 | Publication), data = results_no_one_seed)
print(summary(mod_mod))
print(r.squaredGLMM(mod_mod))

# Modularity linear model
cat("\n--- Modularity linear model ---\n")
mod_mod <- lm(weighted_modularity ~ Latitude, data = results_no_one_seed)
print(summary(mod_mod))

################################################################################
#
# Average variance in latitude per publication
#
################################################################################

cat("\n==== Latitude range per publication ====\n")

results_no_one <- results %>%
  filter(Publication != "One network per publication")

lat_summary <- results_no_one %>%
  group_by(Publication) %>%
  summarise(
    n_networks = n(),
    lat_range = max(Latitude) - min(Latitude),
    .groups = "drop"
  )

print(lat_summary, n = Inf)
print(mean(lat_summary$lat_range))

sink()
write_csv(lat_summary, file.path(output_dir, "latitude_range_per_publication.csv"))
