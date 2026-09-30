library(lme4)
library(lmerTest)
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

meta <- read_csv(
  here("data","general_network_information.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    ID = str_trim(ID),
    TYPE = TYPE %>%
      str_trim() %>%
      str_to_lower() %>%
      str_replace_all(" ", "-") %>%
      str_to_title()
  )

# H2', weighted NODF, and weighted modularity were already calculated
# during the Patefield randomization analysis
indices <- read_csv(
  here("data","all_results_for_patefield_randomization.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    ID = str_trim(ID),
    H2 = as.numeric(H2),
    weighted_NODF = as.numeric(weighted_NODF),
    weighted_modularity = as.numeric(DIRTMod)
  )

results <- meta %>%
  left_join(
    indices,
    by = "ID"
  ) %>%
  mutate(
    Publication = factor(Publication)
  )

################################################################################
#
# Main manuscript analyses:
# Latitude effects on quantitative interaction networks
#
# Evaluate whether relationships between latitude and quantitative network
# structure change after accounting for non-independence among networks from
# the same publication. Analyses are conducted separately for pollination and
# seed-dispersal networks and for specialization, nestedness, and modularity.
#
# Linear models (LMs) and linear mixed models (LMMs) are fitted to the same
# subset of networks, restricted to publications contributing multiple networks.
# LMMs include publication as a random intercept to account for shared
# publication-level variation. Results are reported in Table 2.
#
################################################################################

################################################################################
#
# Restrict analyses to publications contributing multiple networks
#
# Networks from publications contributing only one network cannot contribute
# information about within-publication dependence. These networks are therefore
# excluded so that the LM and LMM for each network type are fitted to the same
# subset of data.
#
################################################################################

results_no_one <- results %>%
  filter(
    Publication !=
      "One network per publication"
  )

results_no_one_pollination <- results_no_one %>%
  filter(
    TYPE == "Pollination"
  )

results_no_one_seed <- results_no_one %>%
  filter(
    TYPE == "Seed-Dispersal"
  )

################################################################################
#
# Extract latitude effects from fitted models
#
# For each LM and LMM, extract the estimated latitude slope, 95% confidence
# interval, and p-value. 
#
################################################################################

get_lm_latitude <- function(model) {
  
  x <- summary(model)$coefficients["Latitude", ]
  
  estimate <- x["Estimate"]
  se <- x["Std. Error"]
  df <- df.residual(model)
  critical_value <- qt(
    0.975,
    df = df
  )
  
  data.frame(
    beta = estimate,
    CI_lower = estimate - critical_value * se,
    CI_upper = estimate + critical_value * se,
    p = x["Pr(>|t|)"]
  )
}

get_lmm_latitude <- function(model) {
  
  x <- coef(summary(model))["Latitude", ]
  
  estimate <- x["Estimate"]
  se <- x["Std. Error"]
  df <- x["df"]
  critical_value <- qt(
    0.975,
    df = df
  )
  
  data.frame(
    beta = estimate,
    CI_lower = estimate - critical_value * se,
    CI_upper = estimate + critical_value * se,
    p = x["Pr(>|t|)"]
  )
}

################################################################################
#
# Fit pollination models
#
# For each quantitative network index, fit:
#
#   1. LM:  Index ~ Latitude
#   2. LMM: Index ~ Latitude + (1 | Publication)
#
# Comparing these models shows how the estimated latitude relationship changes
# after accounting for non-independence among networks from the same publication.
#
################################################################################

# Specialization
mod_H2_lm_poll <- lm(
  H2 ~ Latitude,
  data = results_no_one_pollination
)

mod_H2_lmm_poll <- lmer(
  H2 ~ Latitude + (1 | Publication),
  data = results_no_one_pollination
)

# Nestedness
mod_wNODF_lm_poll <- lm(
  weighted_NODF ~ Latitude,
  data = results_no_one_pollination
)

mod_wNODF_lmm_poll <- lmer(
  weighted_NODF ~ Latitude + (1 | Publication),
  data = results_no_one_pollination
)

# Modularity
mod_mod_lm_poll <- lm(
  weighted_modularity ~ Latitude,
  data = results_no_one_pollination
)

mod_mod_lmm_poll <- lmer(
  weighted_modularity ~ Latitude + (1 | Publication),
  data = results_no_one_pollination
)

################################################################################
#
# Fit seed-dispersal models
#
# Fit the same LM and LMM specifications used for pollination networks to
# specialization, nestedness, and modularity of seed-dispersal networks.
#
################################################################################

# Specialization
mod_H2_lm_seed <- lm(
  H2 ~ Latitude,
  data = results_no_one_seed
)

mod_H2_lmm_seed <- lmer(
  H2 ~ Latitude + (1 | Publication),
  data = results_no_one_seed
)

# Nestedness
mod_wNODF_lm_seed <- lm(
  weighted_NODF ~ Latitude,
  data = results_no_one_seed
)

mod_wNODF_lmm_seed <- lmer(
  weighted_NODF ~ Latitude + (1 | Publication),
  data = results_no_one_seed
)

# Modularity
mod_mod_lm_seed <- lm(
  weighted_modularity ~ Latitude,
  data = results_no_one_seed
)

mod_mod_lmm_seed <- lmer(
  weighted_modularity ~ Latitude + (1 | Publication),
  data = results_no_one_seed
)

################################################################################
#
# Table 2:
# Latitude effects from linear and linear mixed models
#
# Combine the latitude slopes, 95% confidence intervals, and p-values from the
# LM and LMM for each quantitative network index and network type.
#
################################################################################

latitude_results <- bind_rows(
  
  # Pollination
  bind_cols(
    Type = "Pollination",
    Index = "Specialization",
    Model = "LM",
    get_lm_latitude(
      mod_H2_lm_poll
    )
  ),
  
  bind_cols(
    Type = "Pollination",
    Index = "Specialization",
    Model = "LMM",
    get_lmm_latitude(
      mod_H2_lmm_poll
    )
  ),
  
  bind_cols(
    Type = "Pollination",
    Index = "Nestedness",
    Model = "LM",
    get_lm_latitude(
      mod_wNODF_lm_poll
    )
  ),
  
  bind_cols(
    Type = "Pollination",
    Index = "Nestedness",
    Model = "LMM",
    get_lmm_latitude(
      mod_wNODF_lmm_poll
    )
  ),
  
  bind_cols(
    Type = "Pollination",
    Index = "Modularity",
    Model = "LM",
    get_lm_latitude(
      mod_mod_lm_poll
    )
  ),
  
  bind_cols(
    Type = "Pollination",
    Index = "Modularity",
    Model = "LMM",
    get_lmm_latitude(
      mod_mod_lmm_poll
    )
  ),
  
  # Seed-dispersal
  bind_cols(
    Type = "Seed-dispersal",
    Index = "Specialization",
    Model = "LM",
    get_lm_latitude(
      mod_H2_lm_seed
    )
  ),
  
  bind_cols(
    Type = "Seed-dispersal",
    Index = "Specialization",
    Model = "LMM",
    get_lmm_latitude(
      mod_H2_lmm_seed
    )
  ),
  
  bind_cols(
    Type = "Seed-dispersal",
    Index = "Nestedness",
    Model = "LM",
    get_lm_latitude(
      mod_wNODF_lm_seed
    )
  ),
  
  bind_cols(
    Type = "Seed-dispersal",
    Index = "Nestedness",
    Model = "LMM",
    get_lmm_latitude(
      mod_wNODF_lmm_seed
    )
  ),
  
  bind_cols(
    Type = "Seed-dispersal",
    Index = "Modularity",
    Model = "LM",
    get_lm_latitude(
      mod_mod_lm_seed
    )
  ),
  
  bind_cols(
    Type = "Seed-dispersal",
    Index = "Modularity",
    Model = "LMM",
    get_lmm_latitude(
      mod_mod_lmm_seed
    )
  )
)

latitude_results

# Save table
write_csv(latitude_results, file.path(output_dir, "table_2.csv"))
