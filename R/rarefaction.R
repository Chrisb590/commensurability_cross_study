library(bipartite)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(tibble)
library(here)

################################################################################
#
# Data
#
################################################################################

# Directory where all figures and tables are written (created if missing)
output_dir <- here("output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

rare <- read_csv(
  here("data","rarefaction_sensitivity_500_iterations.csv")
)

metadata <- read_csv(
  here("data","general_network_information.csv")
) %>%
  mutate(
    Publication = as.character(Publication)
  )

################################################################################
#
# Main manuscript analyses:
# Rarefaction analysis of quantitative interaction networks
#
# Evaluate whether differences in quantitative network structure within and
# between publications persist after standardizing networks to common sampling
# depths. Analyses are conducted separately for pollination and seed-dispersal
# networks at rarefaction depths of 50, 100, and 200 interaction events. Results
# reported in Table 3 use a rarefaction depth of 100 interaction events.
#
# For each network and sampling depth, network indices are averaged across
# 500 rarefaction iterations before publication-level variation is calculated.
#
################################################################################

################################################################################
#
# Average network indices across rarefaction iterations
#
# For each network and rarefaction depth, calculate the mean specialization,
# nestedness, and modularity across the 500 rarefaction iterations. These mean
# values represent the rarefied network-level estimates used in subsequent
# publication-level comparisons.
#
################################################################################

rare_network_mean <- rare %>%
  filter(
    depth %in% c(
      50,
      100,
      200
    )
  ) %>%
  group_by(
    ID,
    depth
  ) %>%
  summarise(
    H2 = mean(
      H2,
      na.rm = TRUE
    ),
    weighted_NODF = mean(
      weighted_NODF,
      na.rm = TRUE
    ),
    weighted_modularity = mean(
      weighted_modularity,
      na.rm = TRUE
    ),
    n_iterations = n(),
    .groups = "drop"
  )

################################################################################
#
# Define publication groups at each rarefaction depth
#
# A publication contributes to the within-publication comparison only when at
# least two of its networks remain represented at a given rarefaction depth.
#
# If only one network from a publication remains at that depth, it is reassigned
# to the "One network per publication" group. This pooled group therefore
# contains one network from each contributing publication and is used to
# characterize variation between publications.
#
################################################################################

rare_network_mean <- rare_network_mean %>%
  inner_join(
    metadata %>%
      select(
        ID,
        Publication,
        TYPE
      ),
    by = "ID"
  ) %>%
  group_by(
    depth,
    TYPE,
    Publication
  ) %>%
  mutate(
    n_networks_at_depth =
      n_distinct(ID)
  ) %>%
  ungroup() %>%
  mutate(
    Publication = if_else(
      Publication ==
        "One network per publication" |
        n_networks_at_depth == 1,
      "One network per publication",
      Publication
    )
  )

################################################################################
#
# Reshape network indices for publication-level analyses
#
# Convert specialization, nestedness, and modularity to long format so that
# publication-level means and standard deviations can be calculated using the
# same workflow for each network metric.
#
################################################################################

rare_long <- rare_network_mean %>%
  select(
    ID,
    depth,
    Publication,
    TYPE,
    H2,
    weighted_NODF,
    weighted_modularity
  ) %>%
  pivot_longer(
    cols = c(
      H2,
      weighted_NODF,
      weighted_modularity
    ),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = case_when(
      Metric == "H2" ~
        "Specialization (H2)",
      Metric == "weighted_NODF" ~
        "Nestedness (weighted NODF)",
      Metric == "weighted_modularity" ~
        "Modularity",
      TRUE ~
        Metric
    )
  )

################################################################################
#
# Calculate variation within publication groups
#
# For each rarefaction depth, network type, publication group, and metric,
# calculate the mean and SD among networks.
#
# For publications contributing multiple networks, the SD measures variation
# among networks from the same publication. For the pooled "One network per
# publication" group, the SD measures variation among networks originating
# from different publications.
#
################################################################################

rare_pub_stats <- rare_long %>%
  group_by(
    depth,
    TYPE,
    Publication,
    Metric
  ) %>%
  summarise(
    mean = mean(
      Value,
      na.rm = TRUE
    ),
    sd_val = sd(
      Value,
      na.rm = TRUE
    ),
    n = sum(
      !is.na(Value)
    ),
    .groups = "drop"
  )

################################################################################
#
# Tables 3, S4, and S5:
# Compare variation within and between publications after rarefaction
#
# Table 3: rarefaction depth = 100 interactions
# Table S4: rarefaction depth = 50 interactions
# Table S5: rarefaction depth = 200 interactions
#
# For networks originating from different publications, use the SD among
# networks in the pooled "One network per publication" group.
#
# For publications contributing multiple networks, first calculate the SD
# separately within each publication and then average these publication-level
# SDs. Each publication therefore contributes equally to the within-publication
# estimate regardless of how many networks it contains.
#
################################################################################

# Between-publication variation
summary_one_net <- rare_pub_stats %>%
  filter(
    Publication ==
      "One network per publication"
  ) %>%
  select(
    depth,
    TYPE,
    Metric,
    sd_val,
    n
  ) %>%
  rename(
    SD_One_Network_Per_Publication =
      sd_val,
    n_One_Network_Per_Publication =
      n
  )

# Mean within-publication variation
summary_other_unweighted <- rare_pub_stats %>%
  filter(
    Publication !=
      "One network per publication"
  ) %>%
  group_by(
    depth,
    TYPE,
    Metric
  ) %>%
  summarise(
    Avg_SD_Other =
      mean(
        sd_val,
        na.rm = TRUE
      ),
    n_publications_other =
      sum(
        !is.na(sd_val)
      ),
    n_networks_other =
      sum(
        n[
          !is.na(sd_val)
        ]
      ),
    .groups = "drop"
  )

# Combine within- and between-publication estimates
rarefaction_sd_comparison <-
  summary_other_unweighted %>%
  left_join(
    summary_one_net,
    by = c(
      "depth",
      "TYPE",
      "Metric"
    )
  ) %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(
        .x,
        2
      )
    )
  )

################################################################################
#
# Table 3: 100 interactions
#
################################################################################

table_3 <- rarefaction_sd_comparison %>%
  filter(
    depth == 100
  ) %>%
  select(
    -depth
  )

print(
  table_3
)

# Save table
write_csv(table_3, file.path(output_dir, "table_3.csv"))

################################################################################
#
# Table S4: 50 interactions
#
################################################################################

table_s4 <- rarefaction_sd_comparison %>%
  filter(
    depth == 50
  ) %>%
  select(
    -depth
  )

print(
  table_s4
)

# Save table
write_csv(table_s4, file.path(output_dir, "table_s4.csv"))

################################################################################
#
# Table S5: 200 interactions
#
################################################################################

table_s5 <- rarefaction_sd_comparison %>%
  filter(
    depth == 200
  ) %>%
  select(
    -depth
  )

print(
  table_s5
)

# Save table
write_csv(table_s5, file.path(output_dir, "table_s5.csv"))

################################################################################
#
# Table S3:
# Publication-level network indices after rarefaction to 100 interactions
#
# Report the mean and SD of specialization, nestedness, and modularity for
# each publication, using network-level indices averaged across 500
# rarefaction iterations.
#
# Also calculate unweighted averages across publications contributing
# multiple networks. Each publication contributes equally to these averages.
#
################################################################################

# Reshape publication-level statistics into a wide table
publication_table <- rare_pub_stats %>%
  filter(
    depth == 100
  ) %>%
  mutate(
    Metric = recode(
      Metric,
      "Specialization (H2)" = "H2",
      "Nestedness (weighted NODF)" = "wNODF",
      "Modularity" = "Modularity"
    )
  ) %>%
  select(
    TYPE,
    Publication,
    Metric,
    mean,
    sd_val
  ) %>%
  pivot_wider(
    names_from = Metric,
    values_from = c(
      mean,
      sd_val
    ),
    names_glue = "{Metric}_{.value}"
  )

# Count the number of networks in each publication group
publication_counts <- rare_network_mean %>%
  filter(
    depth == 100
  ) %>%
  count(
    TYPE,
    Publication,
    name = "n"
  )

# Combine statistics and network counts
publication_table <- publication_table %>%
  left_join(
    publication_counts,
    by = c(
      "TYPE",
      "Publication"
    )
  ) %>%
  mutate(
    Type = if_else(
      TYPE == "Pollination",
      "PL",
      "SD"
    )
  ) %>%
  select(
    Type,
    Publication,
    n,
    H2_mean,
    H2_sd_val,
    wNODF_mean,
    wNODF_sd_val,
    Modularity_mean,
    Modularity_sd_val
  )

# Separate publications contributing multiple networks
multiple_publications <- publication_table %>%
  filter(
    Publication != "One network per publication"
  ) %>%
  arrange(
    Type,
    Publication
  )

# Calculate unweighted averages across publications
average_publications <- multiple_publications %>%
  group_by(
    Type
  ) %>%
  summarise(
    Publication = paste0(
      "Average unweighted ",
      n(),
      " ",
      first(Type),
      " publications above"
    ),
    n = NA_integer_,
    across(
      c(
        H2_mean,
        H2_sd_val,
        wNODF_mean,
        wNODF_sd_val,
        Modularity_mean,
        Modularity_sd_val
      ),
      ~ mean(
        .x,
        na.rm = TRUE
      )
    ),
    .groups = "drop"
  )

# Separate publications contributing one network
single_publications <- publication_table %>%
  filter(
    Publication == "One network per publication"
  )

# Assemble Table S3 in manuscript order
table_s3 <- bind_rows(
  multiple_publications %>%
    filter(Type == "PL"),
  
  average_publications %>%
    filter(Type == "PL"),
  
  multiple_publications %>%
    filter(Type == "SD"),
  
  average_publications %>%
    filter(Type == "SD"),
  
  single_publications %>%
    filter(Type == "PL"),
  
  single_publications %>%
    filter(Type == "SD")
) %>%
  mutate(
    across(
      where(is.numeric) & !all_of("n"),
      ~ round(
        .x,
        2
      )
    )
  )

print(
  table_s3,
  n = Inf,
  width = Inf
)

# Save table
write_csv(
  table_s3,
  file.path(output_dir, "table_s3.csv")
)