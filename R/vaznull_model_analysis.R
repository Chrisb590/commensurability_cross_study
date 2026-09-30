library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(purrr)
library(readr)
library(patchwork)
library(here)

################################################################################
#
# Data
#
################################################################################

# Directory where all figures and tables are written (created if missing)
output_dir <- here("output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

results_df <- read_csv(
  here("data","all_results_vaznull_randomizations.csv"),
  show_col_types = FALSE
)

metadata <- read_csv(
  here("data","general_network_information.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    Publication =
      as.character(
        Publication
      )
  )

# Convert stored Vázquez null-model values to numeric vectors
results_df <- results_df %>%
  mutate(
    
    H2_null_values =
      str_split(
        H2_null_values,
        ";"
      ) %>%
      map(
        ~ as.numeric(.x)
      ),
    
    weighted_NODF_null_values =
      str_split(
        weighted_NODF_null_values,
        ";"
      ) %>%
      map(
        ~ as.numeric(.x)
      ),
    
    DIRTMod_null_values =
      str_split(
        DIRTMod_null_values,
        ";"
      ) %>%
      map(
        ~ as.numeric(.x)
      )
  )

# Join Vázquez results with network metadata
results_df <- results_df %>%
  inner_join(
    metadata,
    by = "ID"
  )

################################################################################
#
# Main manuscript analyses:
# Vázquez null-model correction of quantitative interaction networks
#
# Evaluate whether publication-specific variation in quantitative network
# topology persists after correction using the Vázquez null model. Analyses
# are conducted for specialization, nestedness, and modularity using 1,000
# Vázquez randomizations for each empirical network.
#
# Figure S5 compares empirical network indices with their Vázquez null-model
# distributions. Table 5 compares variation in Delta Vázquez indices between
# the "One network per publication" and "Multiple networks per publication"
# groups, with publication-level results reported in Table S9. The corresponding
# analyses using z Vázquez indices are reported in Tables S10 and S11.
#
################################################################################

################################################################################
#
# Prepare empirical and Vázquez null-model values for Figure S5
#
# Separate the "One network per publication" group by interaction type for
# plotting. Empirical indices and all 1,000 Vázquez null-model values are
# then converted to a common long format.
#
################################################################################

results_df <- results_df %>%
  mutate(
    Publication = case_when(
      
      Publication ==
        "One network per publication" &
        TYPE ==
        "Pollination" ~
        "One network per\npublication (Pollination)",
      
      Publication ==
        "One network per publication" &
        TYPE ==
        "Seed Dispersal" ~
        "One network per\npublication (Seed Dispersal)",
      
      TRUE ~
        Publication
    )
  )

# Place the "One network per publication" groups first
ordered_publications <- c(
  
  "One network per\npublication (Pollination)",
  "One network per\npublication (Seed Dispersal)",
  
  sort(
    setdiff(
      unique(
        results_df$Publication
      ),
      c(
        "One network per\npublication (Pollination)",
        "One network per\npublication (Seed Dispersal)"
      )
    )
  )
)

results_df$Publication <- factor(
  results_df$Publication,
  levels =
    ordered_publications
)

# Empirical network indices
empirical_long <- results_df %>%
  select(
    ID,
    Publication,
    H2,
    weighted_NODF,
    DIRTMod
  ) %>%
  
  pivot_longer(
    cols = c(
      H2,
      weighted_NODF,
      DIRTMod
    ),
    names_to =
      "metric",
    values_to =
      "value"
  ) %>%
  
  mutate(
    Type =
      "Empirical",
    
    Metric = case_when(
      
      metric ==
        "H2" ~
        "Specialization (H2)",
      
      metric ==
        "weighted_NODF" ~
        "Nestedness (weighted NODF)",
      
      metric ==
        "DIRTMod" ~
        "Modularity (DIRTMod)"
    )
  ) %>%
  
  select(
    ID,
    Publication,
    value,
    Type,
    Metric
  )

# Vázquez null-model indices
null_long <- results_df %>%
  select(
    ID,
    Publication,
    H2_null_values,
    weighted_NODF_null_values,
    DIRTMod_null_values
  ) %>%
  
  pivot_longer(
    cols = c(
      H2_null_values,
      weighted_NODF_null_values,
      DIRTMod_null_values
    ),
    names_to =
      "metric",
    values_to =
      "value_list"
  ) %>%
  
  unnest(
    value_list
  ) %>%
  
  mutate(
    value =
      value_list,
    
    Type =
      "Vázquez null",
    
    Metric = case_when(
      
      metric ==
        "H2_null_values" ~
        "Specialization (H2)",
      
      metric ==
        "weighted_NODF_null_values" ~
        "Nestedness (weighted NODF)",
      
      metric ==
        "DIRTMod_null_values" ~
        "Modularity (DIRTMod)"
    )
  ) %>%
  
  select(
    ID,
    Publication,
    value,
    Type,
    Metric
  )

# Combine empirical and Vázquez null-model indices
df_long <- bind_rows(
  empirical_long,
  null_long
) %>%
  filter(
    !is.na(
      Publication
    )
  )

df_specialization <- df_long %>%
  filter(
    Metric ==
      "Specialization (H2)"
  )

df_nestedness <- df_long %>%
  filter(
    Metric ==
      "Nestedness (weighted NODF)"
  )

df_modularity <- df_long %>%
  filter(
    Metric ==
      "Modularity (DIRTMod)"
  )

# Assign publication-label colours according to network type
publication_types <- results_df %>%
  distinct(
    Publication,
    TYPE
  ) %>%
  mutate(
    label_color =
      if_else(
        TYPE ==
          "Pollination",
        "#e6550d",
        "#41ab5d"
      )
  )

df_specialization <- df_specialization %>%
  left_join(
    publication_types,
    by = "Publication"
  )

df_nestedness <- df_nestedness %>%
  left_join(
    publication_types,
    by = "Publication"
  )

df_modularity <- df_modularity %>%
  left_join(
    publication_types,
    by = "Publication"
  )

################################################################################
#
# Figure S5:
# Empirical networks and their Vázquez null-model replicates
#
# Compare quantitative topological indices for the 279 empirical networks with
# their 1,000 Vázquez null-model replicates, grouped by publication source.
# Text colour distinguishes pollination and seed-dispersal networks.
#
################################################################################

# Specialization
p1 <- ggplot(
  df_specialization,
  aes(
    x =
      Publication,
    y =
      value,
    fill =
      Type
  )
) +
  
  geom_boxplot(
    outlier.size =
      0.6,
    alpha =
      0.8,
    width =
      0.6,
    color =
      "grey30",
    linewidth =
      0.3
  ) +
  
  geom_vline(
    xintercept =
      2.5,
    linetype =
      "dashed",
    color =
      "black",
    linewidth =
      0.5
  ) +
  
  ylim(
    0,
    1
  ) +
  
  labs(
    x =
      NULL,
    y =
      "Specialization (H2')",
    fill =
      NULL
  ) +
  
  scale_fill_manual(
    values = c(
      "Empirical" =
        "#fcae91",
      "Vázquez null" =
        "#9972af"
    )
  ) +
  
  theme_minimal(
    base_size =
      12
  ) +
  
  theme(
    axis.text.x =
      element_blank(),
    
    axis.ticks.x =
      element_blank(),
    
    panel.border =
      element_rect(
        color = "grey50",
        fill = NA,
        linewidth = 0.8
      )
  )

# Nestedness
p2 <- ggplot(
  df_nestedness,
  aes(
    x =
      Publication,
    y =
      value,
    fill =
      Type
  )
) +
  
  geom_boxplot(
    outlier.size =
      0.6,
    alpha =
      0.8,
    width =
      0.6,
    color =
      "grey30",
    linewidth =
      0.3
  ) +
  
  geom_vline(
    xintercept =
      2.5,
    linetype =
      "dashed",
    color =
      "black",
    linewidth =
      0.5
  ) +
  
  ylim(
    0,
    100
  ) +
  
  labs(
    x =
      NULL,
    y =
      "Nestedness (wNODF)",
    fill =
      NULL
  ) +
  
  scale_fill_manual(
    values = c(
      "Empirical" =
        "#fcae91",
      "Vázquez null" =
        "#9972af"
    )
  ) +
  
  theme_minimal(
    base_size =
      12
  ) +
  
  theme(
    axis.text.x =
      element_blank(),
    
    axis.ticks.x =
      element_blank(),
    
    panel.border =
      element_rect(
        color = "grey50",
        fill = NA,
        linewidth = 0.8
      )
  )

# Modularity
p3 <- ggplot(
  df_modularity,
  aes(
    x =
      Publication,
    y =
      value,
    fill =
      Type
  )
) +
  
  geom_boxplot(
    outlier.size =
      0.6,
    alpha =
      0.8,
    width =
      0.6,
    color =
      "grey30",
    linewidth =
      0.3
  ) +
  
  geom_vline(
    xintercept =
      2.5,
    linetype =
      "dashed",
    color =
      "black",
    linewidth =
      0.5
  ) +
  
  ylim(
    0,
    1
  ) +
  
  labs(
    x =
      NULL,
    y =
      "Modularity (DIRTLPAwb+)",
    fill =
      NULL
  ) +
  
  scale_fill_manual(
    values = c(
      "Empirical" =
        "#fcae91",
      "Vázquez null" =
        "#9972af"
    )
  ) +
  
  theme_minimal(
    base_size =
      12
  ) +
  
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1,
        size = 10,
        color =
          df_modularity %>%
          distinct(
            Publication,
            label_color
          ) %>%
          arrange(
            factor(
              Publication,
              levels =
                levels(
                  df_modularity$Publication
                )
            )
          ) %>%
          pull(
            label_color
          )
      ),
    
    panel.border =
      element_rect(
        color = "grey50",
        fill = NA,
        linewidth = 0.8
      )
  )

# Combine specialization, nestedness, and modularity panels
p_combined <- (
  p1 /
    p2 /
    p3 +
    plot_layout(
      guides =
        "collect"
    )
) &
  theme(
    legend.position =
      "top"
  )

print(
  p_combined
)

ggsave(
  file.path(output_dir,"figure_s5.pdf"),
  plot =
    p_combined,
  width =
    13,
  height =
    11
)

################################################################################
#
# Calculate Delta Vázquez and z Vázquez indices
#
# For each empirical network, calculate the mean and SD of its 1,000 Vázquez
# null-model values. Delta Vázquez is the empirical network index minus the
# mean of its null distribution. z Vázquez standardizes this difference by
# the SD of the null distribution.
#
################################################################################

# Collapse the interaction-specific "One network per publication" labels back
# to a common publication category for comparisons among publication groups
norm_pub <- function(x) {
  
  x_chr <- as.character(
    x
  )
  
  ifelse(
    str_detect(
      x_chr,
      regex(
        "^One network per\\s*publication",
        ignore_case = TRUE
      )
    ),
    "One network per publication",
    x_chr
  )
}

results_corrected <- results_df %>%
  mutate(
    
    # Vázquez null-model means
    H2_null_mean =
      map_dbl(
        H2_null_values,
        ~ mean(
          .x,
          na.rm = TRUE
        )
      ),
    
    NODF_null_mean =
      map_dbl(
        weighted_NODF_null_values,
        ~ mean(
          .x,
          na.rm = TRUE
        )
      ),
    
    Mod_null_mean =
      map_dbl(
        DIRTMod_null_values,
        ~ mean(
          .x,
          na.rm = TRUE
        )
      ),
    
    # Vázquez null-model SDs
    H2_null_sd =
      map_dbl(
        H2_null_values,
        ~ sd(
          .x,
          na.rm = TRUE
        )
      ),
    
    NODF_null_sd =
      map_dbl(
        weighted_NODF_null_values,
        ~ sd(
          .x,
          na.rm = TRUE
        )
      ),
    
    Mod_null_sd =
      map_dbl(
        DIRTMod_null_values,
        ~ sd(
          .x,
          na.rm = TRUE
        )
      ),
    
    # Delta Vázquez indices
    H2_delta =
      H2 -
      H2_null_mean,
    
    NODF_delta =
      weighted_NODF -
      NODF_null_mean,
    
    Modularity_delta =
      DIRTMod -
      Mod_null_mean,
    
    # z Vázquez indices
    H2_z =
      if_else(
        !is.na(
          H2_null_sd
        ) &
          H2_null_sd > 0,
        H2_delta /
          H2_null_sd,
        NA_real_
      ),
    
    NODF_z =
      if_else(
        !is.na(
          NODF_null_sd
        ) &
          NODF_null_sd > 0,
        NODF_delta /
          NODF_null_sd,
        NA_real_
      ),
    
    Modularity_z =
      if_else(
        !is.na(
          Mod_null_sd
        ) &
          Mod_null_sd > 0,
        Modularity_delta /
          Mod_null_sd,
        NA_real_
      ),
    
    Publication_clean =
      norm_pub(
        Publication
      )
  )

# Convert empirical, Delta Vázquez, and z Vázquez indices to a common
# long format
bias_long <- results_corrected %>%
  select(
    ID,
    Publication_clean,
    TYPE,
    H2,
    weighted_NODF,
    DIRTMod,
    H2_delta,
    NODF_delta,
    Modularity_delta,
    H2_z,
    NODF_z,
    Modularity_z
  ) %>%
  
  pivot_longer(
    cols = c(
      H2,
      weighted_NODF,
      DIRTMod,
      H2_delta,
      NODF_delta,
      Modularity_delta,
      H2_z,
      NODF_z,
      Modularity_z
    ),
    names_to =
      "Metric",
    values_to =
      "Value"
  ) %>%
  
  mutate(
    ValueType = case_when(
      
      str_detect(
        Metric,
        "_delta$"
      ) ~
        "Delta",
      
      str_detect(
        Metric,
        "_z$"
      ) ~
        "Z",
      
      TRUE ~
        "Raw"
    ),
    
    Metric = case_when(
      
      Metric %in% c(
        "H2",
        "H2_delta",
        "H2_z"
      ) ~
        "Specialization (H2)",
      
      Metric %in% c(
        "weighted_NODF",
        "NODF_delta",
        "NODF_z"
      ) ~
        "Nestedness (weighted NODF)",
      
      Metric %in% c(
        "DIRTMod",
        "Modularity_delta",
        "Modularity_z"
      ) ~
        "Modularity (DIRTMod)"
    )
  )

################################################################################
#
# Calculate variation within publication groups
#
# Calculate the SD of each empirical, Delta Vázquez, and z Vázquez index
# within each publication grouping. These values are used to compare the
# "One network per publication" and "Multiple networks per publication" groups.
#
################################################################################

bias_stats <- bias_long %>%
  group_by(
    TYPE,
    Publication =
      Publication_clean,
    Metric,
    ValueType
  ) %>%
  
  summarise(
    sd_val =
      sd(
        Value,
        na.rm = TRUE
      ),
    
    n =
      n(),
    
    .groups =
      "drop"
  )

# Variation among networks in the "One network per publication" group
summary_one_network <- bias_stats %>%
  filter(
    Publication ==
      "One network per publication"
  ) %>%
  
  select(
    TYPE,
    Metric,
    ValueType,
    sd_val,
    n
  ) %>%
  
  rename(
    SD_One_Network_Per_Publication =
      sd_val,
    
    n_One_Network_Per_Publication =
      n
  )

# Within-publication variation for the "Multiple networks per publication"
# groups. Publication-specific SDs are averaged without weighting by the
# number of networks in each publication.
summary_multiple_networks <- bias_stats %>%
  filter(
    Publication !=
      "One network per publication"
  ) %>%
  
  group_by(
    TYPE,
    Metric,
    ValueType
  ) %>%
  
  summarise(
    SD_Multiple_Networks_Per_Publication =
      mean(
        sd_val,
        na.rm = TRUE
      ),
    
    n_Multiple_Networks_Publications =
      n_distinct(
        Publication
      ),
    
    n_Multiple_Networks =
      sum(
        n,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  )

# Combine variation estimates for the "One network per publication" and
# "Multiple networks per publication" groups
bias_sd_comparison <- summary_multiple_networks %>%
  left_join(
    summary_one_network,
    by = c(
      "TYPE",
      "Metric",
      "ValueType"
    )
  ) %>%
  
  mutate(
    across(
      where(
        is.numeric
      ),
      ~ round(
        .x,
        2
      )
    )
  )

################################################################################
#
# Table 5:
# Variation in Delta Vázquez network indices
#
# Compare SDs of Delta Vázquez indices between the "One network per
# publication" and "Multiple networks per publication" groups. Values for the
# "Multiple networks per publication" group are unweighted averages of
# publication-specific SDs.
#
################################################################################

table5 <- bias_sd_comparison %>%
  filter(
    ValueType ==
      "Delta"
  )

table5

# Save table
write_csv(table5, file.path(output_dir, "table_5.csv"))

################################################################################
#
# Table S10:
# Variation in z Vázquez network indices
#
# Compare SDs of z Vázquez indices between the "One network per publication"
# and "Multiple networks per publication" groups. Values for the "Multiple
# networks per publication" group are unweighted averages of publication-
# specific SDs.
#
################################################################################

table_s10 <- bias_sd_comparison %>%
  filter(
    ValueType ==
      "Z"
  )

table_s10

# Save table
write_csv(table_s10, file.path(output_dir, "table_s10.csv"))

################################################################################
#
# Publication-level summaries for Tables S9 and S11
#
# Calculate the mean, SD, and number of networks for each quantitative network
# index within each publication grouping. Empirical and Delta Vázquez values
# are reported in Table S9, while empirical and z Vázquez values are reported
# in Table S11.
#
################################################################################

pub_metric_summary <- bias_long %>%
  group_by(
    TYPE,
    Publication =
      Publication_clean,
    Metric,
    ValueType
  ) %>%
  
  summarise(
    mean =
      mean(
        Value,
        na.rm = TRUE
      ),
    
    sd =
      sd(
        Value,
        na.rm = TRUE
      ),
    
    n_networks =
      n(),
    
    .groups =
      "drop"
  )

metric_key <- c(
  "Specialization (H2)" =
    "Specialization",
  "Nestedness (weighted NODF)" =
    "Nestedness",
  "Modularity (DIRTMod)" =
    "Modularity"
)

################################################################################
#
# Table S9:
# Publication-level empirical and Delta Vázquez network indices
#
# Report publication-level means and SDs of the empirical quantitative network
# indices and their Delta Vázquez null-corrected values for specialization,
# nestedness, and modularity. Unweighted averages are also calculated across
# the seven pollination and 18 seed-dispersal publications in the "Multiple
# networks per publication" group.
#
################################################################################

table_s9 <- pub_metric_summary %>%
  filter(
    ValueType %in% c(
      "Raw",
      "Delta"
    )
  ) %>%
  
  mutate(
    MetricShort =
      recode(
        Metric,
        !!!metric_key
      )
  ) %>%
  
  select(
    TYPE,
    Publication,
    ValueType,
    MetricShort,
    mean,
    sd,
    n_networks
  ) %>%
  
  pivot_wider(
    names_from = c(
      ValueType,
      MetricShort
    ),
    values_from = c(
      mean,
      sd
    ),
    names_glue =
      "{ValueType}_{MetricShort}_{.value}"
  ) %>%
  
  select(
    TYPE,
    Publication,
    n_networks,
    
    Raw_Specialization_mean,
    Raw_Specialization_sd,
    Delta_Specialization_mean,
    Delta_Specialization_sd,
    
    Raw_Nestedness_mean,
    Raw_Nestedness_sd,
    Delta_Nestedness_mean,
    Delta_Nestedness_sd,
    
    Raw_Modularity_mean,
    Raw_Modularity_sd,
    Delta_Modularity_mean,
    Delta_Modularity_sd
  )

# Unweighted average across publications in the
# "Multiple networks per publication" group
table_s9_average <- table_s9 %>%
  filter(
    Publication !=
      "One network per publication"
  ) %>%
  
  group_by(
    TYPE
  ) %>%
  
  summarise(
    across(
      where(is.numeric) &
        !all_of("n_networks"),
      ~ mean(
        .x,
        na.rm = TRUE
      )
    ),
    
    .groups =
      "drop"
  ) %>%
  
  mutate(
    Publication = case_when(
      
      TYPE ==
        "Pollination" ~
        "Average unweighted 7 PL publications above",
      
      TYPE ==
        "Seed Dispersal" ~
        "Average unweighted 18 SD publications above"
    ),
    
    n_networks =
      NA_real_
  ) %>%
  
  select(
    names(
      table_s9
    )
  )

# Add unweighted averages to Table S9
table_s9 <- bind_rows(
  table_s9,
  table_s9_average
) %>%
  
  mutate(
    row_order = case_when(
      
      TYPE ==
        "Pollination" &
        Publication !=
        "One network per publication" &
        !str_detect(
          Publication,
          "^Average unweighted"
        ) ~
        1,
      
      TYPE ==
        "Pollination" &
        str_detect(
          Publication,
          "^Average unweighted"
        ) ~
        2,
      
      TYPE ==
        "Seed Dispersal" &
        Publication !=
        "One network per publication" &
        !str_detect(
          Publication,
          "^Average unweighted"
        ) ~
        3,
      
      TYPE ==
        "Seed Dispersal" &
        str_detect(
          Publication,
          "^Average unweighted"
        ) ~
        4,
      
      TYPE ==
        "Pollination" &
        Publication ==
        "One network per publication" ~
        5,
      
      TYPE ==
        "Seed Dispersal" &
        Publication ==
        "One network per publication" ~
        6
    )
  ) %>%
  
  arrange(
    row_order,
    Publication
  ) %>%
  
  select(
    -row_order
  ) %>%
  
  mutate(
    across(
      where(
        is.numeric
      ),
      ~ round(
        .x,
        2
      )
    )
  )

table_s9

# Save table
write_csv(table_s9, file.path(output_dir, "table_s9.csv"))

################################################################################
#
# Table S11:
# Publication-level empirical and z Vázquez network indices
#
# Report publication-level means and SDs of the empirical quantitative network
# indices and their z Vázquez null-corrected values for specialization,
# nestedness, and modularity. Unweighted averages are also calculated across
# the seven pollination and 18 seed-dispersal publications in the "Multiple
# networks per publication" group.
#
################################################################################

table_s11 <- pub_metric_summary %>%
  filter(
    ValueType %in% c(
      "Raw",
      "Z"
    )
  ) %>%
  
  mutate(
    MetricShort =
      recode(
        Metric,
        !!!metric_key
      )
  ) %>%
  
  select(
    TYPE,
    Publication,
    ValueType,
    MetricShort,
    mean,
    sd,
    n_networks
  ) %>%
  
  pivot_wider(
    names_from = c(
      ValueType,
      MetricShort
    ),
    values_from = c(
      mean,
      sd
    ),
    names_glue =
      "{ValueType}_{MetricShort}_{.value}"
  ) %>%
  
  select(
    TYPE,
    Publication,
    n_networks,
    
    Raw_Specialization_mean,
    Raw_Specialization_sd,
    Z_Specialization_mean,
    Z_Specialization_sd,
    
    Raw_Nestedness_mean,
    Raw_Nestedness_sd,
    Z_Nestedness_mean,
    Z_Nestedness_sd,
    
    Raw_Modularity_mean,
    Raw_Modularity_sd,
    Z_Modularity_mean,
    Z_Modularity_sd
  )

# Unweighted average across publications in the
# "Multiple networks per publication" group
table_s11_average <- table_s11 %>%
  filter(
    Publication !=
      "One network per publication"
  ) %>%
  
  group_by(
    TYPE
  ) %>%
  
  summarise(
    across(
      where(is.numeric) &
        !all_of("n_networks"),
      ~ mean(
        .x,
        na.rm = TRUE
      )
    ),
    
    .groups =
      "drop"
  ) %>%
  
  mutate(
    Publication = case_when(
      
      TYPE ==
        "Pollination" ~
        "Average unweighted 7 PL publications above",
      
      TYPE ==
        "Seed Dispersal" ~
        "Average unweighted 18 SD publications above"
    ),
    
    n_networks =
      NA_real_
  ) %>%
  
  select(
    names(
      table_s11
    )
  )

# Add unweighted averages to Table S11
table_s11 <- bind_rows(
  table_s11,
  table_s11_average
) %>%
  
  mutate(
    row_order = case_when(
      
      TYPE ==
        "Pollination" &
        Publication !=
        "One network per publication" &
        !str_detect(
          Publication,
          "^Average unweighted"
        ) ~
        1,
      
      TYPE ==
        "Pollination" &
        str_detect(
          Publication,
          "^Average unweighted"
        ) ~
        2,
      
      TYPE ==
        "Seed Dispersal" &
        Publication !=
        "One network per publication" &
        !str_detect(
          Publication,
          "^Average unweighted"
        ) ~
        3,
      
      TYPE ==
        "Seed Dispersal" &
        str_detect(
          Publication,
          "^Average unweighted"
        ) ~
        4,
      
      TYPE ==
        "Pollination" &
        Publication ==
        "One network per publication" ~
        5,
      
      TYPE ==
        "Seed Dispersal" &
        Publication ==
        "One network per publication" ~
        6
    )
  ) %>%
  
  arrange(
    row_order,
    Publication
  ) %>%
  
  select(
    -row_order
  ) %>%
  
  mutate(
    across(
      where(
        is.numeric
      ),
      ~ round(
        .x,
        2
      )
    )
  )

table_s11

# Save table
write_csv(table_s11, file.path(output_dir, "table_s11.csv"))
