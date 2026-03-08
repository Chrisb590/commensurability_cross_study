library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(purrr)
library(readr)
library(patchwork)
library(showtext)
library(here)
showtext_auto()
font_add("courier", "C:/Windows/Fonts/cour.ttf")

################################################################################
#
# Data
#
################################################################################

# Vaznull randomization results
results_df <- read_csv(here("all_results_vaznull_randomizations.csv"))

# Metadata
metadata <- read_csv(here("general_network_information.csv"))
metadata$Publication <- as.character(metadata$Publication)

# Parse Vaznull columns into lists
results_df <- results_df %>%
  mutate(
    H2_null_values = str_split(H2_null_values, ";") %>% map(~ as.numeric(.x)),
    weighted_NODF_null_values = str_split(weighted_NODF_null_values, ";") %>% map(~ as.numeric(.x)),
    DIRTMod_null_values = str_split(DIRTMod_null_values, ";") %>% map(~ as.numeric(.x))
  )

# Join Vaznull results with metadata
results_df <- results_df %>%
  inner_join(metadata, by = "ID") 

# Add network type when "One network per publication"
results_df <- results_df %>%
  mutate(
    Publication = case_when(
      Publication == "One network per publication" & TYPE == "Pollination" ~ "One network per\npublication (Pollination)",
      Publication == "One network per publication" & TYPE == "Seed Dispersal" ~ "One network per\npublication (Seed Dispersal)",
      TRUE ~ Publication
    )
  )

# Order publications so one-network types appear first
ordered_publications <- c(
  "One network per\npublication (Pollination)",
  "One network per\npublication (Seed Dispersal)",
  sort(setdiff(unique(results_df$Publication), c(
    "One network per\npublication (Pollination)",
    "One network per\npublication (Seed Dispersal)"
  )))
)
results_df$Publication <- factor(results_df$Publication, levels = ordered_publications)

# Results in long format
empirical_long <- results_df %>%
  select(ID, Publication, H2, weighted_NODF, DIRTMod) %>%
  pivot_longer(cols = c(H2, weighted_NODF, DIRTMod), names_to = "metric", values_to = "value") %>%
  mutate(
    Type = "Empirical",
    Metric = case_when(
      metric == "H2" ~ "Specialization (H2)",
      metric == "weighted_NODF" ~ "Nestedness (weighted NODF)",
      metric == "DIRTMod" ~ "Modularity (DIRTMod)"
    )
  ) %>%
  select(ID, Publication, value, Type, Metric)

# Long format for all Vaznull values (1000 rows per index)
null_long <- results_df %>%
  select(ID, Publication, H2_null_values, weighted_NODF_null_values, DIRTMod_null_values) %>%
  pivot_longer(cols = c(H2_null_values, weighted_NODF_null_values, DIRTMod_null_values),
               names_to = "metric", values_to = "value_list") %>%
  unnest(value_list) %>%
  mutate(
    value = value_list,
    Type = "Vázquez null",
    Metric = case_when(
      metric == "H2_null_values" ~ "Specialization (H2)",
      metric == "weighted_NODF_null_values" ~ "Nestedness (weighted NODF)",
      metric == "DIRTMod_null_values" ~ "Modularity (DIRTMod)"
    )
  ) %>%
  select(ID, Publication, value, Type, Metric)

# Both long formats of raw data and Vaznull corrected indices
df_long <- bind_rows(empirical_long, null_long) %>%
  filter(!is.na(Publication))

df_specialization <- df_long %>% filter(Metric == "Specialization (H2)")
df_nestedness     <- df_long %>% filter(Metric == "Nestedness (weighted NODF)")
df_modularity     <- df_long %>% filter(Metric == "Modularity (DIRTMod)")

# Colours for the different publication groupings
publication_types <- results_df %>%
  distinct(Publication, TYPE) %>%
  mutate(
    label_color = ifelse(TYPE == "Pollination", "#e6550d", "#41ab5d")
  )

df_specialization <- df_specialization %>%
  left_join(publication_types, by = "Publication")
df_nestedness <- df_nestedness %>%
  left_join(publication_types, by = "Publication")
df_modularity <- df_modularity %>%
  left_join(publication_types, by = "Publication")

################################################################################
#
# Figure S5: Boxplot of raw vs Vaznull model for each publication grouping and 
# the three indices of specialization, nestedness, and modularity
#
################################################################################

fill_scale <- scale_fill_manual(
  values = c("Empirical" = "#fcae91", "Vázquez null" = "#9972af"),
  labels = c(
    "Empirical" = "Empirical",
    "Vázquez null" = "<span style='font-family:courier;'>vaznull</span>"
  )
)

###########################
# Top panel: Specialization
###########################

p1 <- ggplot(df_specialization, aes(x = Publication, y = value, fill = Type)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.8, width = 0.6, color = "grey30", size = 0.3) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black", size = 0.5) +
  ylim(0, 1) +
  labs(y = "Specialization (H2')", title = NULL, x = NULL, fill = NULL) +
  fill_scale +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "grey50", fill = NA, size = 0.8)
  )

##########################
# Middle panel: Nestedness
##########################

p2 <- ggplot(df_nestedness, aes(x = Publication, y = value, fill = Type)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.8, width = 0.6, color = "grey30", size = 0.3) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black", size = 0.5) +
  ylim(0, 100) +
  labs(y = "Nestedness (wNODF)", title = NULL, x = NULL, fill = NULL) +
  fill_scale +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "grey50", fill = NA, size = 0.8)
  )

##########################
# Bottom panel: Modularity
##########################

p3 <- ggplot(df_modularity, aes(x = Publication, y = value, fill = Type)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.8, width = 0.6, color = "grey30", size = 0.3) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "black", size = 0.5) +
  ylim(0, 1) +
  labs(y = "Modularity (DIRTLPAwb+)", title = NULL, x = NULL, fill = NULL) +
  fill_scale +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = df_modularity %>%
                                 distinct(Publication, label_color) %>%
                                 arrange(factor(Publication, levels = levels(df_modularity$Publication))) %>%
                                 pull(label_color)),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "grey50", fill = NA, size = 0.8)
  )

# Combine plots
p_combined <- p1 / p2 / p3 +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "top",
    legend.text = ggtext::element_markdown()
  )

print(p_combined)

ggsave(
  "combined_plot_vaznull.pdf",
  plot = p_combined,
  width = 13,
  height = 11
)

################################################################################
#
# Table 4 and S8: Delta and Z indices of specialization, nestedness, and 
# modularity
#
################################################################################

# Normalize publication label when a variant of one network per publication
norm_pub <- function(x) {
  x_chr <- as.character(x)
  ifelse(
    str_detect(x_chr, regex("^One network per\\s*publication", ignore_case = TRUE)),
    "One network per publication",
    x_chr
  )
}

null_mean <- function(x) mean(x, na.rm = TRUE)
null_sd   <- function(x) sd(x, na.rm = TRUE)

# Calculating Delta and Z indices
results_corrected <- results_df %>%
  mutate(
    # Null means
    H2_null_mean   = map_dbl(H2_null_values, null_mean),
    NODF_null_mean = map_dbl(weighted_NODF_null_values, null_mean),
    Mod_null_mean  = map_dbl(DIRTMod_null_values, null_mean),
    
    # Null SDs
    H2_null_sd   = map_dbl(H2_null_values, null_sd),
    NODF_null_sd = map_dbl(weighted_NODF_null_values, null_sd),
    Mod_null_sd  = map_dbl(DIRTMod_null_values, null_sd),
    
    # Deltas (already what you do)
    H2_delta         = H2 - H2_null_mean,
    NODF_delta       = weighted_NODF - NODF_null_mean,
    Modularity_delta = DIRTMod - Mod_null_mean,
    
    # Z = delta / sd(null)
    H2_z = ifelse(!is.na(H2_null_sd)   & H2_null_sd   > 0, H2_delta         / H2_null_sd,   NA_real_),
    NODF_z = ifelse(!is.na(NODF_null_sd) & NODF_null_sd > 0, NODF_delta       / NODF_null_sd, NA_real_),
    Modularity_z = ifelse(!is.na(Mod_null_sd)  & Mod_null_sd  > 0, Modularity_delta / Mod_null_sd,  NA_real_),
    
    Publication_clean = norm_pub(Publication)
  )

# Long-format of results_corrected
bias_long <- results_corrected %>%
  select(ID, Publication_clean, TYPE,
         H2, weighted_NODF, DIRTMod,
         H2_delta, NODF_delta, Modularity_delta,
         H2_z, NODF_z, Modularity_z) %>%
  pivot_longer(
    cols      = c(H2, weighted_NODF, DIRTMod,
                  H2_delta, NODF_delta, Modularity_delta,
                  H2_z, NODF_z, Modularity_z),
    names_to  = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    ValueType = case_when(
      str_detect(Metric, "_delta$") ~ "Delta",
      str_detect(Metric, "_z$")     ~ "Z",
      TRUE                          ~ "Raw"
    ),
    Metric = case_when(
      Metric %in% c("H2", "H2_delta", "H2_z")                         ~ "Specialization (H2)",
      Metric %in% c("weighted_NODF", "NODF_delta", "NODF_z")          ~ "Nestedness (weighted NODF)",
      Metric %in% c("DIRTMod", "Modularity_delta", "Modularity_z")    ~ "Modularity (DIRTMod)"
    )
  )

# Publication-level standard deviation summaries for the raw and Delta indices
bias_stats <- bias_long %>%
  group_by(TYPE, Publication = Publication_clean, Metric, ValueType) %>%
  summarise(
    sd_val   = sd(Value,   na.rm = TRUE),
    n        = n(),
    .groups  = "drop"
  )

# Standard deviation for all one network per publication by network type and 
# index (both raw and Delta)
summary_one_net <- bias_stats %>%
  filter(Publication == "One network per publication") %>%
  select(TYPE, Metric, ValueType, sd_val, n) %>%
  rename(
    SD_One_Network_Per_Publication = sd_val,
    n_One_Network_Per_Publication  = n
  )

# Standard deviation for publications that produced multiple networks by network 
# type and index (both raw and Delta). Not weighted by number of networks per 
# publication
summary_other_unweighted <- bias_stats %>%
  filter(Publication != "One network per publication") %>%
  group_by(TYPE, Metric, ValueType) %>%
  summarise(
    Avg_SD_Other = mean(sd_val, na.rm = TRUE),
    n_publications_other = n_distinct(Publication),
    n_networks_other     = sum(n, na.rm = TRUE),
    .groups = "drop"
  )

# Merging SD for all publication groupings
bias_sd_comparison <- summary_other_unweighted %>%
  left_join(summary_one_net, by = c("TYPE", "Metric", "ValueType")) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

print(bias_sd_comparison)

################################################################################
#
# Table S7 and S9: Vaznull corrected indices of specialization, nestedness, and 
# modularity
#
################################################################################

# Summarize each index (raw and Delta) at the publication level:
# mean, SD, and number of networks per TYPE × publication × index × value type
pub_metric_summary <- bias_long %>%
  filter(ValueType %in% c("Raw", "Delta", "Z")) %>%
  group_by(TYPE, Publication = Publication_clean, Metric, ValueType) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    sd   = sd(Value,   na.rm = TRUE),
    n_networks = n(),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), round, 2))

metric_key <- c(
  "Specialization (H2)"        = "Specialization",
  "Nestedness (weighted NODF)" = "Nestedness",
  "Modularity (DIRTMod)"       = "Modularity"
)

pub_metric_summary_12 <- pub_metric_summary %>%
  mutate(
    MetricShort = recode(Metric, !!!metric_key),
    vt = ValueType  # Raw / Delta / Z
  ) %>%
  select(TYPE, Publication, vt, MetricShort, mean, sd, n_networks) %>%
  tidyr::pivot_wider(
    names_from  = c(vt, MetricShort),
    values_from = c(mean, sd),
    names_glue  = "{vt}_{MetricShort}_{.value}"
  ) %>%
  select(
    TYPE, Publication, n_networks,
    
    Raw_Specialization_mean,   Raw_Specialization_sd,
    Raw_Nestedness_mean,       Raw_Nestedness_sd,
    Raw_Modularity_mean,       Raw_Modularity_sd,
    
    Delta_Specialization_mean, Delta_Specialization_sd,
    Delta_Nestedness_mean,     Delta_Nestedness_sd,
    Delta_Modularity_mean,     Delta_Modularity_sd,
    
    Z_Specialization_mean,     Z_Specialization_sd,
    Z_Nestedness_mean,         Z_Nestedness_sd,
    Z_Modularity_mean,         Z_Modularity_sd
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

pub_metric_summary_12

################################################################################
#
# Table S10: % difference between each empirical network and average of 1000 
# Patefield null network models for specialization, nestedness, and modularity
#
################################################################################

percent_difference_for_each_network <- results_df %>%
  mutate(
    h2_pct = 100 * (map_dbl(H2_null_values, mean) - H2) / H2,
    wNODF_pct = 100 * (map_dbl(weighted_NODF_null_values, mean) - weighted_NODF) / weighted_NODF,
    DIRTLP_pct = 100 * (map_dbl(DIRTMod_null_values, mean) - DIRTMod) / DIRTMod
  ) %>%
  select(
    ID,
    TYPE,
    `h2' (%)` = h2_pct,
    `wNODF (%)` = wNODF_pct,
    `DIRTLP (%)` = DIRTLP_pct
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

percent_difference_for_each_network
