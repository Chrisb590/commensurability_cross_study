library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(here)
library(readr)

################################################################################
#
# Data
#
################################################################################

# Location of data
metadata <- read_csv(here("general_network_information.csv"))
fricke_metadata <- read_csv(here("fricke_metadata.csv"))
pollination_metadata <- read_csv(here("pollination_sampling_metadata.csv"))

# Making sure Fricke networks are labelled the same as pollination
fricke_metadata$net_id <- gsub(" ", "_", fricke_metadata$net_id)
fricke_metadata$net_id <- gsub("-", "_", fricke_metadata$net_id)
fricke_metadata$net_id <- gsub("&", "and", fricke_metadata$net_id)
fricke_metadata$net_id <- gsub("\\.", "", fricke_metadata$net_id)

# Combining metadata
both_combined <- left_join(metadata,fricke_metadata,by=c("ID"="net_id"))
both_combined <- both_combined %>%
  rows_update(
    pollination_metadata %>% select(any_of(names(both_combined))),
    by = "ID"
  )

################################################################################
#
# Fixing strings for sampling_method and sampling_unit
#
################################################################################

# Normalize spacing
both_combined$sampling_method <- str_squish(both_combined$sampling_method)
both_combined$sampling_unit   <- str_squish(both_combined$sampling_unit)

# Simple spelling fixes (sampling_method)
method_recode <- c(
  "mistnetting"         = "mist netting",
  "mitnetting"          = "mist netting",
  "general observations"= "general obs",
  "focal observations"  = "focal obs",
  "plot observation"    = "plot obs"
)
for (pat in names(method_recode)) {
  both_combined$sampling_method <- gsub(
    pat, method_recode[[pat]],
    both_combined$sampling_method,
    ignore.case = TRUE
  )
}

# Simple spelling fixes (sampling_unit)
unit_recode <- c(
  "fecals with species"                 = "fecal with species",
  "bin fecal with species"              = "fecal with species",
  "feeding visit, fecal with feces"     = "feeding visits, fecal with species"
)
for (pat in names(unit_recode)) {
  both_combined$sampling_unit <- gsub(
    pat, unit_recode[[pat]],
    both_combined$sampling_unit,
    ignore.case = TRUE
  )
}

# If sampling_unit is exactly "visits", relabel to "frugivore visits"
both_combined$sampling_unit <- ifelse(
  str_detect(both_combined$sampling_unit, regex("^visits$", ignore_case = TRUE)),
  "frugivore visits",
  both_combined$sampling_unit
)

# Rename "frugivore visits" to "feeding visits"
both_combined$sampling_unit <- gsub(
  "^frugivore\\s+visits$",
  "feeding visits",
  both_combined$sampling_unit,
  ignore.case = TRUE
)

# Seed-dispersal specific rule: focal obs -> Focal obs (flora)
rows_seed <- grepl("(?i)^\\s*seed\\s+dispersal\\s*$", both_combined$TYPE, perl = TRUE)
both_combined$sampling_method[rows_seed] <- gsub(
  "(?i)\\bfocal\\s*obs\\b(?!\\s*\\(\\s*flora\\s*\\))",
  "Focal obs (flora)",
  both_combined$sampling_method[rows_seed],
  perl = TRUE
)

# Also catch the case where the whole string is exactly "focal obs"
both_combined$sampling_method <- ifelse(
  str_detect(both_combined$sampling_method, regex("^focal\\s*obs$", ignore_case = TRUE)),
  "Focal obs (flora)",
  both_combined$sampling_method
)

# Replace ", transect walk" / ", transect-walk" -> ", transect"
both_combined$sampling_method <- str_replace_all(
  both_combined$sampling_method,
  regex(",\\s*transect\\s*-?\\s*walk\\b", ignore_case = TRUE),
  ", transect"
)

# Ensure bare "focal" not followed by "obs" becomes "focal obs"
both_combined$sampling_method <- gsub(
  "\\bfocal\\b(?!\\s*obs\\b)",
  "focal obs",
  both_combined$sampling_method,
  ignore.case = TRUE,
  perl = TRUE
)

# Convert "targeted (fauna|flora)" -> "targeted fauna|flora"
both_combined$sampling_method <- gsub(
  "\\btargeted\\s*\\(\\s*(fauna|flora)\\s*\\)",
  "targeted \\1",
  both_combined$sampling_method,
  ignore.case = TRUE
)

# Drop ", incidental" after Focal obs (flora)
both_combined$sampling_method <- gsub(
  "(?i)\\bfocal\\s*obs\\s*\\(\\s*flora\\s*\\)\\s*,\\s*incidental\\b",
  "Focal obs (flora)",
  both_combined$sampling_method,
  perl = TRUE
)

# Final whitespace tidy (once)
both_combined$sampling_method <- str_squish(both_combined$sampling_method)
both_combined$sampling_unit   <- str_squish(both_combined$sampling_unit)

# Capitalizing the first letter
cap_first <- function(x) sub("^\\s*([[:alpha:]])", "\\U\\1", x, perl = TRUE)
both_combined$sampling_method <- cap_first(both_combined$sampling_method)
both_combined$sampling_unit   <- cap_first(both_combined$sampling_unit)

# Standardize network type to exactly two categories (Pollination vs.
# Seed-dispersal)
both2 <- both_combined %>%
  mutate(
    TYPE2 = ifelse(TYPE == "Pollination", "Pollination", "Seed Dispersal")
  )

################################################################################
#
# Publication collapsed tables for percentage of sample method and edge weight
#
################################################################################

# Explicit ordering of network types
wanted_types <- c("Pollination", "Seed Dispersal")

# Create publication-collapsed sampling-method data (one row per publication 
# × method × network type)
method_collapsed_long <- both2 %>%
  filter(!is.na(sampling_method), str_squish(sampling_method) != "") %>%
  mutate(
    sampling_method = str_squish(sampling_method),
    pub_key = if_else(
      grepl("^\\s*one\\s+network\\s+per\\s+publication\\s*$",
            Publication, ignore.case = TRUE),
      paste0("SINGLE:", ID),
      Publication
    )
  ) %>%
  distinct(TYPE2, pub_key, sampling_method)

# Summarize sampling methods at the publication-collapsed level
pct_collapsed_method <- method_collapsed_long %>%
  count(TYPE2, sampling_method, name = "n_pub_method") %>%
  group_by(TYPE2) %>%
  mutate(
    total_pub = sum(n_pub_method),
    pct = 100 * n_pub_method / total_pub
  ) %>%
  ungroup() %>%
  arrange(TYPE2, desc(pct))

# Sampling-method summary for pollination networks
pct_collapsed_method_pollination <- pct_collapsed_method %>%
  filter(TYPE2 == "Pollination") %>%
  select(-TYPE2)

# Sampling-method summary for seed-dispersal networks
pct_collapsed_method_seed <- pct_collapsed_method %>%
  filter(TYPE2 == "Seed Dispersal") %>%
  select(-TYPE2)

# Create publication-collapsed edge weight data (one row per publication × 
# method × network type)
unit_collapsed_long <- both2 %>%
  filter(!is.na(sampling_unit), str_squish(sampling_unit) != "") %>%
  mutate(
    sampling_unit = str_squish(sampling_unit),
    pub_key = if_else(
      grepl("^\\s*one\\s+network\\s+per\\s+publication\\s*$",
            Publication, ignore.case = TRUE),
      paste0("SINGLE:", ID),
      Publication
    )
  ) %>%
  distinct(TYPE2, pub_key, sampling_unit)

# Summarize edge weights at the publication-collapsed level
pct_collapsed_unit <- unit_collapsed_long %>%
  count(TYPE2, sampling_unit, name = "n_pub_unit") %>%
  group_by(TYPE2) %>%
  mutate(
    total_pub = sum(n_pub_unit),
    pct = 100 * n_pub_unit / total_pub
  ) %>%
  ungroup() %>%
  arrange(TYPE2, desc(pct))

# Edge weight summary for pollination networks
pct_collapsed_unit_pollination <- pct_collapsed_unit %>%
  filter(TYPE2 == "Pollination") %>%
  select(-TYPE2)

# Edge weight summary for seed-dispersal networks
pct_collapsed_unit_seed <- pct_collapsed_unit %>%
  filter(TYPE2 == "Seed Dispersal") %>%
  select(-TYPE2)

# Print tables
print(pct_collapsed_method_pollination, n = Inf)
print(pct_collapsed_method_seed, n = Inf)
print(pct_collapsed_unit_pollination, n = Inf)
print(pct_collapsed_unit_seed, n = Inf)

################################################################################
#
# Figure 4
#
################################################################################

fill_vals  <- c("Pollination" = "#fdbb84", "Seed Dispersal" = "#a1dab4")
color_vals <- c("Pollination" = "#d95f02", "Seed Dispersal" = "#1b9e77")

# Prepare sampling-method figure
plot_df_method <- both2 %>%
  filter(!is.na(sampling_method), str_squish(sampling_method) != "") %>%
  mutate(
    sampling_method = str_squish(sampling_method),
    pub_key = if_else(
      grepl("^\\s*one\\s+network\\s+per\\s+publication\\s*$",
            Publication, ignore.case = TRUE),
      paste0("SINGLE:", ID),
      Publication
    ),
    TYPE2 = as.character(TYPE2)
  ) %>%
  distinct(pub_key, TYPE2, sampling_method) %>%
  mutate(method_label = str_wrap(sampling_method, width = 30))

# Prepare edge weight figure
plot_df_unit <- both2 %>%
  filter(!is.na(sampling_unit), str_squish(sampling_unit) != "") %>%
  mutate(
    sampling_unit = str_squish(sampling_unit),
    pub_key = if_else(
      grepl("^\\s*one\\s+network\\s+per\\s+publication\\s*$",
            Publication, ignore.case = TRUE),
      paste0("SINGLE:", ID),
      Publication
    ),
    TYPE2 = as.character(TYPE2)
  ) %>%
  distinct(pub_key, TYPE2, sampling_unit) %>%
  mutate(unit_label = str_wrap(sampling_unit, width = 30))

# Ordering by publication collapsed frequency 
method_order <- plot_df_method %>% count(method_label, name = "n") %>% arrange(desc(n)) %>% pull(method_label)
unit_order   <- plot_df_unit   %>% count(unit_label,   name = "n") %>% arrange(desc(n)) %>% pull(unit_label)

# Plot 1: Sampling Method
p1 <- ggplot(
  plot_df_method,
  aes(
    x = factor(method_label, levels = method_order),
    fill  = factor(TYPE2, levels = wanted_types),
    color = factor(TYPE2, levels = wanted_types)
  )
) +
  geom_bar(width = 0.8, size = 0.6) +
  scale_fill_manual(
    values = fill_vals,
    limits = wanted_types,
    drop = FALSE,
    na.value = "grey80",
    labels = c("Pollination", "Seed-dispersal"),
    name   = NULL
  ) +
  scale_color_manual(
    values = color_vals,
    limits = wanted_types,
    drop = FALSE,
    na.value = NA,
    labels = c("Pollination", "Seed-dispersal"),
    name   = NULL
  ) +
  guides(fill = guide_legend(override.aes = list(color = color_vals, fill = fill_vals))) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16, color = "black"),
    axis.text.y  = element_text(size = 20, color = "black"),
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(face = "bold", color = "black"),
    axis.title.y = element_text(face = "bold", color = "black"),
    plot.margin  = margin(10, 10, 20, 10),
    panel.border = element_rect(colour = "grey", fill = NA, linewidth = 0.1)
  ) +
  xlab("Sampling method") +
  ylab("Network count\n(Publication-collapsed)") +
  coord_cartesian(clip = "off")

# Plot 2: Edge weight
p2 <- ggplot(
  plot_df_unit,
  aes(
    x = factor(unit_label, levels = unit_order),
    fill  = factor(TYPE2, levels = wanted_types),
    color = factor(TYPE2, levels = wanted_types)
  )
) +
  geom_bar(width = 0.8, size = 0.6) +
  scale_fill_manual(
    values = fill_vals,
    limits = wanted_types,
    drop = FALSE,
    na.value = "grey80",
    labels = c("Pollination", "Seed-dispersal"),
    name   = NULL
  ) +
  scale_color_manual(
    values = color_vals,
    limits = wanted_types,
    drop = FALSE,
    na.value = NA,
    labels = c("Pollination", "Seed-dispersal"),
    name   = NULL
  ) +
  guides(fill = guide_legend(override.aes = list(color = color_vals, fill = fill_vals))) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16, color = "black"),
    axis.text.y  = element_text(size = 20, color = "black"),
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(face = "bold", color = "black"),
    axis.title.y = element_text(face = "bold", color = "black"),
    plot.margin  = margin(10, 10, 20, 10),
    panel.border = element_rect(colour = "grey", fill = NA, linewidth = 0.1)
  ) +
  xlab("Edge weight") +
  ylab("Network count\n(Publication-collapsed)") +
  coord_cartesian(clip = "off")

# Combine plots with single legend
combined_plot <- (p1 / p2) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "top",
    axis.text.y  = element_text(size = 20, color = "black"),
    legend.text  = element_text(size = 18),
    plot.title   = element_text(size = 16, face = "bold"),
    axis.text.x  = element_text(size = 15, angle = 45,
                                hjust = 1, vjust = 1, color = "black"),
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black")
  )
print(combined_plot)

ggsave(
  "sampling_unit_and_sampling_method_publication_collapsed.pdf",
  plot = combined_plot,
  device = cairo_pdf,
  width = 15,
  height = 11,
  units = "in"
)
