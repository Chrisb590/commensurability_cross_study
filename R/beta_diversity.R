library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(here)
library(readr)

################################################################################
#
# Data
#
################################################################################

metadata <- read_csv(here("general_network_information.csv"))

# Directories where beta-diversity (for rows and columns) CSVs are stored
row_dir <- here("networks", "row_names")
col_dir <- here("networks", "column_names")

# Interaction types for file names
interaction_types <- c("Seed_Dispersal", "Pollination")

# Names
display_map <- c(
  "Seed_Dispersal" = "Seed-dispersal",
  "Pollination"    = "Pollination"
)

################################################################################
#
# Getting pairwise beta-diversity information
#
################################################################################

# Function to extract pairwise beta-diversity and classify pairs by publication 
# group
get_pairwise_beta <- function(beta_file, interaction_key, data_type) {
  
  beta_mat <- as.matrix(read.csv(beta_file, row.names = 1, check.names = FALSE))
  ids <- rownames(beta_mat)
  
  pairwise_df <- data.frame(
    ID1 = character(),
    ID2 = character(),
    beta = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(length(ids) - 1)) {
    for (j in (i + 1):length(ids)) {
      pairwise_df <- rbind(
        pairwise_df,
        data.frame(
          ID1 = ids[i],
          ID2 = ids[j],
          beta = beta_mat[i, j],
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  pairwise_df %>%
    mutate(
      ID1_clean = sub("(_row_names|_column_names)$", "", ID1),
      ID2_clean = sub("(_row_names|_column_names)$", "", ID2)
    ) %>%
    left_join(metadata %>% select(ID, Publication), by = c("ID1_clean" = "ID")) %>%
    rename(Pub1 = Publication) %>%
    left_join(metadata %>% select(ID, Publication), by = c("ID2_clean" = "ID")) %>%
    rename(Pub2 = Publication) %>%
    mutate(
      Group = case_when(
        Pub1 == "One network per publication" & Pub2 == "One network per publication" ~ "One network per publication",
        Pub1 != "One network per publication" & Pub2 != "One network per publication" & Pub1 == Pub2 ~ "Within Publication",
        TRUE ~ NA_character_
      ),
      Interaction_key = interaction_key,
      DataType = data_type
    ) %>%
    filter(!is.na(Group)) %>%
    select(ID1, ID2, beta, Pub1, Pub2, Group, Interaction_key, DataType)
}

# Storage for pairwise beta-diversity
all_pairwise <- list()

# For each interaction type, read available row- and column-level beta-diversity 
# files and store pairwise results
for (intr in interaction_types) {
  row_file <- file.path(row_dir, paste0(intr, "_beta_diversity_row.csv"))
  col_file <- file.path(col_dir, paste0(intr, "_beta_diversity_column.csv"))
  
  if (file.exists(row_file)) {
    all_pairwise[[paste0(intr, "_row")]] <- get_pairwise_beta(row_file, intr, "row")
  }
  if (file.exists(col_file)) {
    all_pairwise[[paste0(intr, "_col")]] <- get_pairwise_beta(col_file, intr, "column")
  }
}

# x-axis in the boxplot
combined_df <- bind_rows(all_pairwise) %>%
  mutate(
    Interaction = factor(Interaction_key, levels = names(display_map), labels = unname(display_map)),
    
    Combo = case_when(
      Group == "One network per publication" & DataType == "row"    ~ "Row - One per pub",
      Group == "One network per publication" & DataType == "column" ~ "Column - One per pub",
      Group == "Within Publication"          & DataType == "row"    ~ "Row - Within pub",
      Group == "Within Publication"          & DataType == "column" ~ "Column - Within pub",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Combo)) %>%
  mutate(ComboInteraction = interaction(Combo, Interaction_key, drop = TRUE))

# Ensure all 'One network per publication' levels come first
levels_one <- combined_df %>%
  filter(Group == "One network per publication") %>%
  distinct(ComboInteraction) %>%
  pull()

# Multiple networks per publication levels come second
levels_within <- combined_df %>%
  filter(Group == "Within Publication") %>%
  distinct(ComboInteraction) %>%
  pull()

# Setting the publication levels
combo_levels <- c(levels_one, levels_within)
combined_df <- combined_df %>%
  mutate(ComboInteraction = factor(ComboInteraction, levels = combo_levels))

# Separator position between one-per-pub and within-pub groups
separator_position <- length(levels_one) + 0.5

# Simplified axis labels with line breaks (Plant vs Animal guild)
combo_labels <- combined_df %>%
  distinct(ComboInteraction, DataType) %>%
  mutate(ComboLabel = ifelse(DataType == "row", "Plant\nguild", "Animal\nguild")) %>%
  select(ComboInteraction, ComboLabel) %>%
  deframe()

# Counts for number of counts above each boxplot
combo_counts <- combined_df %>%
  count(ComboInteraction) %>%
  mutate(label = paste0("'(' * italic(n) * '=' * ", n, " * ')'"))

################################################################################
#
# Beta-diversity plot (Figure 3)
#
################################################################################

fill_vals <- c("Seed_Dispersal" = "#a1dab4", "Pollination" = "#fdbb84")
col_vals  <- c("Seed_Dispersal" = "#1b9e77", "Pollination" = "#d95f02")

p_beta_dot <- ggplot(combined_df, aes(x = ComboInteraction, y = beta)) +
  geom_jitter(aes(color = Interaction_key),
              width = 0.2, size = 1.8, alpha = 0.7) +
  geom_boxplot(aes(fill = Interaction_key),
               width = 0.5, alpha = 0.4,
               outlier.shape = NA,
               color = "black", size = 0.6) +
  geom_text(
    data = combo_counts,
    aes(
      x = ComboInteraction,
      y = max(combined_df$beta, na.rm = TRUE) + 0.05,
      label = label
    ),
    inherit.aes = FALSE,
    size = 4.5,
    parse = TRUE
  ) +
  scale_fill_manual(
    values = fill_vals,
    breaks = names(display_map),
    labels = unname(display_map)
  ) +
  scale_color_manual(
    values = col_vals,
    breaks = names(display_map),
    labels = unname(display_map)
  ) +
  scale_x_discrete(labels = combo_labels) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(y = "Beta-diversity (Jaccard)", x = NULL) +
  theme_bw() +
  theme(
    axis.text.x     = element_text(size = 16, color = "black", lineheight = 0.9),
    axis.text.y     = element_text(size = 18, color = "black"),
    axis.title.y    = element_text(size = 18, face = "bold"),
    legend.position = "top",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 18),
    panel.border    = element_rect(colour = "grey", fill = NA, size = 0.4)
  ) +
  geom_vline(xintercept = separator_position,
             linetype = "dashed", size = 0.7) +
  annotate(
    "text",
    x = separator_position - 0.8,
    y = max(combined_df$beta, na.rm = TRUE) + 0.12,
    label = "One network per publication",
    hjust = 1, size = 5.5, fontface = "bold"
  ) +
  annotate(
    "text",
    x = separator_position + 0.8,
    y = max(combined_df$beta, na.rm = TRUE) + 0.12,
    label = "Multiple networks per publication",
    hjust = 0, size = 5.5, fontface = "bold"
  )
print(p_beta_dot)

# Save final plot as PDF
ggsave(
  "beta_diversity.pdf",
  device = cairo_pdf,
  width = 12,
  height = 5,
  units = "in"
)

################################################################################
#
# Beta-diversity table
#
################################################################################

# The two sets
group_map <- c(
  "Within Publication"          = "Within publication",
  "One network per publication" = "Between publications"
)

# Pairwsie-beta diversity between all networks
combined_df_processed <- combined_df %>%
  mutate(
    Set   = recode(Group, !!!group_map),
    Guild = if_else(DataType == "row", "Plant guild", "Animal guild"),
    Interaction = if ("Interaction_key" %in% names(.)) Interaction_key else Interaction  # safety
  ) %>%
  filter(!is.na(Set))

# Median beta-diversity within/between by interaction type and guild
beta_medians_tbl <- combined_df_processed %>%
  group_by(Set, Interaction, Guild) %>%
  summarise(
    median_beta = median(beta, na.rm = TRUE),
    n_pairs     = sum(!is.na(beta)),
    .groups     = "drop"
  )

# Number of networks within/between by interaction type and guild
network_counts_tbl <- combined_df_processed %>%
  select(Set, Interaction, ID1, ID2) %>%
  pivot_longer(cols = c(ID1, ID2), values_to = "ID_with_suffix") %>%
  mutate(
    ID = sub("(_row_names|_column_names)$", "", ID_with_suffix)
  ) %>%
  distinct(Set, Interaction, ID) %>%
  count(Set, Interaction, name = "n_networks")

# Join network counts onto the medians table
final_within_between_tbl <- beta_medians_tbl %>%
  left_join(network_counts_tbl, by = c("Set", "Interaction")) %>%
  mutate(
    median_beta = round(median_beta, 2),
    n_pairs     = as.integer(n_pairs),
    n_networks  = as.integer(n_networks)
  ) %>%
  arrange(Set, Guild, Interaction)

print(final_within_between_tbl)
