library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(tidyr)
library(here)
library(readr)

################################################################################
#
# Data
#
################################################################################

metadata <- read_csv(here("general_network_information.csv"))
fricke_metadata <- read_csv(here("fricke_metadata.csv"))
pollination_metadata <- read_csv(here("pollination_sampling_metadata.csv"))

# Making sure Fricke networks are labelled the same as pollination 
fricke_metadata$net_id <- gsub(" ", "_", fricke_metadata$net_id)
fricke_metadata$net_id <- gsub("-", "_", fricke_metadata$net_id)
fricke_metadata$net_id <- gsub("&", "and", fricke_metadata$net_id)
fricke_metadata$net_id <- gsub("\\.", "", fricke_metadata$net_id)

# Combining metadata of seed-dispersal and pollination
both_combined <- left_join(metadata,fricke_metadata,by=c("ID"="net_id"))

# Overwrite metadata fields for pollination networks using curated pollination 
# metadata
both_combined <- both_combined %>%
  rows_update(
    pollination_metadata %>% select(any_of(names(both_combined))),
    by = "ID"
  )

# Network directories
net_dir <- "networks"
row_dir <- file.path(net_dir, "row_names")
col_dir <- file.path(net_dir, "column_names")
dir.create(row_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(col_dir, recursive = TRUE, showWarnings = FALSE)

# Network CSVs and their network IDs
all_csv_paths <- list.files(net_dir, pattern = "\\.csv$", full.names = TRUE)
all_csv_ids   <- sub("\\.csv$", "", basename(all_csv_paths))

# Paths for networks retained in both_combined
use_ids   <- both_combined$ID
use_paths <- file.path(net_dir, paste0(use_ids, ".csv"))

################################################################################
#
# Counting unidentified species
#
################################################################################

# Function to flag labels indicating unidentified taxa based on common 
# placeholder patterns
is_unidentified <- function(x) {
  grepl("unidsp|\\d", x, ignore.case = TRUE)
}

# Extract row and column species names from a network CSV
extract_species_names <- function(path) {
  df <- read_csv(path)
  
  list(
    row_species = trimws(as.character(df[[1]])),
    col_species = trimws(colnames(df)[-1])
  )
}

# Summarize which networks each species appears in by collapsing network IDs
collapse_networks <- function(df) {
  df %>%
    group_by(species) %>%
    summarise(networks = paste(sort(unique(ID)), collapse = "; "), .groups = "drop") %>%
    arrange(species)
}

# Storage for number of nodes down to species and number of nodes not down to
# species
counts_df <- data.frame(
  ID = use_ids,
  identified_row_species_n    = integer(length(use_ids)),
  unidentified_row_species_n  = integer(length(use_ids)),
  identified_col_species_n    = integer(length(use_ids)),
  unidentified_col_species_n  = integer(length(use_ids)),
  stringsAsFactors = FALSE
)

# Storage for the names of nodes not identified down to species (and which 
# they appear in networks)
name_listing_unid <- list(
  row = data.frame(species = character(), ID = character(), stringsAsFactors = FALSE),
  col = data.frame(species = character(), ID = character(), stringsAsFactors = FALSE)
)

# Storage for the names of nodes identified down to species (and which they 
# appear in networks)
name_listing_ident <- list(
  row = data.frame(species = character(), ID = character(), stringsAsFactors = FALSE),
  col = data.frame(species = character(), ID = character(), stringsAsFactors = FALSE)
)

# Iterate over networks to count identified vs unidentified node labels 
# (row/col) and record their names by network ID
for (i in seq_along(use_paths)) {
  id <- use_ids[i]
  sp <- extract_species_names(use_paths[i])
  
  row_sp <- sp$row_species
  col_sp <- sp$col_species
  
  # Clean empties
  row_sp <- row_sp[!is.na(row_sp) & nzchar(row_sp)]
  col_sp <- col_sp[!is.na(col_sp) & nzchar(col_sp)]
  
  # Totals per axis
  row_total <- length(row_sp)
  col_total <- length(col_sp)
  
  # Unidentified counts
  row_unid <- sum(is_unidentified(row_sp), na.rm = TRUE)
  col_unid <- sum(is_unidentified(col_sp), na.rm = TRUE)
  
  # Identified = total - unidentified
  row_id <- row_total - row_unid
  col_id <- col_total - col_unid
  
  counts_df$identified_row_species_n[i]   <- row_id
  counts_df$unidentified_row_species_n[i] <- row_unid
  counts_df$identified_col_species_n[i]   <- col_id
  counts_df$unidentified_col_species_n[i] <- col_unid
  
  if (row_total > 0) {
    # rows unidentified
    ru <- row_sp[is_unidentified(row_sp)]
    if (length(ru)) {
      name_listing_unid$row <- rbind(name_listing_unid$row, data.frame(species = ru, ID = id, stringsAsFactors = FALSE))
    }
    # rows identified
    ri <- row_sp[!is_unidentified(row_sp)]
    if (length(ri)) {
      name_listing_ident$row <- rbind(name_listing_ident$row, data.frame(species = ri, ID = id, stringsAsFactors = FALSE))
    }
  }
  if (col_total > 0) {
    # cols unidentified
    cu <- col_sp[is_unidentified(col_sp)]
    if (length(cu)) {
      name_listing_unid$col <- rbind(name_listing_unid$col, data.frame(species = cu, ID = id, stringsAsFactors = FALSE))
    }
    # cols identified
    ci <- col_sp[!is_unidentified(col_sp)]
    if (length(ci)) {
      name_listing_ident$col <- rbind(name_listing_ident$col, data.frame(species = ci, ID = id, stringsAsFactors = FALSE))
    }
  }
}

# Combine information about all networks (i.e., metadata) and the number of 
# identified and non-identified nodes down to species
both_combined <- both_combined %>%
  left_join(counts_df, by = "ID") %>%
  mutate(
    identified_row_species_n    = ifelse(is.na(identified_row_species_n),    0L, identified_row_species_n),
    unidentified_row_species_n  = ifelse(is.na(unidentified_row_species_n),  0L, unidentified_row_species_n),
    identified_col_species_n    = ifelse(is.na(identified_col_species_n),    0L, identified_col_species_n),
    unidentified_col_species_n  = ifelse(is.na(unidentified_col_species_n),  0L, unidentified_col_species_n)
  )

# Compute per-network % unidentified (rows+cols combined)
pct_df <- both_combined %>%
  mutate(
    row_total = identified_row_species_n + unidentified_row_species_n,
    col_total = identified_col_species_n + unidentified_col_species_n,
    species_total = row_total + col_total,
    unidentified_total = unidentified_row_species_n + unidentified_col_species_n,
    unidentified_pct = ifelse(species_total > 0, 100 * unidentified_total / species_total, NA_real_)
  ) %>%
  filter(is.finite(unidentified_pct))

################################################################################
#
# Histogram of the number of networks and their degree (%) of nodes not 
# identified to species level (Figure S2)
#
################################################################################

# Define breaks, labels, and colors
breaks <- seq(0, 100, by = 10)
labels <- paste0(head(breaks, -1), "–", tail(breaks, -1), "%")
fill_vals  <- c("Pollination" = "#fdbb84", "Seed Dispersal" = "#a1dab4")
color_vals <- c("Pollination" = "#d95f02", "Seed Dispersal" = "#1b9e77")

# Bining data
binned_all <- pct_df %>%
  mutate(
    bin = case_when(
      unidentified_pct == 0 ~ "0%",
      TRUE ~ as.character(cut(unidentified_pct,
                              breaks = breaks,
                              right = FALSE,
                              labels = labels,
                              include.lowest = FALSE))
    ),
    bin = factor(bin, levels = c("0%", labels))
  ) %>%
  group_by(TYPE, bin) %>%
  tally(name = "n") %>%
  ungroup() %>%
  complete(bin, TYPE, fill = list(n = 0))

# Keep only bins that have any observations across both types
present_bins_tbl <- binned_all %>%
  group_by(bin) %>% 
  summarise(total = sum(n), .groups = "drop") %>%
  filter(total > 0)

# Restrict binned data to observed bins and retain network type for plotting
plot_df <- binned_all %>%
  semi_join(present_bins_tbl, by = "bin") %>%
  mutate(
    TYPE = factor(TYPE, levels = c("Seed Dispersal", "Pollination")),
    bin  = forcats::fct_drop(bin)
  )

# Histogram plot
p_combined <- ggplot(plot_df, aes(x = bin, y = n, fill = TYPE)) +
  geom_col(
    aes(color = ifelse(n > 0, as.character(TYPE), NA)),
    width = 0.9, 
    position = position_dodge2(preserve = "single")
  ) +
  scale_fill_manual(
    values = fill_vals,
    breaks = c("Pollination", "Seed Dispersal"),
    labels = c("Pollination", "Seed-dispersal"),
    name   = NULL
  ) +
  scale_color_manual(
    values = color_vals,
    guide  = "none"
  ) +
  labs(x = "Taxonomically unresolved nodes to species (%) per network", y = "Number of networks") +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "top",
    axis.text.x     = element_text(size = 12, colour = "black"),
    axis.text.y     = element_text(size = 12, colour = "black"),
    axis.title.x    = element_text(face = "bold"),
    axis.title.y    = element_text(face = "bold"),
    axis.ticks      = element_line(colour = "black"),            # tick marks black
    panel.border    = element_rect(colour = "grey60", fill = NA, linewidth = 0.5)
  )

print(p_combined)

ggsave("unidentified_species_percent.pdf", plot = last_plot(), device = "pdf", width = 10, height = 5, units = "in", dpi = 600)

################################################################################
#
# Figure 5: Boxplot of publication groupings and the percentage of nodes lacking
# species level taxonomic identification (%) per network
#
################################################################################

# Average unidentified per publication grouping
avg_by_publication <- both_combined %>%
  mutate(
    row_total = identified_row_species_n + unidentified_row_species_n,
    col_total = identified_col_species_n + unidentified_col_species_n,
    species_total = row_total + col_total,
    unidentified_total = unidentified_row_species_n + unidentified_col_species_n,
    unidentified_pct = ifelse(species_total > 0, 100 * unidentified_total / species_total, NA_real_),
    TYPE2 = ifelse(TYPE == "Pollination", "Pollination", "Seed Dispersal"),
    Publication_split = ifelse(
      grepl("^\\s*one\\s+network\\s+per\\s+publication\\s*$", Publication, ignore.case = TRUE),
      paste0("One network per publication — ", TYPE2),
      Publication
    )
  ) %>%
  filter(is.finite(unidentified_pct)) %>%
  group_by(Publication_split) %>%
  summarise(
    n_networks           = n(),
    avg_unidentified_pct = round(mean(unidentified_pct), 1),
    sd_unidentified_pct  = round(sd(unidentified_pct), 1),
    types_present        = paste(sort(unique(TYPE2)), collapse = " & "),
    .groups = "drop"
  ) %>%
  arrange(desc(avg_unidentified_pct))
avg_by_publication

# Fix accents in Publication names (UTF-8 safe)
both_combined$Publication <- enc2utf8(both_combined$Publication)
both_combined <- both_combined %>%
  mutate(
    Publication = Publication %>%
      str_replace_all(fixed("Hernandez-Montero"), "Hernández-Montero") %>%
      str_replace_all(fixed("Garcıa"),            "García") %>%
      str_replace_all(fixed("Quitian"),           "Quitián") %>%
      str_replace_all(c(
        "´a" = "á", "´e" = "é", "´i" = "í", "´ı" = "í", "´o" = "ó", "´u" = "ú",
        "´A" = "Á", "´E" = "É", "´I" = "Í", "´O" = "Ó", "´U" = "Ú"
      ))
  )

# Build plotting dataset
both_combined_for_plot <- both_combined %>%
  mutate(
    row_total          = identified_row_species_n + unidentified_row_species_n,
    col_total          = identified_col_species_n + unidentified_col_species_n,
    species_total      = row_total + col_total,
    unidentified_total = unidentified_row_species_n + unidentified_col_species_n,
    unidentified_pct   = if_else(species_total > 0,
                                 100 * unidentified_total / species_total,
                                 NA_real_),
    TYPE2 = if_else(TYPE == "Pollination", "Pollination", "Seed Dispersal")
  ) %>%
  filter(is.finite(unidentified_pct)) %>%
  mutate(
    is_single = grepl("^\\s*one\\s+network\\s+per\\s+publication\\s*$",
                      Publication, ignore.case = TRUE),
    Publication_key  = if_else(is_single,
                               paste0("One network per publication — ", TYPE2),
                               Publication),
    Publication_tick = if_else(is_single,
                               "One network per publication",
                               Publication),
    pub_group = if_else(is_single, "Single", "Multiple"),
    TYPE2 = factor(TYPE2, levels = c("Seed Dispersal", "Pollination"))
  )

# Medians per x-level (Publication_key) for ordering
pub_summ <- both_combined_for_plot %>%
  group_by(Publication_key, Publication_tick, pub_group) %>%
  summarise(med = median(unidentified_pct), .groups = "drop")

# Order multi-network publications by descending median unidentified percentage
multi_ord <- pub_summ %>%
  filter(pub_group == "Multiple") %>%
  arrange(desc(med)) %>%
  pull(Publication_key)

# Order the two single-network-per-publication ticks explicitly (Pollination 
# first, then Seed Dispersal)
single_ord <- pub_summ %>%
  filter(pub_group == "Single") %>%
  mutate(TYPE2 = if_else(grepl("—\\s*Pollination$", Publication_key),
                         "Pollination", "Seed Dispersal")) %>%
  arrange(factor(TYPE2, levels = c("Pollination", "Seed Dispersal"))) %>%
  pull(Publication_key)

# Combine ordering so multi-publication entries appear first, followed by the 
# two single-publication ticks
levels_vec <- c(multi_ord, single_ord)

# Apply the custom x-axis ordering to the plotting data
both_combined_for_plot <- both_combined_for_plot %>%
  mutate(Publication_key = factor(Publication_key, levels = levels_vec))

# Position of the separator line between multi- and single-publication groups
split_pos <- length(multi_ord)

# Tick labels with [n=..] per tick
label_vec <- both_combined_for_plot %>%
  count(Publication_key, Publication_tick, name = "n_networks") %>%
  mutate(label = paste0(Publication_tick, " [n=", n_networks, "]")) %>%
  { setNames(.$label, .$Publication_key) }

# Plot
p_unid <- ggplot(both_combined_for_plot,
                 aes(x = Publication_key, y = unidentified_pct,
                     fill = TYPE2, color = TYPE2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4, size = 0.6,
               position = position_dodge2(preserve = "single")) +
  geom_point(size = 1.8, alpha = 0.7,
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  coord_flip() +
  scale_fill_manual(
    values = c("Seed Dispersal" = "#a1dab4",
               "Pollination"   = "#fdbb84"),
    breaks = c("Pollination", "Seed Dispersal"),
    labels = c("Pollination", "Seed-dispersal")
  ) +
  scale_color_manual(
    values = c("Seed Dispersal" = "#1b9e77",
               "Pollination"   = "#d95f02"),
    breaks = c("Pollination", "Seed Dispersal"),
    labels = c("Pollination", "Seed-dispersal")
  ) +
  labs(
    x = NULL,
    y = "Nodes lacking species level identification (%) per network",
    fill = "Network type", color = "Network type"
  ) +
  scale_x_discrete(labels = label_vec) +
  theme_bw(base_size = 14) +
  geom_vline(xintercept = split_pos + 0.5, linetype = "dashed", color = "black")

print(p_unid)

ggsave("unidentified_species_percent_per_publication.pdf", plot = last_plot(), device = "pdf", width = 9, height = 7, units = "in", dpi = 600)

################################################################################
#
# Percent unidentified taxonomy (e.g., family, genus) for plant and animal 
# guilds of pollination and seed-dispersal (Table S2)
#
################################################################################

seed_dispersal_row <- read_csv(here("networks", "row_names", "seed_dispersal_unique_unidentified_row_species_names.csv"))
seed_dispersal_column <- read_csv(here("networks", "column_names", "seed_dispersal_unique_unidentified_column_species_names.csv"))
pollination_row <- read_csv(here("networks", "row_names", "pollination_unique_unidentified_row_species_names.csv"))
pollination_column <- read_csv(here("networks", "column_names", "pollination_unique_unidentified_column_species_names.csv"))

taxonomy_summary <- function(df, filename) {

  df <- subset(df, Taxonomy != "")
  
  counts <- table(df$Taxonomy)
  percentages <- prop.table(counts) * 100
  
  data.frame(
    File       = filename,
    Taxonomy   = names(counts),
    Count      = as.vector(counts),
    Percentage = round(as.vector(percentages), 2)
  )
}

# Apply to each file
out_seed_row     <- taxonomy_summary(seed_dispersal_row, "seed_dispersal_row")
out_seed_column  <- taxonomy_summary(seed_dispersal_column, "seed_dispersal_column")
out_poll_row     <- taxonomy_summary(pollination_row, "pollination_row")
out_poll_column  <- taxonomy_summary(pollination_column, "pollination_column")
all_summaries <- rbind(out_seed_row, out_seed_column, out_poll_row, out_poll_column)

all_summaries
