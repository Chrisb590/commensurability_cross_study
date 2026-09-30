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

# Directory where all figures and tables are written (created if missing)
output_dir <- here("output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

metadata <- read_csv(
  here("data","general_network_information.csv"),
  show_col_types = FALSE
)

fricke_metadata <- read_csv(
  here("data","fricke_metadata.csv"),
  show_col_types = FALSE
)

pollination_metadata <- read_csv(
  here("data","pollination_sampling_metadata.csv"),
  show_col_types = FALSE
)

# Standardize Fricke network IDs to match IDs used in the main metadata
fricke_metadata <- fricke_metadata %>%
  mutate(
    net_id = net_id %>%
      str_replace_all(" ", "_") %>%
      str_replace_all("-", "_") %>%
      str_replace_all("&", "and") %>%
      str_replace_all("\\.", "")
  )

# Combine seed-dispersal metadata with the main network metadata
both_combined <- metadata %>%
  left_join(
    fricke_metadata,
    by = c(
      "ID" = "net_id"
    )
  )

# Replace metadata fields for pollination networks with curated pollination
# metadata
both_combined <- both_combined %>%
  rows_update(
    pollination_metadata %>%
      select(
        any_of(
          names(both_combined)
        )
      ),
    by = "ID"
  )

# Standardize interaction-type labels for matching curated taxonomic lists
both_combined <- both_combined %>%
  mutate(
    Interaction_key = case_when(
      
      str_detect(
        TYPE,
        regex(
          "pollination",
          ignore_case = TRUE
        )
      ) ~
        "Pollination",
      
      str_detect(
        TYPE,
        regex(
          "seed",
          ignore_case = TRUE
        )
      ) ~
        "Seed_Dispersal",
      
      TRUE ~
        NA_character_
    )
  )

# Standardize accented publication names used in Figure 5
both_combined$Publication <-
  enc2utf8(
    both_combined$Publication
  )

both_combined <- both_combined %>%
  mutate(
    Publication = Publication %>%
      str_replace_all(
        fixed("Hernandez-Montero"),
        "Hernández-Montero"
      ) %>%
      str_replace_all(
        fixed("Garcıa"),
        "García"
      ) %>%
      str_replace_all(
        fixed("Quitian"),
        "Quitián"
      ) %>%
      str_replace_all(
        c(
          "´a" = "á",
          "´e" = "é",
          "´i" = "í",
          "´ı" = "í",
          "´o" = "ó",
          "´u" = "ú",
          "´A" = "Á",
          "´E" = "É",
          "´I" = "Í",
          "´O" = "Ó",
          "´U" = "Ú"
        )
      )
  )

# Network and curated taxonomic-list directories
net_dir <- here("data","networks")

row_dir <- here(
  "data",
  "networks",
  "row_names"
)

col_dir <- here(
  "data",
  "networks",
  "column_names"
)

# Networks included in the analysis
use_ids <- both_combined$ID

use_paths <- file.path(
  net_dir,
  paste0(
    use_ids,
    ".csv"
  )
)

################################################################################
#
# Main manuscript analyses:
# Taxonomic resolution of quantitative interaction networks
#
# Evaluate the taxonomic resolution of nodes across pollination and
# seed-dispersal networks using manually curated lists of nodes not identified
# to species level. For each network, calculate the number and percentage of
# unresolved nodes across the plant and animal guilds.
#
# The distribution of unresolved nodes across all networks is reported in
# Figure S2, while variation among publication groups is shown in Figure 5.
# The taxonomic resolution of unresolved nodes is summarized in Table S2.
#
################################################################################

################################################################################
#
# Read curated lists of nodes not identified to species level
#
# These files contain the manually curated node names that are not resolved to
# species level, together with their taxonomic resolution where available.
#
# The species column is used to identify unresolved nodes throughout the
# analysis, while the Taxonomy column is used to summarize their taxonomic
# resolution in Table S2.
#
################################################################################

unresolved_pollination_row <- read_csv(
  file.path(
    row_dir,
    "pollination_unique_unidentified_row_species_names.csv"
  ),
  show_col_types = FALSE
)

unresolved_pollination_column <- read_csv(
  file.path(
    col_dir,
    "pollination_unique_unidentified_column_species_names.csv"
  ),
  show_col_types = FALSE
)

unresolved_seed_row <- read_csv(
  file.path(
    row_dir,
    "seed_dispersal_unique_unidentified_row_species_names.csv"
  ),
  show_col_types = FALSE
)

unresolved_seed_column <- read_csv(
  file.path(
    col_dir,
    "seed_dispersal_unique_unidentified_column_species_names.csv"
  ),
  show_col_types = FALSE
)

# Create lookup lists from the manually curated species names
unresolved <- list(
  
  Pollination_row =
    unresolved_pollination_row %>%
    pull(species) %>%
    str_squish() %>%
    unique(),
  
  Pollination_column =
    unresolved_pollination_column %>%
    pull(species) %>%
    str_squish() %>%
    unique(),
  
  Seed_Dispersal_row =
    unresolved_seed_row %>%
    pull(species) %>%
    str_squish() %>%
    unique(),
  
  Seed_Dispersal_column =
    unresolved_seed_column %>%
    pull(species) %>%
    str_squish() %>%
    unique()
)

################################################################################
#
# Extract node names from quantitative interaction networks
#
# Row names represent the plant guild and column names represent the animal
# guild. Node labels are cleaned before comparison with the curated unresolved
# taxonomic lists.
#
################################################################################

extract_species_names <- function(path) {
  
  df <- read_csv(
    path,
    show_col_types = FALSE,
    name_repair = "minimal"
  )
  
  list(
    row_species =
      str_squish(
        as.character(
          df[[1]]
        )
      ),
    
    col_species =
      str_squish(
        colnames(df)[-1]
      )
  )
}

################################################################################
#
# Count identified and unresolved nodes within each network
#
# For each network, classify every plant and animal node using the appropriate
# curated unresolved-node list and calculate the number resolved and unresolved
# to species level within each guild.
#
################################################################################

counts_df <- data.frame(
  
  ID =
    use_ids,
  
  identified_row_species_n =
    integer(
      length(use_ids)
    ),
  
  unidentified_row_species_n =
    integer(
      length(use_ids)
    ),
  
  identified_col_species_n =
    integer(
      length(use_ids)
    ),
  
  unidentified_col_species_n =
    integer(
      length(use_ids)
    ),
  
  stringsAsFactors =
    FALSE
)

for (i in seq_along(use_paths)) {
  
  id <- use_ids[i]
  
  interaction_key <-
    both_combined$Interaction_key[
      match(
        id,
        both_combined$ID
      )
    ]
  
  sp <- extract_species_names(
    use_paths[i]
  )
  
  row_sp <- sp$row_species
  col_sp <- sp$col_species
  
  # Remove missing or blank node labels
  row_sp <- row_sp[
    !is.na(row_sp) &
      nzchar(row_sp)
  ]
  
  col_sp <- col_sp[
    !is.na(col_sp) &
      nzchar(col_sp)
  ]
  
  # Select the appropriate curated unresolved-node lists
  if (interaction_key == "Pollination") {
    
    unresolved_row <-
      unresolved$Pollination_row
    
    unresolved_col <-
      unresolved$Pollination_column
    
  } else if (interaction_key == "Seed_Dispersal") {
    
    unresolved_row <-
      unresolved$Seed_Dispersal_row
    
    unresolved_col <-
      unresolved$Seed_Dispersal_column
    
  }
  
  # Identify unresolved nodes
  row_unidentified <-
    tolower(
      str_squish(
        row_sp
      )
    ) %in%
    tolower(
      str_squish(
        unresolved_row
      )
    )
  
  col_unidentified <-
    tolower(
      str_squish(
        col_sp
      )
    ) %in%
    tolower(
      str_squish(
        unresolved_col
      )
    )
  
  # Count identified and unresolved nodes
  counts_df$unidentified_row_species_n[i] <-
    sum(
      row_unidentified
    )
  
  counts_df$identified_row_species_n[i] <-
    sum(
      !row_unidentified
    )
  
  counts_df$unidentified_col_species_n[i] <-
    sum(
      col_unidentified
    )
  
  counts_df$identified_col_species_n[i] <-
    sum(
      !col_unidentified
    )
}

################################################################################
#
# Calculate percentage of unresolved nodes per network
#
# Join node counts to the network metadata and calculate the percentage of all
# plant and animal nodes in each network that are not identified to species
# level.
#
################################################################################

both_combined <- both_combined %>%
  left_join(
    counts_df,
    by = "ID"
  )

pct_df <- both_combined %>%
  mutate(
    
    row_total =
      identified_row_species_n +
      unidentified_row_species_n,
    
    col_total =
      identified_col_species_n +
      unidentified_col_species_n,
    
    species_total =
      row_total +
      col_total,
    
    unidentified_total =
      unidentified_row_species_n +
      unidentified_col_species_n,
    
    unidentified_pct =
      if_else(
        species_total > 0,
        100 *
          unidentified_total /
          species_total,
        NA_real_
      ),
    
    TYPE_plot =
      case_when(
        
        Interaction_key ==
          "Pollination" ~
          "Pollination",
        
        Interaction_key ==
          "Seed_Dispersal" ~
          "Seed Dispersal"
      )
  ) %>%
  
  filter(
    is.finite(
      unidentified_pct
    )
  )

################################################################################
#
# Figure S2:
# Distribution of taxonomically unresolved nodes among networks
#
# Group networks according to the percentage of plant and animal nodes that are
# not identified to species level and compare the resulting distributions for
# pollination and seed-dispersal networks. Zero-percent networks are shown as a
# separate category. All remaining percentage intervals are open on the right.
#
################################################################################

breaks <- seq(
  0,
  100,
  by = 10
)

labels <- paste0(
  head(
    breaks,
    -1
  ),
  "–",
  tail(
    breaks,
    -1
  ),
  "%"
)

fill_vals <- c(
  "Pollination" =
    "#fdbb84",
  "Seed Dispersal" =
    "#a1dab4"
)

color_vals <- c(
  "Pollination" =
    "#d95f02",
  "Seed Dispersal" =
    "#1b9e77"
)

# Assign networks to right-open percentage bins
binned_all <- pct_df %>%
  mutate(
    
    bin = case_when(
      
      unidentified_pct == 0 ~
        "0%",
      
      TRUE ~
        as.character(
          cut(
            unidentified_pct,
            breaks =
              breaks,
            right =
              FALSE,
            labels =
              labels,
            include.lowest =
              FALSE
          )
        )
    ),
    
    bin = factor(
      bin,
      levels = c(
        "0%",
        labels
      )
    )
  ) %>%
  
  group_by(
    TYPE_plot,
    bin
  ) %>%
  
  tally(
    name = "n"
  ) %>%
  
  ungroup() %>%
  
  complete(
    bin,
    TYPE_plot,
    fill = list(
      n = 0
    )
  )

# Retain only bins containing at least one network
present_bins_tbl <- binned_all %>%
  group_by(
    bin
  ) %>%
  summarise(
    total =
      sum(n),
    .groups =
      "drop"
  ) %>%
  filter(
    total > 0
  )

plot_df <- binned_all %>%
  semi_join(
    present_bins_tbl,
    by = "bin"
  ) %>%
  mutate(
    
    TYPE_plot = factor(
      TYPE_plot,
      levels = c(
        "Seed Dispersal",
        "Pollination"
      )
    ),
    
    bin =
      forcats::fct_drop(
        bin
      )
  )

p_combined <- ggplot(
  plot_df,
  aes(
    x =
      bin,
    y =
      n,
    fill =
      TYPE_plot
  )
) +
  
  geom_col(
    aes(
      color = ifelse(
        n > 0,
        as.character(
          TYPE_plot
        ),
        NA
      )
    ),
    width =
      0.9,
    position =
      position_dodge2(
        preserve =
          "single"
      )
  ) +
  
  scale_fill_manual(
    values =
      fill_vals,
    breaks = c(
      "Pollination",
      "Seed Dispersal"
    ),
    labels = c(
      "Pollination",
      "Seed-dispersal"
    ),
    name =
      NULL
  ) +
  
  scale_color_manual(
    values =
      color_vals,
    guide =
      "none"
  ) +
  
  labs(
    x =
      "Nodes lacking species-level identification (%) per network",
    y =
      "Number of networks"
  ) +
  
  theme_minimal(
    base_size =
      16
  ) +
  
  theme(
    
    legend.position =
      "top",
    
    axis.text.x =
      element_text(
        size = 12,
        colour = "black"
      ),
    
    axis.text.y =
      element_text(
        size = 12,
        colour = "black"
      ),
    
    axis.title.x =
      element_text(
        face = "bold"
      ),
    
    axis.title.y =
      element_text(
        face = "bold"
      ),
    
    axis.ticks =
      element_line(
        colour = "black"
      ),
    
    panel.border =
      element_rect(
        colour = "grey60",
        fill = NA,
        linewidth = 0.5
      )
  )

print(
  p_combined
)

ggsave(
  file.path(output_dir,"figure_s2.pdf"),
  plot =
    p_combined,
  device =
    "pdf",
  width =
    10,
  height =
    5,
  units =
    "in",
  dpi =
    600
)

################################################################################
#
# Figure 5:
# Taxonomic resolution among publication groups
#
# Compare the percentage of unresolved nodes among publication groups.
# Publications contributing multiple networks are shown separately, while
# networks originating from single-network publications are pooled separately
# for pollination and seed-dispersal networks.
#
################################################################################

figure5_df <- pct_df %>%
  mutate(
    
    is_single =
      str_detect(
        Publication,
        regex(
          "^\\s*one\\s+network\\s+per\\s+publication\\s*$",
          ignore_case = TRUE
        )
      ),
    
    Publication_key =
      if_else(
        is_single,
        paste0(
          "One network per publication — ",
          TYPE_plot
        ),
        Publication
      ),
    
    Publication_tick =
      if_else(
        is_single,
        "One network per publication",
        Publication
      ),
    
    pub_group =
      if_else(
        is_single,
        "Single",
        "Multiple"
      ),
    
    TYPE_plot = factor(
      TYPE_plot,
      levels = c(
        "Seed Dispersal",
        "Pollination"
      )
    )
  )

################################################################################
#
# Publication-level descriptive statistics for Figure 5
#
# Calculate the number of networks and the mean and SD percentage of unresolved
# nodes represented by each publication grouping.
#
################################################################################

avg_by_publication <- figure5_df %>%
  group_by(
    Publication_key
  ) %>%
  summarise(
    
    n_networks =
      n(),
    
    avg_unidentified_pct =
      round(
        mean(
          unidentified_pct,
          na.rm = TRUE
        ),
        1
      ),
    
    sd_unidentified_pct =
      round(
        sd(
          unidentified_pct,
          na.rm = TRUE
        ),
        1
      ),
    
    types_present =
      paste(
        sort(
          unique(
            as.character(
              TYPE_plot
            )
          )
        ),
        collapse = " & "
      ),
    
    .groups =
      "drop"
  ) %>%
  
  arrange(
    desc(
      avg_unidentified_pct
    )
  )

avg_by_publication

# Order publication groups by median percentage unresolved
pub_summ <- figure5_df %>%
  group_by(
    Publication_key,
    Publication_tick,
    pub_group
  ) %>%
  summarise(
    med =
      median(
        unidentified_pct,
        na.rm = TRUE
      ),
    .groups =
      "drop"
  )

# Multi-network publications are ordered by descending median
multi_ord <- pub_summ %>%
  filter(
    pub_group ==
      "Multiple"
  ) %>%
  arrange(
    desc(
      med
    )
  ) %>%
  pull(
    Publication_key
  )

# Pooled single-network publication groups are shown last
single_ord <- pub_summ %>%
  filter(
    pub_group ==
      "Single"
  ) %>%
  mutate(
    TYPE_plot =
      if_else(
        str_detect(
          Publication_key,
          "—\\s*Pollination$"
        ),
        "Pollination",
        "Seed Dispersal"
      )
  ) %>%
  arrange(
    factor(
      TYPE_plot,
      levels = c(
        "Pollination",
        "Seed Dispersal"
      )
    )
  ) %>%
  pull(
    Publication_key
  )

levels_vec <- c(
  multi_ord,
  single_ord
)

figure5_df <- figure5_df %>%
  mutate(
    Publication_key = factor(
      Publication_key,
      levels =
        levels_vec
    )
  )

# Position separating multiple- and single-network publication groups
split_pos <-
  length(
    multi_ord
  )

# Add the number of networks represented by each publication group
label_vec <- figure5_df %>%
  count(
    Publication_key,
    Publication_tick,
    name =
      "n_networks"
  ) %>%
  mutate(
    label =
      paste0(
        Publication_tick,
        " [n=",
        n_networks,
        "]"
      )
  ) %>%
  {
    setNames(
      .$label,
      .$Publication_key
    )
  }

p_unid <- ggplot(
  figure5_df,
  aes(
    x =
      Publication_key,
    y =
      unidentified_pct,
    fill =
      TYPE_plot,
    color =
      TYPE_plot
  )
) +
  
  geom_boxplot(
    outlier.shape =
      NA,
    alpha =
      0.4,
    size =
      0.6,
    position =
      position_dodge2(
        preserve =
          "single"
      )
  ) +
  
  geom_point(
    size =
      1.8,
    alpha =
      0.7,
    position =
      position_jitterdodge(
        jitter.width =
          0.2,
        dodge.width =
          0.8
      )
  ) +
  
  coord_flip() +
  
  scale_fill_manual(
    values = c(
      "Seed Dispersal" =
        "#a1dab4",
      "Pollination" =
        "#fdbb84"
    ),
    breaks = c(
      "Pollination",
      "Seed Dispersal"
    ),
    labels = c(
      "Pollination",
      "Seed-dispersal"
    )
  ) +
  
  scale_color_manual(
    values = c(
      "Seed Dispersal" =
        "#1b9e77",
      "Pollination" =
        "#d95f02"
    ),
    breaks = c(
      "Pollination",
      "Seed Dispersal"
    ),
    labels = c(
      "Pollination",
      "Seed-dispersal"
    )
  ) +
  
  labs(
    x =
      NULL,
    y =
      "Nodes lacking species-level identification (%) per network",
    fill =
      "Network type",
    color =
      "Network type"
  ) +
  
  scale_x_discrete(
    labels =
      label_vec
  ) +
  
  theme_bw(
    base_size =
      14
  ) +
  
  geom_vline(
    xintercept =
      split_pos + 0.5,
    linetype =
      "dashed",
    color =
      "black"
  )

print(
  p_unid
)

ggsave(
  file.path(output_dir,"figure_5.pdf"),
  plot =
    p_unid,
  device =
    "pdf",
  width =
    9,
  height =
    7,
  units =
    "in",
  dpi =
    600
)

################################################################################
#
# Table S2:
# Taxonomic resolution of unresolved nodes
#
# Summarize the taxonomic classifications assigned during manual curation to
# nodes that were not identified to species level, separately for the plant and
# animal guilds of pollination and seed-dispersal networks.
#
################################################################################

taxonomy_summary <- function(
    df,
    network_type,
    guild
) {
  
  if (!"Taxonomy" %in% names(df)) {
    
    stop(
      paste0(
        "A Taxonomy column is required for ",
        network_type,
        " - ",
        guild
      )
    )
  }
  
  df %>%
    mutate(
      Taxonomy =
        str_squish(
          as.character(
            Taxonomy
          )
        )
    ) %>%
    
    filter(
      !is.na(Taxonomy),
      Taxonomy != ""
    ) %>%
    
    count(
      Taxonomy,
      name =
        "Count"
    ) %>%
    
    mutate(
      Percentage =
        round(
          100 *
            Count /
            sum(Count),
          2
        ),
      
      Network_type =
        network_type,
      
      Guild =
        guild
    ) %>%
    
    select(
      Network_type,
      Guild,
      Taxonomy,
      Count,
      Percentage
    )
}

all_summaries <- bind_rows(
  
  taxonomy_summary(
    unresolved_pollination_row,
    "Pollination",
    "Plant"
  ),
  
  taxonomy_summary(
    unresolved_pollination_column,
    "Pollination",
    "Animal"
  ),
  
  taxonomy_summary(
    unresolved_seed_row,
    "Seed-dispersal",
    "Plant"
  ),
  
  taxonomy_summary(
    unresolved_seed_column,
    "Seed-dispersal",
    "Animal"
  )
) %>%
  
  arrange(
    Network_type,
    Guild,
    desc(
      Percentage
    )
  )

all_summaries

# Save table
write_csv(all_summaries, file.path(output_dir, "table_s2.csv"))
