library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
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

metadata <- read_csv(here("data", "general_network_information.csv"))

# Raw network directory
network_dir <- here("data","networks")

# Directories where beta-diversity CSVs are stored
row_dir <- here("data","networks", "row_names")
col_dir <- here("data","networks", "column_names")

dir.create(row_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(col_dir, recursive = TRUE, showWarnings = FALSE)

# Interaction types
interaction_types <- c(
  "Seed_Dispersal",
  "Pollination"
)

# Names for plotting
display_map <- c(
  "Seed_Dispersal" = "Seed-dispersal",
  "Pollination"    = "Pollination"
)

# Standardize interaction-type label in metadata
metadata <- metadata %>%
  mutate(
    Interaction_key = case_when(
      str_detect(
        TYPE,
        regex("pollination", ignore_case = TRUE)
      ) ~ "Pollination",
      
      str_detect(
        TYPE,
        regex("seed", ignore_case = TRUE)
      ) ~ "Seed_Dispersal",
      
      TRUE ~ NA_character_
    )
  )

################################################################################
#
# Main manuscript analyses:
# Beta-diversity among quantitative interaction networks
#
# Calculate beta-diversity using only nodes identified to species level,
# compare beta-diversity within and between publications, generate Figure 3,
# and calculate the exact median beta-diversity values reported in the main
# manuscript.
#
################################################################################

################################################################################
#
# Read and clean curated lists of unresolved node names
#
# This function reads a CSV containing node names previously identified as
# unresolved to species level, extracts the species-name column, removes
# missing or blank entries, standardizes whitespace, and returns unique names.
#
################################################################################

read_unresolved_names <- function(path) {
  
  x <- read_csv(
    path,
    show_col_types = FALSE
  )
  
  species_col <- intersect(
    c("species", "Species"),
    names(x)
  )
  
  if (length(species_col) == 0) {
    
    stop(
      paste0(
        "No species column found in ",
        path
      )
    )
  }
  
  out <- as.character(
    x[[species_col[1]]]
  )
  
  out <- str_squish(out)
  
  out <- out[
    !is.na(out) &
      nzchar(out)
  ]
  
  unique(out)
}

unresolved <- list(
  
  Seed_Dispersal_row =
    read_unresolved_names(
      here(row_dir,
        "seed_dispersal_unique_unidentified_row_species_names.csv"
      )
    ),
  
  Seed_Dispersal_column =
    read_unresolved_names(
      here(col_dir,
        "seed_dispersal_unique_unidentified_column_species_names.csv"
      )
    ),
  
  Pollination_row =
    read_unresolved_names(
      here(row_dir,
        "pollination_unique_unidentified_row_species_names.csv"
      )
    ),
  
  Pollination_column =
    read_unresolved_names(
      here(col_dir,
        "pollination_unique_unidentified_column_species_names.csv"
      )
    )
)

################################################################################
#
# Extract species-level nodes from a single network
#
# For a given network and guild, read the node names, clean the labels, and
# remove any nodes listed in the corresponding curated unresolved-taxa list.
# Matching to unresolved taxa is case-insensitive.
#
################################################################################

extract_species_only <- function(
    id,
    interaction_key,
    data_type
) {
  
  path <- file.path(
    network_dir,
    paste0(id, ".csv")
  )
  
  if (!file.exists(path)) {
    
    warning(
      "Network file not found: ",
      path
    )
    
    return(character())
  }
  
  x <- read_csv(
    path,
    show_col_types = FALSE,
    name_repair = "minimal"
  )
  
  # Rows = plant guild
  if (data_type == "row") {
    
    species <- as.character(
      x[[1]]
    )
  }
  
  # Columns = animal guild
  if (data_type == "column") {
    
    species <- colnames(x)[-1]
  }
  
  # Clean names
  species <- str_squish(species)
  
  species <- species[
    !is.na(species) &
      nzchar(species)
  ]
  
  # Select appropriate unresolved-node list
  lookup_name <- paste(
    interaction_key,
    data_type,
    sep = "_"
  )
  
  unresolved_here <-
    unresolved[[lookup_name]]
  
  # Remove unresolved nodes
  # Matching is case-insensitive but otherwise exact
  species <- species[
    !tolower(species) %in%
      tolower(unresolved_here)
  ]
  
  unique(species)
}

################################################################################
#
# Calculate pairwise Jaccard beta-diversity from supplied node sets
#
# This function converts a list of node sets into a presence/absence matrix
# and calculates pairwise Jaccard dissimilarity among networks.
#
################################################################################

jaccard_matrix_from_species_sets <- function(
    species_sets
) {
  
  ids <- names(species_sets)
  
  all_species <- sort(
    unique(
      unlist(
        species_sets,
        use.names = FALSE
      )
    )
  )
  
  if (length(all_species) == 0) {
    
    stop(
      "No species-level taxa remain."
    )
  }
  
  presence <- matrix(
    0L,
    nrow = length(species_sets),
    ncol = length(all_species),
    dimnames = list(
      ids,
      all_species
    )
  )
  
  for (i in seq_along(species_sets)) {
    
    spp <- species_sets[[i]]
    
    if (length(spp) > 0) {
      
      presence[
        i,
        match(
          spp,
          all_species
        )
      ] <- 1L
    }
  }
  
  # Number of shared species
  intersection_n <-
    tcrossprod(presence)
  
  # Species richness of each network
  species_n <-
    rowSums(presence)
  
  # Number of species in union
  union_n <- outer(
    species_n,
    species_n,
    "+"
  ) - intersection_n
  
  # Jaccard dissimilarity
  jaccard <-
    1 -
    intersection_n /
    union_n
  
  diag(jaccard) <- 0
  
  dimnames(jaccard) <- list(
    ids,
    ids
  )
  
  jaccard
}

################################################################################
#
# Build species-only beta-diversity matrices
#
# For each interaction type and guild, extract only nodes identified to species,
# remove networks with no species-level taxa remaining, and calculate pairwise
# Jaccard beta-diversity using the filtered node sets.
#
################################################################################

build_species_beta <- function(
    interaction_key,
    data_type
) {
  
  ids <- metadata %>%
    filter(
      Interaction_key ==
        interaction_key
    ) %>%
    pull(ID) %>%
    unique()
  
  species_sets <- setNames(
    
    lapply(
      ids,
      function(id) {
        
        extract_species_only(
          id = id,
          interaction_key = interaction_key,
          data_type = data_type
        )
      }
    ),
    
    ids
  )
  
  # Remove a network only if no species-level taxa
  # remain within this guild
  keep <-
    lengths(species_sets) > 0
  
  
  if (any(!keep)) {
    
    warning(
      sum(!keep),
      " network(s) removed because no species-level taxa remained for ",
      interaction_key,
      " / ",
      data_type
    )
  }
  
  species_sets <-
    species_sets[keep]
  
  beta_mat <-
    jaccard_matrix_from_species_sets(
      species_sets
    )
  
  # Preserve naming convention used by original
  # beta-diversity workflow
  suffix <- ifelse(
    data_type == "row",
    "_row_names",
    "_column_names"
  )
  
  rownames(beta_mat) <-
    paste0(
      rownames(beta_mat),
      suffix
    )
  
  colnames(beta_mat) <-
    paste0(
      colnames(beta_mat),
      suffix
    )
  
  beta_mat
}

################################################################################
#
# Build beta-diversity matrices using only nodes identified to species level
#
# Four Jaccard beta-diversity matrices are generated:
#   1. Pollination - plant guild
#   2. Pollination - pollinator guild
#   3. Seed dispersal - plant guild
#   4. Seed dispersal - seed-disperser guild
#
################################################################################

beta_species_mats <- list()

for (intr in interaction_types) {
  
  beta_species_mats[[
    paste0(
      intr,
      "_row"
    )
  ]] <- build_species_beta(
    interaction_key = intr,
    data_type = "row"
  )
  
  beta_species_mats[[
    paste0(
      intr,
      "_column"
    )
  ]] <- build_species_beta(
    interaction_key = intr,
    data_type = "column"
  )
}

################################################################################
#
# Convert beta-diversity matrices to pairwise network comparisons
#
# Extract the unique pairwise Jaccard dissimilarities from each beta-diversity
# matrix and join publication metadata to each network pair. Pairs are then
# classified according to whether the networks come from the same publication
# or represent publications containing a single network.
#
################################################################################

get_pairwise_beta <- function(
    beta_mat,
    interaction_key,
    data_type
) {
  
  beta_mat <- as.matrix(beta_mat)
  
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
      ID1_clean = sub(
        "(_row_names|_column_names)$",
        "",
        ID1
      ),
      
      ID2_clean = sub(
        "(_row_names|_column_names)$",
        "",
        ID2
      )
    ) %>%
    
    left_join(
      metadata %>%
        select(
          ID,
          Publication
        ),
      by = c(
        "ID1_clean" = "ID"
      )
    ) %>%
    
    rename(
      Pub1 = Publication
    ) %>%
    
    left_join(
      metadata %>%
        select(
          ID,
          Publication
        ),
      by = c(
        "ID2_clean" = "ID"
      )
    ) %>%
    
    rename(
      Pub2 = Publication
    ) %>%
    
    mutate(
      
      Group = case_when(
        
        Pub1 ==
          "One network per publication" &
          Pub2 ==
          "One network per publication" ~
          "One network per publication",
        
        Pub1 !=
          "One network per publication" &
          Pub2 !=
          "One network per publication" &
          Pub1 == Pub2 ~
          "Within Publication",
        
        TRUE ~ NA_character_
      ),
      
      Interaction_key =
        interaction_key,
      
      DataType =
        data_type
    ) %>%
    
    filter(
      !is.na(Group)
    ) %>%
    
    select(
      ID1,
      ID2,
      beta,
      Pub1,
      Pub2,
      Group,
      Interaction_key,
      DataType
    )
}

################################################################################
#
# Convert species-only beta-diversity matrices to pairwise comparisons
#
################################################################################

all_pairwise <- list()

for (intr in interaction_types) {
  
  all_pairwise[[
    paste0(
      intr,
      "_row"
    )
  ]] <- get_pairwise_beta(
    beta_species_mats[[
      paste0(
        intr,
        "_row"
      )
    ]],
    intr,
    "row"
  )
  
  all_pairwise[[
    paste0(
      intr,
      "_column"
    )
  ]] <- get_pairwise_beta(
    beta_species_mats[[
      paste0(
        intr,
        "_column"
      )
    ]],
    intr,
    "column"
  )
}

################################################################################
#
# Prepare beta-diversity plot
#
################################################################################

combined_df <- bind_rows(
  all_pairwise
) %>%
  
  mutate(
    
    Interaction = factor(
      Interaction_key,
      levels = names(display_map),
      labels = unname(display_map)
    ),
    
    Combo = case_when(
      
      Group ==
        "One network per publication" &
        DataType == "row" ~
        "Row - One per pub",
      
      Group ==
        "One network per publication" &
        DataType == "column" ~
        "Column - One per pub",
      
      Group ==
        "Within Publication" &
        DataType == "row" ~
        "Row - Within pub",
      
      Group ==
        "Within Publication" &
        DataType == "column" ~
        "Column - Within pub",
      
      TRUE ~ NA_character_
    )
  ) %>%
  
  filter(
    !is.na(Combo)
  ) %>%
  
  mutate(
    ComboInteraction = interaction(
      Combo,
      Interaction_key,
      drop = TRUE
    )
  )

# Ensure all "One network per publication"
# levels come first
levels_one <- combined_df %>%
  
  filter(
    Group ==
      "One network per publication"
  ) %>%
  
  distinct(
    ComboInteraction
  ) %>%
  
  pull()

# Multiple networks per publication come second
levels_within <- combined_df %>%
  
  filter(
    Group ==
      "Within Publication"
  ) %>%
  
  distinct(
    ComboInteraction
  ) %>%
  
  pull()

combo_levels <- c(
  levels_one,
  levels_within
)

combined_df <- combined_df %>%
  
  mutate(
    ComboInteraction = factor(
      ComboInteraction,
      levels = combo_levels
    )
  )

# Separator between publication groups
separator_position <-
  length(levels_one) + 0.5

# Simplified x-axis labels
combo_labels <- combined_df %>%
  
  distinct(
    ComboInteraction,
    DataType
  ) %>%
  
  mutate(
    ComboLabel = ifelse(
      DataType == "row",
      "Plant\nguild",
      "Animal\nguild"
    )
  ) %>%
  
  select(
    ComboInteraction,
    ComboLabel
  ) %>%
  
  deframe()

# Number of pairwise comparisons above each box
combo_counts <- combined_df %>%
  
  count(
    ComboInteraction
  ) %>%
  
  mutate(
    label = paste0(
      "'(' * italic(n) * '=' * ",
      n,
      " * ')'"
    )
  )

################################################################################
#
# Beta-diversity plot (Figure 3)
#
################################################################################

fill_vals <- c(
  "Seed_Dispersal" = "#a1dab4",
  "Pollination" = "#fdbb84"
)

col_vals <- c(
  "Seed_Dispersal" = "#1b9e77",
  "Pollination" = "#d95f02"
)

p_beta_dot <- ggplot(
  combined_df,
  aes(
    x = ComboInteraction,
    y = beta
  )
) +
  
  geom_jitter(
    aes(
      color = Interaction_key
    ),
    width = 0.2,
    size = 1.8,
    alpha = 0.7
  ) +
  
  geom_boxplot(
    aes(
      fill = Interaction_key
    ),
    width = 0.5,
    alpha = 0.4,
    outlier.shape = NA,
    color = "black",
    size = 0.6
  ) +
  
  geom_text(
    data = combo_counts,
    aes(
      x = ComboInteraction,
      y =
        max(
          combined_df$beta,
          na.rm = TRUE
        ) + 0.05,
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
  
  scale_x_discrete(
    labels = combo_labels
  ) +
  
  scale_y_continuous(
    breaks = seq(
      0,
      1,
      by = 0.2
    )
  ) +
  
  labs(
    y = "Beta-diversity (Jaccard)",
    x = NULL
  ) +
  
  theme_bw() +
  
  theme(
    
    axis.text.x =
      element_text(
        size = 16,
        color = "black",
        lineheight = 0.9
      ),
    
    axis.text.y =
      element_text(
        size = 18,
        color = "black"
      ),
    
    axis.title.y =
      element_text(
        size = 18,
        face = "bold"
      ),
    
    legend.position =
      "top",
    
    legend.title =
      element_blank(),
    
    legend.text =
      element_text(
        size = 18
      ),
    
    panel.border =
      element_rect(
        colour = "grey",
        fill = NA,
        size = 0.4
      )
  ) +
  
  geom_vline(
    xintercept =
      separator_position,
    linetype =
      "dashed",
    size =
      0.7
  ) +
  
  annotate(
    "text",
    x =
      separator_position - 0.8,
    y =
      max(
        combined_df$beta,
        na.rm = TRUE
      ) + 0.12,
    label =
      "One network per publication",
    hjust =
      1,
    size =
      5.5,
    fontface =
      "bold"
  ) +
  
  annotate(
    "text",
    x =
      separator_position + 0.8,
    y =
      max(
        combined_df$beta,
        na.rm = TRUE
      ) + 0.12,
    label =
      "Multiple networks per publication",
    hjust =
      0,
    size =
      5.5,
    fontface =
      "bold"
  )

print(
  p_beta_dot
)

# Save final plot as PDF
ggsave(
  file.path(output_dir, "figure_3.pdf"),
  device = cairo_pdf,
  width = 12,
  height = 5,
  units = "in"
)

################################################################################
#
# Calculate exact median beta-diversity values reported in Figure 3
#
# For each interaction type and guild, calculate the median pairwise Jaccard
# beta-diversity corresponding to the distributions shown in Figure 3.
# Also report the number of pairwise comparisons and unique networks
# contributing to each group.
#
################################################################################

group_map <- c(
  "Within Publication" =
    "Within publication",
  
  "One network per publication" =
    "Between publications"
)

combined_df_processed <- combined_df %>%
  
  mutate(
    
    Set = recode(
      Group,
      !!!group_map
    ),
    
    Guild = if_else(
      DataType == "row",
      "Plant guild",
      "Animal guild"
    ),
    
    Interaction =
      Interaction_key
  ) %>%
  
  filter(
    !is.na(Set)
  )

# Median beta-diversity within/between
beta_medians_tbl <- combined_df_processed %>%
  
  group_by(
    Set,
    Interaction,
    Guild
  ) %>%
  
  summarise(
    
    median_beta =
      median(
        beta,
        na.rm = TRUE
      ),
    
    n_pairs =
      sum(
        !is.na(beta)
      ),
    
    .groups =
      "drop"
  )

# Number of networks within/between
network_counts_tbl <- combined_df_processed %>%
  
  select(
    Set,
    Interaction,
    ID1,
    ID2
  ) %>%
  
  pivot_longer(
    cols = c(
      ID1,
      ID2
    ),
    values_to =
      "ID_with_suffix"
  ) %>%
  
  mutate(
    ID = sub(
      "(_row_names|_column_names)$",
      "",
      ID_with_suffix
    )
  ) %>%
  
  distinct(
    Set,
    Interaction,
    ID
  ) %>%
  
  count(
    Set,
    Interaction,
    name =
      "n_networks"
  )

final_within_between_tbl <-
  beta_medians_tbl %>%
  
  left_join(
    network_counts_tbl,
    by = c(
      "Set",
      "Interaction"
    )
  ) %>%
  
  mutate(
    
    median_beta =
      round(
        median_beta,
        2
      ),
    
    n_pairs =
      as.integer(
        n_pairs
      ),
    
    n_networks =
      as.integer(
        n_networks
      )
  ) %>%
  
  arrange(
    Set,
    Guild,
    Interaction
  )

print(
  final_within_between_tbl
)

################################################################################
#
# Section S2:
# Geographic distance and publication identity as predictors of beta-diversity
#
# Evaluate whether pairwise beta-diversity among networks is associated with
# geographic distance and/or whether networks originate from the same
# publication. Analyses are conducted separately for the plant and animal
# guilds of pollination and seed-dispersal networks, using only nodes
# identified to species level.
#
# For each interaction type and guild, three multiple regressions on distance
# matrices (MRMs) are fitted:
#
#   1. Geographic distance + publication identity
#   2. Geographic distance only
#   3. Publication identity only
#
################################################################################


################################################################################
#
# Prepare metadata for MRM analyses
#
# Convert latitude and longitude to numeric values for calculating geographic
# distances. Assign each network a publication identity for determining whether
# pairs of networks originate from the same or different publications.
#
# Networks labelled "One network per publication" are assigned unique
# publication identities because each represents a different publication.
#
################################################################################

metadata_mrm <- metadata %>%
  
  mutate(
    
    Latitude_mrm =
      as.numeric(
        Latitude
      ),
    
    Longitude_mrm =
      as.numeric(
        Longitude
      ),
    
    Publication_ID =
      case_when(
        
        str_detect(
          Publication,
          regex(
            "^\\s*one\\s+network\\s+per\\s+publication\\s*$",
            ignore_case = TRUE
          )
        ) ~
          paste0(
            "Singleton_",
            ID
          ),
        
        is.na(Publication) |
          Publication == "" ~
          paste0(
            "Unknown_",
            ID
          ),
        
        TRUE ~
          Publication
      )
  )

################################################################################
#
# Calculate pairwise geographic distance among networks
#
# Calculate the geographic distance in kilometres between each pair of network
# locations from their latitude and longitude coordinates. Geographic distance
# is used as one of the pairwise predictors in the MRM analyses.
#
# Distances are calculated along the curved surface of the Earth using the
# Haversine formula.
#
################################################################################

haversine_distance <- function(
    longitude,
    latitude,
    ids
) {
  
  n <- length(ids)
  
  output <- matrix(
    0,
    nrow = n,
    ncol = n,
    dimnames = list(
      ids,
      ids
    )
  )
  
  # Mean radius of the Earth in kilometres, used to convert angular
  # separation between coordinates into geographic distance
  earth_radius_km <-
    6371.0088
  
  for (i in seq_len(n - 1)) {
    
    for (j in (i + 1):n) {
      
      # Convert latitude from degrees to radians
      lat1 <-
        latitude[i] *
        pi / 180
      
      lat2 <-
        latitude[j] *
        pi / 180
      
      # Difference in latitude and longitude between locations
      delta_lat <-
        (
          latitude[j] -
            latitude[i]
        ) *
        pi / 180
      
      delta_lon <-
        (
          longitude[j] -
            longitude[i]
        ) *
        pi / 180
      
      # Haversine calculation
      a <-
        sin(
          delta_lat / 2
        )^2 +
        cos(lat1) *
        cos(lat2) *
        sin(
          delta_lon / 2
        )^2
      
      # Constrain the intermediate value to [0, 1] to avoid numerical
      # rounding errors before calculating angular distance
      a <- min(
        1,
        max(
          0,
          a
        )
      )
      
      angular_distance <-
        2 *
        atan2(
          sqrt(a),
          sqrt(1 - a)
        )
      
      geographic_distance <-
        earth_radius_km *
        angular_distance
      
      
      output[i, j] <-
        geographic_distance
      
      output[j, i] <-
        geographic_distance
    }
  }
  
  as.dist(output)
}

################################################################################
#
# Create pairwise publication-identity predictor
#
# For each pair of networks, indicate whether they originate from the same
# publication (0) or different publications (1). This binary pairwise variable
# is used as the publication-identity predictor in the Section S2 MRM analyses.
#
################################################################################

publication_difference <- function(
    publication_ids,
    ids
) {
  
  output <- outer(
    publication_ids,
    publication_ids,
    FUN = "!="
  ) * 1
  
  dimnames(output) <- list(
    ids,
    ids
  )
  
  as.dist(output)
}

################################################################################
#
# Fit MRM models for one interaction type and guild
#
# For a specified interaction type and guild:
#
#   1. Match the beta-diversity matrix with network metadata.
#   2. Calculate pairwise geographic distance.
#   3. Create the pairwise publication-identity predictor.
#   4. Fit the combined, geographic-only, and publication-only MRM models.
#
# Beta-diversity matrices include only nodes identified to species level.
#
################################################################################

run_mrm <- function(
    interaction_key,
    data_type,
    nperm = 9999
) {
  
  ##########################################################################
  # Select beta-diversity matrix
  ##########################################################################
  
  beta_mat <-
    beta_species_mats[[
      paste0(
        interaction_key,
        "_",
        data_type
      )
    ]]
  
  # Remove row/column suffixes to recover network IDs
  beta_ids <- sub(
    "(_row_names|_column_names)$",
    "",
    rownames(beta_mat)
  )
  
  ##########################################################################
  # Match network metadata to beta-diversity matrix
  ##########################################################################
  
  # Retain networks in the beta-diversity matrix with valid geographic coordinates
  dat <- metadata_mrm %>%
    filter(
      ID %in% beta_ids,
      Interaction_key == interaction_key,
      is.finite(Latitude_mrm),
      is.finite(Longitude_mrm)
    )
  
  # Preserve the network ordering used in the beta-diversity matrix
  valid_ids <- beta_ids[
    beta_ids %in% dat$ID
  ]
  
  dat <- dat[
    match(
      valid_ids,
      dat$ID
    ),
  ]
  
  ##########################################################################
  # Prepare pairwise beta-diversity response
  ##########################################################################
  
  # Identify the corresponding row or column names in the beta-diversity
  # matrix
  
  suffix <- ifelse(
    data_type == "row",
    "_row_names",
    "_column_names"
  )
  
  matrix_ids <- paste0(
    valid_ids,
    suffix
  )
  
  # Retain networks with complete metadata
  beta_sub <- beta_mat[
    matrix_ids,
    matrix_ids,
    drop = FALSE
  ]
  
  # Replace suffixed names with clean network IDs
  dimnames(beta_sub) <- list(
    valid_ids,
    valid_ids
  )
  
  # Convert beta-diversity matrix to a distance object for MRM
  beta_dist <-
    as.dist(
      beta_sub
    )
  
  ##########################################################################
  # Prepare pairwise geographic-distance predictor
  ##########################################################################
  
  geo_dist_km <- haversine_distance(
    longitude =
      dat$Longitude_mrm,
    latitude =
      dat$Latitude_mrm,
    ids =
      valid_ids
  )
  
  # Log-transform geographic distance to reduce the influence of very large
  # distances while retaining comparisons with zero geographic distance
  log_geo_dist <- as.dist(
    log1p(
      as.matrix(
        geo_dist_km
      )
    )
  )
  
  ##########################################################################
  # Prepare pairwise publication-identity predictor
  ##########################################################################
  
  publication_diff <- publication_difference(
    publication_ids =
      dat$Publication_ID,
    ids =
      valid_ids
  )
  
  ##########################################################################
  # Fit MRM models
  ##########################################################################
  
  # Model 1: geographic distance + publication identity
  fit_combined <- ecodist::MRM(
    beta_dist ~
      log_geo_dist +
      publication_diff,
    nperm = nperm
  )
  
  # Model 2: geographic distance only
  fit_geography <- ecodist::MRM(
    beta_dist ~
      log_geo_dist,
    nperm = nperm
  )
  
  # Model 3: publication identity only
  fit_publication <- ecodist::MRM(
    beta_dist ~
      publication_diff,
    nperm = nperm
  )
  
  ##########################################################################
  # Return fitted models and variables used in the analyses
  ##########################################################################
  
  list(
    
    fit_combined =
      fit_combined,
    
    fit_geography =
      fit_geography,
    
    fit_publication =
      fit_publication,
    
    metadata =
      dat,
    
    beta_dist =
      beta_dist,
    
    geo_dist_km =
      geo_dist_km,
    
    log_geo_dist =
      log_geo_dist,
    
    publication_diff =
      publication_diff
  )
}

################################################################################
#
# Run MRM analyses for each interaction type and guild
#
# The three MRM models are fitted separately for:
#
#   1. Pollination - plant guild
#   2. Pollination - pollinator guild
#   3. Seed dispersal - plant guild
#   4. Seed dispersal - seed-disperser guild
#
################################################################################

set.seed(
  12345
)

mrm_results <- list(
  
  Pollination_plant =
    run_mrm(
      interaction_key =
        "Pollination",
      data_type =
        "row"
    ),
  
  Pollination_pollinator =
    run_mrm(
      interaction_key =
        "Pollination",
      data_type =
        "column"
    ),
  
  Seed_Dispersal_plant =
    run_mrm(
      interaction_key =
        "Seed_Dispersal",
      data_type =
        "row"
    ),
  
  Seed_Dispersal_seed_disperser =
    run_mrm(
      interaction_key =
        "Seed_Dispersal",
      data_type =
        "column"
    )
)

################################################################################
#
# Summarize MRM results (Table S1)
#
# For each interaction type and guild, report the number of networks and the
# variance in beta-diversity explained by geographic distance, publication
# identity, and both predictors together. Incremental R2 values represent the
# additional variance explained by each predictor after accounting for the
# other.
#
################################################################################

get_mrm_R2 <- function(fit) {
  
  as.numeric(
    fit$r.squared[1]
  )
}

extract_mrm_R2 <- function(
    object,
    interaction,
    guild
) {
  
  combined_R2 <-
    get_mrm_R2(
      object$fit_combined
    )
  
  geography_R2 <-
    get_mrm_R2(
      object$fit_geography
    )
  
  publication_R2 <-
    get_mrm_R2(
      object$fit_publication
    )
  
  tibble(
    
    Interaction =
      interaction,
    
    Guild =
      guild,
    
    n_networks =
      nrow(
        object$metadata
      ),
    
    Geographic_R2 =
      geography_R2,
    
    Publication_R2 =
      publication_R2,
    
    Combined_R2 =
      combined_R2,
    
    Geography_added_R2 =
      combined_R2 -
      publication_R2,
    
    Publication_added_R2 =
      combined_R2 -
      geography_R2
  )
}

mrm_R2_summary <- bind_rows(
  
  extract_mrm_R2(
    mrm_results$Pollination_plant,
    "Pollination",
    "Plant"
  ),
  
  extract_mrm_R2(
    mrm_results$Pollination_pollinator,
    "Pollination",
    "Pollinator"
  ),
  
  extract_mrm_R2(
    mrm_results$Seed_Dispersal_plant,
    "Seed dispersal",
    "Plant"
  ),
  
  extract_mrm_R2(
    mrm_results$Seed_Dispersal_seed_disperser,
    "Seed dispersal",
    "Seed disperser"
  )
) %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 2)
    )
  )

print(
  mrm_R2_summary,
  width = Inf
)

# Save table
write_csv(mrm_R2_summary, file.path(output_dir, "table_s1.csv"))

