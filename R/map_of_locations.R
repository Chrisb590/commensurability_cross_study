library(dplyr)
library(ggplot2)
library(readr)
library(here)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

################################################################################
#
# Data
#
################################################################################

# Directory where all figures and tables are written (created if missing)
output_dir <- here("output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

metadata <- read_csv(
  here("data", "general_network_information.csv"),
  show_col_types = FALSE
)

# Raw network directory
network_dir <- here("data", "networks")

# Identify networks included in the analysis
network_files <- list.files(
  network_dir,
  pattern = "\\.csv$",
  full.names = FALSE
)

network_ids <- sub(
  "\\.csv$",
  "",
  network_files
)

# Retain metadata for networks included in the analysis
metadata <- metadata %>%
  filter(
    ID %in% network_ids
  )

################################################################################
#
# Geographic distribution of quantitative interaction networks (Figure S1)
#
# Map the geographic locations of pollination and seed-dispersal networks.
# Networks sharing the same coordinates and interaction type are pooled,
# with point size indicating the number of networks at each location.
#
################################################################################

################################################################################
#
# Prepare network metadata
#
# Exclude host-parasite networks and standardize interaction-type labels.
#
################################################################################

filtered_metadata <- metadata %>%
  filter(
    TYPE != "Host-Parasite"
  ) %>%
  mutate(
    TYPE = if_else(
      TYPE == "Seed Dispersal",
      "Seed-dispersal",
      TYPE
    )
  )

################################################################################
#
# Count networks at each geographic location
#
# For each interaction type, count the number of networks sharing the
# same latitude and longitude.
#
################################################################################

location_counts <- filtered_metadata %>%
  group_by(
    Latitude,
    Longitude,
    TYPE
  ) %>%
  summarise(
    n_networks = n(),
    .groups = "drop"
  )

################################################################################
#
# Prepare world map
#
################################################################################

world <- ne_countries(
  scale = "medium",
  returnclass = "sf"
)

################################################################################
#
# Map formatting
#
################################################################################

fill_colors <- c(
  "Seed-dispersal" = "#a1dab4",
  "Pollination" = "#fdbb84"
)

border_colors <- c(
  "Seed-dispersal" = "#1b9e77",
  "Pollination" = "#d95f02"
)

################################################################################
#
# Geographic distribution of networks
#
################################################################################

p_map <- ggplot() +
  
  # Ocean background
  geom_rect(
    aes(
      xmin = -180,
      xmax = 180,
      ymin = -90,
      ymax = 90
    ),
    fill = "#deebf7",
    color = NA
  ) +
  
  # Land outlines
  geom_sf(
    data = world,
    fill = "white",
    color = "gray70",
    linewidth = 0.2
  ) +
  
  # Bounding box
  annotate(
    "rect",
    xmin = -180,
    xmax = 180,
    ymin = -90,
    ymax = 90,
    fill = NA,
    color = "black",
    linewidth = 0.7
  ) +
  
  # Network sampling locations
  geom_point(
    data = location_counts,
    aes(
      x = Longitude,
      y = Latitude,
      size = n_networks,
      fill = TYPE,
      color = TYPE
    ),
    shape = 21,
    stroke = 0.7,
    alpha = 0.85
  ) +
  
  # Point sizes
  scale_size_continuous(
    name = "Number of Networks",
    range = c(
      2,
      10
    ),
    limits = c(
      1,
      NA
    ),
    breaks = c(
      1,
      5,
      10,
      20,
      40
    )
  ) +
  
  # Interaction-type colours
  scale_fill_manual(
    name = "Interaction Type",
    values = fill_colors
  ) +
  
  scale_color_manual(
    name = "Interaction Type",
    values = border_colors
  ) +
  
  # Legend formatting
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,
        size = 5,
        stroke = 0.7
      ),
      title.position = "top"
    ),
    color = "none",
    size = guide_legend(
      order = 1
    )
  ) +
  
  # Coordinates and labels
  coord_sf(
    expand = FALSE
  ) +
  
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  
  # Theme
  theme_minimal(
    base_size = 14
  ) +
  
  theme(
    panel.background =
      element_rect(
        fill = NA,
        color = NA
      ),
    
    plot.background =
      element_rect(
        fill = "white",
        color = NA
      ),
    
    panel.grid.major =
      element_line(
        color = "white"
      ),
    
    panel.grid.minor =
      element_blank(),
    
    legend.position =
      "right",
    
    legend.title =
      element_text(
        face = "bold",
        size = 12
      ),
    
    legend.text =
      element_text(
        size = 11
      ),
    
    plot.title =
      element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      ),
    
    axis.title =
      element_text(
        size = 14
      ),
    
    axis.text =
      element_text(
        size = 12
      ),
    
    plot.margin =
      margin(
        5,
        5,
        5,
        5
      )
  )

print(
  p_map
)

# Save final plot as PDF
ggsave(
  file.path(output_dir, "figure_s1.pdf"),
  plot = p_map,
  device = "pdf",
  width = 12,
  height = 6,
  units = "in",
  dpi = 600,
  limitsize = FALSE
)
