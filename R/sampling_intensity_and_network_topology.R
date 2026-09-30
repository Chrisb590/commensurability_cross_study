library(dplyr)
library(ggplot2)
library(ggpmisc)
library(here)
library(patchwork)
library(readr)
library(stringr)

################################################################################
#
# Data
#
################################################################################

# Directory where all figures and tables are written (created if missing)
output_dir <- here("output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

network_directory <- here("data","networks")

metadata <- read_csv(
  here("data","general_network_information.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    ID = str_trim(ID),
    TYPE = TYPE %>%
      str_trim() %>%
      str_to_lower() %>%
      str_replace_all(" ", "-") %>%
      str_to_title()
  )

# The empirical topological indices were previously calculated for the
# Patefield null-model analysis.
topological_indices <- read_csv(
  here("data","all_results_for_patefield_randomization.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    ID = str_trim(ID),
    H2 = as.numeric(H2),
    weighted_NODF = as.numeric(weighted_NODF),
    weighted_modularity = as.numeric(DIRTMod)
  )

# Read an interaction matrix and ensure that all entries are numeric.
read_network_csv <- function(path) {
  network <- read.csv(
    path,
    row.names = 1,
    check.names = FALSE
  ) %>%
    as.matrix()
  
  suppressWarnings(
    storage.mode(network) <- "double"
  )
  
  network[is.na(network)] <- 0
  
  if (is.null(rownames(network))) {
    rownames(network) <- paste0("r", seq_len(nrow(network)))
  }
  
  if (is.null(colnames(network))) {
    colnames(network) <- paste0("c", seq_len(ncol(network)))
  }
  
  network
}

# Calculate the quantities required for sampling intensity and network size.
network_files <- list.files(
  network_directory,
  pattern = "\\.csv$",
  full.names = TRUE
)

network_dimensions <- lapply(
  network_files,
  function(path) {
    network <- read_network_csv(path)
    
    data.frame(
      ID = sub("\\.csv$", "", basename(path)),
      Plants = nrow(network),
      Animals = ncol(network),
      E = sum(network),
      stringsAsFactors = FALSE
    )
  }
) %>%
  bind_rows()

# Combine network dimensions, empirical topological indices, and network type.
# Sampling intensity follows Eq. 1:
# sqrt(total edge weight / [number of plant nodes x number of animal nodes]).
results <- network_dimensions %>%
  mutate(
    ID = str_trim(ID)
  ) %>%
  left_join(
    topological_indices,
    by = "ID"
  ) %>%
  left_join(
    metadata,
    by = "ID"
  ) %>%
  mutate(
    TYPE = factor(
      TYPE,
      levels = c(
        "Pollination",
        "Seed-Dispersal"
      )
    ),
    NetworkSize = Plants + Animals,
    SamplingIntensity = sqrt(E / (Plants * Animals)),
    SI_log = log(SamplingIntensity),
    N_log = log(NetworkSize)
  ) %>%
  filter(
    is.finite(SI_log),
    is.finite(N_log)
  )

################################################################################
#
# Main manuscript analysis:
# Figure 6 - Sampling intensity, network size, and network topology
#
# Panel A compares sampling intensity with network size and overlays the
# expected N^(-1/2) scaling from Eq. 2. Panels B-D compare sampling intensity
# with specialization (H2'), nestedness (wNODF), and modularity (DIRTLPAwb+),
# respectively, across the 279 quantitative interaction networks.
#
################################################################################

################################################################################
#
# Figure formatting
#
################################################################################

network_type_breaks <- c(
  "Pollination",
  "Seed-Dispersal"
)

network_type_labels <- c(
  "Pollination" = "Pollination",
  "Seed-Dispersal" = "Seed-dispersal"
)

network_fill_scale <- scale_fill_manual(
  values = c(
    "Pollination" = "#fdbb84",
    "Seed-Dispersal" = "#a1dab4"
  ),
  breaks = network_type_breaks,
  labels = network_type_labels,
  drop = FALSE
)

network_colour_scale <- scale_color_manual(
  values = c(
    "Pollination" = "#d95f02",
    "Seed-Dispersal" = "#1b9e77"
  ),
  breaks = network_type_breaks,
  guide = "none",
  drop = FALSE
)

theme_figure6 <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "grey",
      fill = NA,
      linewidth = 0.5
    ),
    legend.position = "top",
    legend.title = element_blank(),
    plot.margin = margin(6, 6, 6, 6),
    aspect.ratio = 1
  )

x_label_log_sampling_intensity <- expression("ln(Sampling intensity)")
x_label_log_network_size <- expression("ln(Network size)")

# Create a scatterplot with one linear model fitted across both network types.
scatter_lm_onefit <- function(
    data,
    x_variable,
    y_variable,
    y_label,
    x_label = x_label_log_sampling_intensity,
    text_size = 4,
    x_position = "right",
    y_position = 0.9
) {
  ggplot(
    data,
    aes_string(
      x = x_variable,
      y = y_variable
    )
  ) +
    geom_point(
      aes(
        fill = TYPE,
        color = TYPE
      ),
      shape = 21,
      size = 3,
      na.rm = TRUE
    ) +
    geom_smooth(
      aes(group = 1),
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      color = "black"
    ) +
    stat_poly_eq(
      aes_string(
        x = x_variable,
        y = y_variable,
        label = "paste(..rr.label.., ..p.value.label.., sep = ' ~~~ ')"
      ),
      formula = y ~ x,
      parse = TRUE,
      label.x.npc = x_position,
      label.y.npc = y_position,
      size = text_size,
      inherit.aes = FALSE,
      color = "black"
    ) +
    labs(
      x = x_label,
      y = y_label,
      fill = NULL
    ) +
    network_fill_scale +
    network_colour_scale +
    guides(
      fill = guide_legend(
        override.aes = list(
          color = c(
            "#d95f02",
            "#1b9e77"
          )
        )
      )
    ) +
    theme_figure6
}

################################################################################
#
# Figure 6A:
# Sampling intensity versus network size
#
################################################################################

scaling_model <- lm(
  SI_log ~ offset(-0.5 * N_log),
  data = results
)

k_hat <- as.numeric(
  exp(
    coef(scaling_model)[1]
  )
)

predicted_log_sampling_intensity <- predict(
  scaling_model
)

scaling_r_squared <- cor(
  results$SI_log,
  predicted_log_sampling_intensity
)^2

scaling_curve <- data.frame(
  NetworkSize = seq(
    min(results$NetworkSize, na.rm = TRUE),
    max(results$NetworkSize, na.rm = TRUE),
    length.out = 400
  )
) %>%
  mutate(
    SamplingIntensity = k_hat / sqrt(NetworkSize),
    N_log = log(NetworkSize),
    SI_log = log(SamplingIntensity)
  )

k_label <- format(
  round(k_hat, 2),
  nsmall = 2
)

r_squared_label <- format(
  round(scaling_r_squared, 2),
  nsmall = 2
)

scaling_label <- bquote(
  atop(
    "ln(Sampling intensity)" ==
      ln(k) -
      0.5 %.%
      "ln(Network size)",
    k ==
      .(k_label) ~~
      italic(R)^2 ==
      .(r_squared_label)
  )
) %>%
  as.expression() %>%
  as.character()

p_a <- ggplot(
  results,
  aes(
    x = N_log,
    y = SI_log
  )
) +
  geom_point(
    aes(
      fill = TYPE,
      color = TYPE
    ),
    shape = 21,
    size = 3,
    na.rm = TRUE
  ) +
  geom_line(
    data = scaling_curve,
    aes(
      x = N_log,
      y = SI_log
    ),
    color = "black",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = scaling_label,
    parse = TRUE,
    hjust = 1.05,
    vjust = 1.8,
    size = 4
  ) +
  labs(
    x = x_label_log_network_size,
    y = x_label_log_sampling_intensity,
    fill = NULL
  ) +
  network_fill_scale +
  network_colour_scale +
  guides(
    fill = guide_legend(
      override.aes = list(
        color = c(
          "#d95f02",
          "#1b9e77"
        )
      )
    )
  ) +
  theme_figure6 +
  coord_cartesian(
    clip = "off"
  )

################################################################################
#
# Figure 6B-D:
# Topological indices versus sampling intensity
#
################################################################################

p_b <- scatter_lm_onefit(
  results,
  "SI_log",
  "H2",
  y_label = bquote(
    "Specialization (" *
      H[2] *
      "')"
  )
)

p_c <- scatter_lm_onefit(
  results,
  "SI_log",
  "weighted_NODF",
  y_label = "Weighted nestedness (wNODF)"
)

p_d <- scatter_lm_onefit(
  results,
  "SI_log",
  "weighted_modularity",
  y_label = "Weighted modularity (DIRTLPAwb+)"
)

################################################################################
#
# Assemble and save Figure 6
#
################################################################################

figure6 <- (
  p_a |
    p_b
) /
  (
    p_c |
      p_d
  ) +
  plot_layout(
    guides = "collect",
    axes = "collect"
  ) +
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.background = element_rect(
      colour = "grey",
      fill = NA,
      linewidth = 0.5
    ),
    plot.tag = element_text(
      face = "bold"
    )
  )

print(
  figure6
)

ggsave(
  file.path(output_dir,"figure_6.pdf"),
  plot = figure6,
  width = 9.2,
  height = 9.2,
  units = "in"
)
