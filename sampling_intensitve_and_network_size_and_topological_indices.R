library(ggplot2)
library(patchwork)
library(ggpmisc)
library(readr)
library(dplyr)
library(stringr)
library(here)

################################################################################
#
# Data
#
################################################################################

net_dir <- here("networks")

meta <- read_csv(here("general_network_information.csv"), show_col_types = FALSE) %>%
  mutate(
    ID = trimws(ID),
    TYPE = trimws(TYPE) %>%
      tolower() %>%
      str_replace_all(" ", "-") %>%
      str_to_title()
  )

# Reading network csv
read_network_csv <- function(path) {
  x <- read.csv(path, row.names = 1, check.names = FALSE)
  x <- as.matrix(x)
  suppressWarnings(storage.mode(x) <- "double")
  x[is.na(x)] <- 0
  if (is.null(rownames(x))) rownames(x) <- paste0("r", seq_len(nrow(x)))
  if (is.null(colnames(x))) colnames(x) <- paste0("c", seq_len(ncol(x)))
  x
}

# Load networks
csv_files <- list.files(net_dir, pattern = "\\.csv$", full.names = TRUE)

# Since H2', NODF, weighted modularity is already calculated in Patefield 
# randomization, I can just read them in
indices <- read_csv("all_results_for_patefield_randomization.csv", show_col_types = FALSE) %>%
  transmute(
    ID = str_trim(ID),
    H2 = as.numeric(H2),
    weighted_NODF = as.numeric(weighted_NODF),
    weighted_modularity = as.numeric(DIRTMod)  
  )

# Compute additional network indices
size_list <- lapply(csv_files, function(fp) {
  web_raw <- read_network_csv(fp)
  
  Plants  <- nrow(web_raw)
  Animals <- ncol(web_raw)
  
  L <- sum(web_raw > 0)  
  E <- sum(web_raw)      
  
  denom <- Plants * Animals
  SI <- sqrt(E / denom)
  C  <- L / denom
  
  data.frame(
    file = basename(fp),
    ID = sub("\\.csv$", "", basename(fp)),
    Plants = Plants,
    Animals = Animals,
    L = L,
    E = E,
    SamplingIntensity = SI,
    Connectance = C,
    stringsAsFactors = FALSE
  )
})

sizes <- do.call(rbind, size_list)

# Joining all indices for each network
results <- sizes %>%
  mutate(ID = str_trim(ID)) %>%
  left_join(indices, by = "ID")

# Join metadata and network indices
results <- results %>%
  select(-matches("^TYPE$"), everything()) %>%
  left_join(meta %>% select(ID, TYPE, Publication, Latitude), by = "ID") %>%
  mutate(
    TYPE = factor(TYPE, levels = c("Pollination","Seed-Dispersal")),
    NetworkSize = Plants + Animals,
    SI_log = log(SamplingIntensity),
    N_log  = log(NetworkSize)
  ) %>%
  filter(is.finite(SI_log), is.finite(N_log))

################################################################################
#
# Figure 6
#
################################################################################

# Aesthetics
theme_clean <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    panel.border = element_rect(colour = "grey", fill = NA, linewidth = 0.5),
    legend.position = "top",
    legend.title = element_blank(),
    plot.margin = margin(6, 6, 6, 6)
  )

square_panel <- theme(
  aspect.ratio = 1
)

# Labels for axes
x_lab_short_logSI   <- expression("ln(Sampling intensity)")
x_lab_short_logN    <- expression("ln(Network size)")

# Function for scatter plot with LM annotation
scatter_lm_onefit <- function(df, xvar, yvar, ylab, x_label_expr = x_lab_short_logSI,
                              text_size = 4, x_pos = "right", y_pos = 0.9) {
  ggplot(df, aes_string(x = xvar, y = yvar)) +
    geom_point(aes(fill = TYPE, color = TYPE), shape = 21, size = 3, na.rm = TRUE) +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black") +
    stat_poly_eq(
      aes_string(
        x = xvar, y = yvar,
        label = "paste(..rr.label.., ..p.value.label.., sep = ' ~~~ ')"
      ),
      formula = y ~ x, parse = TRUE,
      label.x.npc = x_pos, label.y.npc = y_pos, size = text_size,
      inherit.aes = FALSE, color = "black"
    ) +
    labs(x = x_label_expr, y = ylab, fill = NULL) +
    scale_fill_manual(
      values = c(
        "Seed-Dispersal" = "#a1dab4",
        "Pollination"    = "#fdbb84"
      ),
      breaks = c("Pollination", "Seed-Dispersal"),
      labels = c(
        "Pollination"    = "Pollination",
        "Seed-Dispersal" = "Seed-dispersal"
      ),
      drop = FALSE
    ) +
    scale_color_manual(
      values = c(
        "Seed-Dispersal" = "#1b9e77",
        "Pollination"    = "#d95f02"
      ),
      breaks = c("Pollination", "Seed-Dispersal"),
      guide = "none",
      drop = FALSE
    ) +
    guides(fill = guide_legend(override.aes = list(color = c("#d95f02", "#1b9e77")))) +
    theme_clean + square_panel
}

################################################
# Sampling intensity vs network size
################################################

fit_df <- results %>%
  filter(is.finite(SI_log), is.finite(N_log))

# Equation
fit_lm_fixed <- lm(SI_log ~ offset(-0.5 * N_log), data = fit_df)
ln_k_hat <- coef(fit_lm_fixed)[1]
k_hat <- exp(ln_k_hat)
fit_df$yhat <- predict(fit_lm_fixed)
R2_reg <- cor(fit_df$SI_log, fit_df$yhat)^2

# Smooth curve in log–log space
curve_df <- data.frame(
  NetworkSize = seq(min(results$NetworkSize, na.rm = TRUE),
                    max(results$NetworkSize, na.rm = TRUE),
                    length.out = 400)
) %>%
  mutate(SI = k_hat / sqrt(NetworkSize),
         N_log = log(NetworkSize),
         SI_log = log(SI))

# Annotation label
k_hat <- as.numeric(exp(coef(fit_lm_fixed)[1]))
k_fmt  <- format(round(k_hat, 2), nsmall = 2)
R2_fmt <- format(round(R2_reg, 2), nsmall = 2)

# Label to go on plot
label_expr <- bquote(
  atop(
    "ln(Sampling intensity)" == ln(k) - 0.5 %.% "ln(Network size)",
    k == .(k_fmt) ~~ italic(R)^2 == .(R2_fmt)
  )
)

label_all <- as.character(as.expression(label_expr))

p_top <- ggplot(results, aes(x = N_log, y = SI_log)) +
  geom_point(aes(fill = TYPE, color = TYPE), shape = 21, size = 3, na.rm = TRUE) +
  geom_line(data = curve_df, aes(x = N_log, y = SI_log), color = "black", linewidth = 1) +
  annotate("text",
           x = Inf, y = Inf,
           label = label_all, parse = TRUE,
           hjust = 1.05, vjust = 1.8, size = 4) +
  labs(x = x_lab_short_logN, y = x_lab_short_logSI, fill = NULL) +
  scale_fill_manual(
    values = c(
      "Seed-Dispersal" = "#a1dab4",
      "Pollination"    = "#fdbb84"
    ),
    breaks = c("Pollination", "Seed-Dispersal"),
    labels = c(
      "Pollination"    = "Pollination",
      "Seed-Dispersal" = "Seed-dispersal"
    ),
    drop = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Seed-Dispersal" = "#1b9e77",
      "Pollination"    = "#d95f02"
    ),
    breaks = c("Pollination", "Seed-Dispersal"),
    guide = "none",
    drop = FALSE
  ) +
  guides(fill = guide_legend(override.aes = list(color = c("#d95f02", "#1b9e77")))) +
  theme_clean + 
  square_panel +
  coord_cartesian(clip = "off")

#############################################################
# H2', wNODF, DIRTLPAwb+ vs sampling intensity
#############################################################

p_h2_log <- scatter_lm_onefit(results, "SI_log", "H2",
                                 ylab = bquote("Specialization (" * H[2] * "')"))
p_wNODF_log <- scatter_lm_onefit(results, "SI_log", "weighted_NODF",
                                 ylab = "Weighted nestedness (wNODF)")
p_mod_log <- scatter_lm_onefit(results, "SI_log", "weighted_modularity",
                                 ylab = "Weighted modularity (DIRTLPAwb+)")

##########################################################
# All panels together
##########################################################

final_plot_log <- (p_top | p_h2_log) /
  (p_wNODF_log | p_mod_log) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.background = element_rect(colour = "grey", fill = NA, linewidth = 0.5)
  )

print(final_plot_log)

ggsave(
  "sampling_intensity_topological_indices.pdf",
  plot = last_plot(),
  width = 9.2,
  height = 9.2,
  units = "in"
)
