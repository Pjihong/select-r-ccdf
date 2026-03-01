# ============================================================
#  visualize_simulation.R
#  Visualization of simulation results for rgev11 model
#
#  Usage:
#    source("inst/examples/visualize_simulation.R")
#
#  Output:
#    figures/fig1_freq_by_method.png
#    figures/fig2_bias_rmse_summary.png
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(grid)

# ------------------------------------------------------------
# 0. User settings — edit here
# ------------------------------------------------------------

data_dir  <- "E:/R simulation/simu_260227/simu_result"  # CSV 저장 경로
ts        <- "20260301_151000"                           # 타임스탬프 prefix
plot_dir  <- "figures"                                   # 그림 저장 경로 (리포지토리 루트 기준)
truer     <- 4                                           # 실제 r 값

# ------------------------------------------------------------
# 1. Read CSV files
# ------------------------------------------------------------

kpa_map <- c(
  "neg0p35"    = -0.35,
  "neg0p2"     = -0.20,
  "neg1eneg04" = -1e-4,
  "0p2"        =  0.20,
  "0p35"       =  0.35
)

nsam   <- c(30, 50, 80)
all_df <- list()

for (n in nsam) {
  for (kpa_str in names(kpa_map)) {
    fname <- file.path(data_dir,
                       paste0(ts, "_rgev11_n", n, "_kpa", kpa_str, "_r_est.csv"))
    if (!file.exists(fname)) {
      warning("File not found: ", fname)
      next
    }
    tmp       <- read.csv(fname)   # columns: irun, ed, ccdf, spacing
    tmp$n     <- n
    tmp$kpa   <- kpa_map[kpa_str]
    all_df[[paste0("n", n, "_", kpa_str)]] <- tmp
  }
}

df_raw <- dplyr::bind_rows(all_df)
cat("Loaded", nrow(df_raw), "rows from", length(all_df), "files\n")

# ------------------------------------------------------------
# 2. Reshape to long format
# ------------------------------------------------------------

df_long <- df_raw %>%
  tidyr::pivot_longer(
    cols      = c(ed, ccdf, spacing),
    names_to  = "method",
    values_to = "r_est"
  ) %>%
  dplyr::mutate(
    Method_plot = dplyr::case_when(
      method == "ed"      ~ "ED",
      method == "ccdf"    ~ "CCDF",
      method == "spacing" ~ "Spacings"
    ),
    Method_plot = factor(Method_plot, levels = c("ED", "CCDF", "Spacings")),
    n_label = factor(paste0("n = ", n), levels = paste0("n = ", nsam)),
    Scenario = dplyr::case_when(
      kpa == -0.35  ~ "WA-1 (k = -0.35)",
      kpa == -0.20  ~ "WA-2 (k = -0.20)",
      kpa == -1e-4  ~ "WA-3 (k =~ 0)",
      kpa ==  0.20  ~ "WA-4 (k =  0.20)",
      kpa ==  0.35  ~ "WA-5 (k =  0.35)"
    ),
    Scenario = factor(Scenario, levels = c(
      "WA-1 (k = -0.35)", "WA-2 (k = -0.20)", "WA-3 (k =~ 0)",
      "WA-4 (k =  0.20)", "WA-5 (k =  0.35)"
    )),
    r_est = as.integer(r_est)
  ) %>%
  dplyr::filter(!is.na(r_est))

# ------------------------------------------------------------
# 3. Frequency table & Bias/RMSE summary
# ------------------------------------------------------------

df_freq <- df_long %>%
  dplyr::group_by(Scenario, n_label, Method_plot, r_est) %>%
  dplyr::summarise(Frequency = dplyr::n(), .groups = "drop")

df_summary <- df_long %>%
  dplyr::group_by(Scenario, n_label, Method_plot) %>%
  dplyr::summarise(
    Bias = mean(r_est - truer, na.rm = TRUE),
    RMSE = sqrt(mean((r_est - truer)^2, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(cols = c(Bias, RMSE),
                      names_to = "Metric", values_to = "Value")

# ------------------------------------------------------------
# 4. Color palette and theme
# ------------------------------------------------------------

method_colors <- c(ED = "#2E86AB", CCDF = "#A23B72", Spacings = "#F18F01")

custom_theme <- ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    panel.grid.major   = ggplot2::element_blank(),
    panel.grid.minor   = ggplot2::element_blank(),
    strip.background   = ggplot2::element_rect(fill = "gray90", color = "black"),
    legend.position    = "top",
    legend.title       = ggplot2::element_text(size = 14, face = "bold"),
    legend.text        = ggplot2::element_text(size = 13),
    axis.title         = ggplot2::element_text(size = 13, face = "bold"),
    axis.text          = ggplot2::element_text(size = 11),
    strip.text         = ggplot2::element_text(size = 12, face = "bold"),
    panel.spacing.y    = grid::unit(1.2, "lines"),
    panel.spacing.x    = grid::unit(0.3, "lines"),
    panel.grid.major.x = ggplot2::element_line(color = "gray75", linewidth = 0.5,
                                                linetype = "dashed"),
    panel.grid.major.y = ggplot2::element_line(color = "gray85", linewidth = 0.4,
                                                linetype = "dashed"),
    panel.border       = ggplot2::element_rect(fill = NA, color = "black",
                                                linewidth = 0.8)
  )

# ------------------------------------------------------------
# 5. Figure 1 — Frequency bar chart
# ------------------------------------------------------------

p_freq <- ggplot2::ggplot(df_freq,
           ggplot2::aes(x = factor(r_est), y = Frequency, fill = Method_plot)) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge2(padding = 0.05, preserve = "single"),
    width = 0.9, alpha = 0.9
  ) +
  ggh4x::facet_grid2(
    rows   = ggplot2::vars(Scenario),
    cols   = ggplot2::vars(n_label),
    scales = "free_x", space = "free_x", axes = "all"
  ) +
  ggplot2::scale_fill_manual(values = method_colors, drop = FALSE) +
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::scale_x_discrete(expand = c(0, 0)) +
  ggplot2::labs(
    x     = "Estimated r",
    y     = "Frequency (out of 200 simulations)",
    fill  = "Method",
    title = "Distribution of Estimated r — rgev11 model (true r = 4)"
  ) +
  custom_theme

# ------------------------------------------------------------
# 6. Figure 2 — Bias / RMSE line plot
# ------------------------------------------------------------

p_summary <- ggplot2::ggplot(df_summary,
             ggplot2::aes(x = n_label, y = Value,
                          color = Method_plot, group = Method_plot)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 3) +
  ggplot2::facet_grid(Metric ~ Scenario, scales = "free_y") +
  ggplot2::scale_color_manual(values = method_colors) +
  ggplot2::labs(
    x     = "Sample Size",
    y     = "Value",
    color = "Method",
    title = "Bias and RMSE by Scenario and Sample Size"
  ) +
  custom_theme

# ------------------------------------------------------------
# 7. Save figures
# ------------------------------------------------------------

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

ggplot2::ggsave(file.path(plot_dir, "fig1_freq_by_method.png"),
                plot = p_freq,    width = 18, height = 14, dpi = 200)
ggplot2::ggsave(file.path(plot_dir, "fig2_bias_rmse_summary.png"),
                plot = p_summary, width = 16, height = 7,  dpi = 200)

cat("Figures saved to:", plot_dir, "\n")
