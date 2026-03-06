# ============================================================
# visualize_simulation_all.R
# paper-style matched to your RDS plotting code
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(grid)

# ============================================================
# 0. 전역 설정
# ============================================================

DATA_DIR <- "E:/R simulation/simu_260227/simu_result"
PLOT_DIR <- "E:/R simulation/simu_260227/simu_result/plots"

TRUE_R <- 4
MAXR   <- 8
SCALE  <- 5

KPA_MAP <- c(
  "neg0p35"    = -0.35,
  "neg0p2"     = -0.20,
  "neg1eneg04" = -1e-4,
  "0p2"        =  0.20,
  "0p35"       =  0.35
)

SCENARIO_LEVELS <- c(
  "PAR1 (k = -0.35)",
  "PAR2 (k = -0.20)",
  "PAR3 (k =~ 0)",
  "PAR4 (k =  0.20)",
  "PAR5 (k =  0.35)"
)

METHOD_COLORS <- c(
  ED       = "#2E86AB",
  CCDF     = "#A23B72",
  Spacings = "#F18F01"
)

if (!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, recursive = TRUE)

# ============================================================
# 1. 논문 스타일 크기 설정
#    (아래 RDS 코드와 동일하게 맞춤)
# ============================================================

lg_legend_title <- 26
lg_legend_text  <- 24
lg_axis_title   <- 26
lg_axis_text    <- 24
lg_axis_text_y  <- 20
lg_strip_text   <- 26
lg_spacing_y    <- 1.8
lg_spacing_x    <- 0.2

# ============================================================
# 2. 파일 진단
# ============================================================

check_files <- function(ts, nsam) {
  cat(sprintf("---- [진단] ts = %s ----\n", ts))
  
  found  <- 0L
  all_na <- 0L
  
  for (n in nsam) {
    for (key in names(KPA_MAP)) {
      
      path <- file.path(
        DATA_DIR,
        sprintf("%s_rgev11_n%d_kpa%s_r_est.csv", ts, n, key)
      )
      
      if (!file.exists(path)) {
        cat(sprintf("  [MISSING]  n=%3d  kpa=%-12s\n", n, key))
        next
      }
      
      found <- found + 1L
      tmp   <- read.csv(path)
      
      na_pct <- round(
        colMeans(is.na(tmp[, c("ed", "ccdf", "spacing")])) * 100
      )
      
      warn <- if (all(na_pct == 100)) {
        "  *** 전부 NA — 재실행 필요 ***"
      } else {
        ""
      }
      
      cat(sprintf(
        "  [OK]  n=%3d  kpa=%-12s  NA(ed=%d%% ccdf=%d%% spacing=%d%%)%s\n",
        n, key, na_pct["ed"], na_pct["ccdf"], na_pct["spacing"], warn
      ))
      
      if (all(na_pct == 100)) all_na <- all_na + 1L
    }
  }
  
  cat(sprintf("  => 발견: %d / 전부NA: %d\n\n", found, all_na))
}

# ============================================================
# 3. 데이터 로딩 + long 변환 + 빈도 집계
# ============================================================

load_freq <- function(ts, nsam) {
  
  raw_list <- lapply(nsam, function(n) {
    lapply(names(KPA_MAP), function(key) {
      
      path <- file.path(
        DATA_DIR,
        sprintf("%s_rgev11_n%d_kpa%s_r_est.csv", ts, n, key)
      )
      
      if (!file.exists(path)) return(NULL)
      
      tmp <- read.csv(path)
      tmp$n   <- n
      tmp$kpa <- KPA_MAP[[key]]
      tmp
    })
  })
  
  df_raw <- bind_rows(
    Filter(Negate(is.null), unlist(raw_list, recursive = FALSE))
  )
  
  if (nrow(df_raw) == 0) {
    message("[경고] 읽힌 데이터 없음 (ts = ", ts, ")")
    return(NULL)
  }
  
  cat(sprintf(
    "  Loaded %d rows | %d files\n",
    nrow(df_raw),
    sum(sapply(unlist(raw_list, recursive = FALSE), Negate(is.null)))
  ))
  
  df_long <- df_raw %>%
    pivot_longer(
      cols      = c(ed, ccdf, spacing),
      names_to  = "method",
      values_to = "r_est"
    ) %>%
    mutate(
      Method_plot = factor(
        case_when(
          method == "ed"      ~ "ED",
          method == "ccdf"    ~ "CCDF",
          method == "spacing" ~ "Spacings"
        ),
        levels = c("ED", "CCDF", "Spacings")
      ),
      
      n_label = factor(
        paste0("n=", n),
        levels = paste0("n=", nsam)
      ),
      
      Scenario = factor(
        case_when(
          kpa == -0.35 ~ "PAR1 (k = -0.35)",
          kpa == -0.20 ~ "PAR2 (k = -0.20)",
          kpa == -1e-4 ~ "PAR3 (k =~ 0)",
          kpa ==  0.20 ~ "PAR4 (k =  0.20)",
          kpa ==  0.35 ~ "PAR5 (k =  0.35)"
        ),
        levels = SCENARIO_LEVELS
      ),
      
      r_est = as.integer(r_est)
    ) %>%
    filter(!is.na(r_est))
  
  if (nrow(df_long) == 0) {
    message("[경고] NA 제거 후 데이터 0행 — CSV 값이 전부 NA (ts = ", ts, ")")
    return(NULL)
  }
  
  df_freq <- df_long %>%
    group_by(Scenario, n_label, Method_plot, r_est) %>%
    summarise(Frequency = n() * SCALE, .groups = "drop") %>%
    complete(
      Scenario,
      n_label,
      Method_plot,
      r_est = 0:MAXR,
      fill = list(Frequency = 0)
    ) %>%
    arrange(Scenario, n_label, Method_plot, r_est)
  
  df_freq
}

# ============================================================
# 4. 공통 테마
#    - 위 RDS 코드 스타일과 맞춤
# ============================================================

custom_theme <- theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(
      fill  = "gray90",
      color = "black"
    )
  )

base_theme <- custom_theme +
  theme(
    legend.position = "top",
    legend.title = element_text(size = lg_legend_title, face = "bold"),
    legend.text  = element_text(size = lg_legend_text),
    
    axis.title.x = element_text(size = lg_axis_title, face = "bold"),
    axis.title.y = element_text(size = lg_axis_title, face = "bold"),
    
    axis.text.x  = element_text(size = lg_axis_text, face = "bold"),
    axis.text.y  = element_text(size = lg_axis_text_y),
    
    strip.text.x = element_text(size = lg_strip_text, face = "bold"),
    strip.text.y = element_text(size = lg_strip_text, face = "bold", angle = 0),
    
    # facet row label이 오른쪽에 오도록 했을 때도 가로 + 굵게 유지
    strip.text.y.right = element_text(
      size  = lg_strip_text,
      face  = "bold",
      angle = 0
    ),
    
    panel.spacing.y = grid::unit(lg_spacing_y, "lines"),
    panel.spacing.x = grid::unit(lg_spacing_x, "lines"),
    
    panel.grid.major.x = element_line(
      color     = "gray70",
      linewidth = 0.7,
      linetype  = "dashed"
    ),
    panel.grid.major.y = element_line(
      color     = "gray85",
      linewidth = 0.5,
      linetype  = "dashed"
    ),
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      fill      = NA,
      color     = "black",
      linewidth = 1
    ),
    
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length = unit(0.16, "cm")
  )

# ============================================================
# 5. 플롯 생성
#    - 오른쪽 y축 Frequency 제거
#    - Scenario strip 글자 가로 + bold
# ============================================================

build_plot <- function(df, subtitle_txt = NULL) {
  
  if (is.null(df) || nrow(df) == 0) {
    stop("build_plot: 데이터가 비어 있습니다.")
  }
  
  ymax <- max(df$Frequency, na.rm = TRUE)
  ymax <- ceiling(ymax / 50) * 50
  
  ggplot(
    df,
    aes(x = r_est, y = Frequency, fill = Method_plot)
  ) +
    geom_col(
      position = position_dodge2(
        padding  = 0.001,
        preserve = "single"
      ),
      width = 0.95,
      alpha = 0.90
    ) +
    geom_vline(
      xintercept = TRUE_R,
      color      = "red",
      linewidth  = 1.0,
      linetype   = "dashed"
    ) +
    ggh4x::facet_grid2(
      rows   = vars(Scenario),
      cols   = vars(n_label),
      scales = "fixed",
      axes   = "all",
      switch = NULL
    ) +
    scale_fill_manual(values = METHOD_COLORS, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, ymax),
      breaks = pretty(c(0, ymax), n = 5),
      expand = c(0, 0),
      position = "left"
    ) +
    scale_x_continuous(
      breaks = 0:MAXR,
      limits = c(0, MAXR),
      expand = c(0, 0)
    ) +
    labs(
      x        = "Estimated r",
      y        = "Frequency",
      fill     = "Method",
      subtitle = subtitle_txt
    ) +
    base_theme
}

# ============================================================
# 6. 플롯 저장
#    - 저장 크기도 위 코드 느낌에 맞춤
# ============================================================

save_plot <- function(p, filename, width = 20, height = 16) {
  out <- file.path(PLOT_DIR, filename)
  ggsave(
    filename = out,
    plot     = p,
    width    = width,
    height   = height,
    units    = "in",
    dpi      = 250
  )
  cat(sprintf("  Saved: %s\n\n", filename))
}

# ============================================================
# 7. fig 설정
# ============================================================

fig_configs <- list(
  fig1 = list(
    ts       = "20260301_151000",
    nsam     = c(30, 50, 80),
    subtitle = NULL,
    file     = "NEW_fig1_freq_plot.png",
    width    = 20,
    height   = 16
  ),
  fig2 = list(
    ts       = "20260304_164031",
    nsam     = c(30, 50, 80),
    subtitle = NULL,
    file     = "NEW_fig2_freq_plot.png",
    width    = 20,
    height   = 16
  ),
  fig3 = list(
    ts       = "20260304_225145",
    nsam     = c(30, 50, 80, 110),
    subtitle = NULL,
    file     = "NEW_fig3_freq_plot.png",
    width    = 20,
    height   = 16
  )
)

# ============================================================
# 8. 메인 실행
# ============================================================

for (fig_name in names(fig_configs)) {
  
  cfg <- fig_configs[[fig_name]]
  
  cat(sprintf("========== %s ==========\n", fig_name))
  
  check_files(cfg$ts, cfg$nsam)
  
  df <- load_freq(cfg$ts, cfg$nsam)
  
  if (!is.null(df)) {
    p <- build_plot(df, cfg$subtitle)
    save_plot(p, cfg$file, width = cfg$width, height = cfg$height)
  } else {
    cat(sprintf("  [skip] %s: 유효한 데이터 없음\n\n", fig_name))
  }
}

cat("============================================================\n")
cat("All done. Figures saved to:", PLOT_DIR, "\n")
cat("============================================================\n")