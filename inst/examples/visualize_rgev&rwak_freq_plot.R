####################################################
### 패키지 로드
####################################################
library(dplyr)    # 데이터 전처리
library(tidyr)    # wide -> long 변환
library(ggplot2)  # 시각화
library(ggh4x)    # facet_grid2
library(grid)     # unit

####################################################
### 0. 사용자 설정
####################################################
# RDS 파일 경로
file_rgev <- "E:/R simulation/simu_260227/simu_result/result.rsel.rgev0.8.rds"
file_rwak <- "E:/R simulation/simu_260227/simu_result/result.rsel.rwak0.8.rds"

# 저장 폴더
plot_dir <- "E:/R simulation/simu_260227/simu_result/plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# true r, max r
true_r <- 5
max_r  <- 10

####################################################
### 1. 공통 스타일 설정
####################################################
lg_legend_title <- 26
lg_legend_text  <- 24
lg_axis_title   <- 26
lg_axis_text    <- 24
lg_axis_text_y  <- 20
lg_strip_text   <- 26
lg_spacing_y    <- 1.8
lg_spacing_x    <- 0.2

method_colors <- c(
  ED       = "#2E86AB",
  CCDF     = "#A23B72",
  Spacings = "#F18F01"
)

custom_theme <- theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(
      fill  = "gray90",
      color = "black"
    )
  )

open_plot_window <- function(width = 16, height = 10, title = "") {
  os <- .Platform$OS.type
  
  if (os == "windows") {
    windows(width = width, height = height)
  } else if (Sys.info()[["sysname"]] == "Darwin") {
    quartz(width = width, height = height)
  } else {
    x11(width = width, height = height)
  }
}

####################################################
### 2. RDS 읽기 + 전처리 함수
####################################################
prepare_freq_data <- function(file_path,
                              dataset_type = c("rgev", "rwak")) {
  
  dataset_type <- match.arg(dataset_type)
  
  # 데이터 읽기
  df <- readRDS(file_path)
  df <- as.data.frame(df)
  
  # 열 이름이 없거나 비어 있으면 강제 부여
  if (is.null(colnames(df)) || any(colnames(df) == "")) {
    colnames(df) <- c("n", "ipar", "method", paste0("ry", 1:10))
  }
  
  # 필요한 r 컬럼
  r_cols <- paste0("ry", 1:10)
  
  # long 변환
  freq_combined <- df %>%
    filter(!is.na(n)) %>%
    pivot_longer(
      cols      = all_of(r_cols),
      names_to  = "r_value",
      values_to = "Frequency"
    ) %>%
    mutate(
      r_value = as.integer(gsub("ry", "", r_value)),
      n       = paste0("n=", n)
    )
  
  # method 라벨
  freq_combined <- freq_combined %>%
    mutate(
      Method_plot = case_when(
        method == 1 ~ "ED",
        method == 2 ~ "CCDF",
        method == 3 ~ "Spacings",
        TRUE        ~ NA_character_
      ),
      Method_plot = factor(
        Method_plot,
        levels = c("ED", "CCDF", "Spacings")
      )
    )
  
  # 시나리오 라벨
  if (dataset_type == "rgev") {
    freq_combined <- freq_combined %>%
      mutate(
        Scenario_par = paste0("par", ipar),
        Scenario_par = factor(
          Scenario_par,
          levels = paste0("par", sort(unique(ipar)))
        )
      )
  }
  
  if (dataset_type == "rwak") {
    freq_combined <- freq_combined %>%
      mutate(
        Scenario_par = paste0("WA-", ipar),
        Scenario_par = factor(
          Scenario_par,
          levels = paste0("WA-", sort(unique(ipar)))
        )
      )
  }
  
  return(freq_combined)
}

####################################################
### 3. 플롯 함수
### - 비율(%) 기준으로 그림
### - true r = 5 위치에 점선 수직선 추가
####################################################
plot_all_scenarios_combined <- function(data,
                                        truer = 5,
                                        max_r = 10,
                                        plot_title = NULL) {
  
  plot_data <- data %>%
    dplyr::filter(!(r_value %in% c(0, 1, 2))) %>%
    dplyr::filter(r_value <= max_r) %>%
    dplyr::filter(!(n %in% c("n=30", "n=50", "n=80") & r_value == 10)) %>%
    dplyr::filter(!(n == "n=30" & Method_plot == "ED")) %>%
    dplyr::filter(!is.na(Method_plot))
  
  n_levels      <- c("n=30", "n=50", "n=80", "n=110")
  scen_levels   <- levels(data$Scenario_par)
  method_levels <- c("ED", "CCDF", "Spacings")
  r_levels_num  <- sort(unique(plot_data$r_value))
  r_levels_chr  <- as.character(r_levels_num)
  
  plot_data$n            <- factor(plot_data$n, levels = n_levels)
  plot_data$Method_plot  <- factor(plot_data$Method_plot, levels = method_levels)
  plot_data$Scenario_par <- factor(plot_data$Scenario_par, levels = scen_levels)
  plot_data$r_value      <- factor(plot_data$r_value, levels = r_levels_chr)
  
  # 중복행 대비 집계
  plot_data <- plot_data %>%
    dplyr::group_by(Scenario_par, n, Method_plot, r_value) %>%
    dplyr::summarise(
      Frequency = sum(Frequency, na.rm = TRUE),
      .groups   = "drop"
    )
  
  # 각 (Scenario_par, n, Method_plot) 내에서 비율(%) 계산
  plot_data <- plot_data %>%
    dplyr::group_by(Scenario_par, n, Method_plot) %>%
    dplyr::mutate(
      Prop = 100 * Frequency / sum(Frequency, na.rm = TRUE)
    ) %>%
    dplyr::ungroup()
  
  # discrete x축에서 true r 위치
  true_r_pos <- match(as.character(truer), r_levels_chr)
  
  p <- ggplot(plot_data, aes(x = r_value, y = Prop, fill = Method_plot)) +
    geom_col(
      position = position_dodge2(padding = 0.001, preserve = "single"),
      width    = 0.95,
      alpha    = 0.90
    ) +
    geom_vline(
      xintercept = true_r_pos,
      linetype   = "dashed",   # 점선으로 변경
      linewidth  = 1.2,
      color      = "red"
    ) +
    ggh4x::facet_grid2(
      rows   = vars(Scenario_par),
      cols   = vars(n),
      scales = "free_x",
      space  = "free_x",
      axes   = "all"
    ) +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    scale_y_continuous(
      limits   = c(0, 100),
      breaks   = seq(0, 100, 20),
      expand   = c(0, 0),
      position = "left"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    labs(
      x     = "Estimated r",
      y     = "Percentage (%)",
      fill  = "Method",
      title = plot_title
    ) +
    custom_theme +
    theme(
      legend.position = "top",
      legend.title    = element_text(size = lg_legend_title, face = "bold"),
      legend.text     = element_text(size = lg_legend_text),
      axis.title.x    = element_text(size = lg_axis_title, face = "bold"),
      axis.title.y    = element_text(size = lg_axis_title, face = "bold"),
      axis.text.x     = element_text(size = lg_axis_text, face = "bold"),
      axis.text.y     = element_text(size = lg_axis_text_y),
      strip.text.x    = element_text(size = lg_strip_text, face = "bold"),
      strip.text.y    = element_text(size = lg_strip_text, face = "bold", angle = 0),
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
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
    )
  
  return(p)
}

####################################################
### 4. rgev 처리
####################################################
cat("===== rgev plot 생성 시작 =====\n")

freq_rgev <- prepare_freq_data(
  file_path    = file_rgev,
  dataset_type = "rgev"
)

p_rgev <- plot_all_scenarios_combined(
  data       = freq_rgev,
  truer      = true_r,
  max_r      = max_r,
  plot_title = NULL
)

open_plot_window(
  width  = 18,
  height = 12,
  title  = "rgev plot"
)
print(p_rgev)

ggsave(
  filename = "E:/R simulation/simu_260227/simu_result/plots/freq_by_rgev_0.8_plot.png",
  plot     = p_rgev,
  width    = 20,
  height   = 16,
  units    = "in",
  dpi      = 300
)

cat("===== rgev plot 저장 완료 =====\n")

####################################################
### 5. rwak 처리
####################################################
cat("===== rwak plot 생성 시작 =====\n")

freq_rwak <- prepare_freq_data(
  file_path    = file_rwak,
  dataset_type = "rwak"
)

p_rwak <- plot_all_scenarios_combined(
  data       = freq_rwak,
  truer      = true_r,
  max_r      = max_r,
  plot_title = NULL
)

open_plot_window(
  width  = 18,
  height = 12,
  title  = "rwak plot"
)
print(p_rwak)

ggsave(
  filename = "E:/R simulation/simu_260227/simu_result/plots/freq_by_rwak_0.8_plot.png",
  plot     = p_rwak,
  width    = 20,
  height   = 16,
  units    = "in",
  dpi      = 300
)

cat("===== rwak plot 저장 완료 =====\n")
