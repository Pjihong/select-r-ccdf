# ============================================================
#  selectrgev 패키지 설치 및 시뮬레이션 실행 스크립트
#  GitHub: Pjihong/select-r-ccdf
#  저장 위치: E:/R simulation/simu_260227/simu_result/
# ============================================================

# ------------------------------------------------------------
# 0. 필요 패키지 설치 및 로드
# ------------------------------------------------------------

auto_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

auto_install("remotes")
auto_install("goftest")
auto_install("Rsolnp")
auto_install("lmomco")
auto_install("eva")
auto_install("extRemes")

# ------------------------------------------------------------
# 1. GitHub에서 selectrgev 패키지 설치
# ------------------------------------------------------------

remotes::install_github("Pjihong/select-r-ccdf",
                        force   = TRUE,
                        upgrade = "never")

# ------------------------------------------------------------
# 2. 패키지 로드
# ------------------------------------------------------------

library(selectrgev)
library(goftest)
library(Rsolnp)
library(lmomco)
library(eva)
library(extRemes)

# ------------------------------------------------------------
# 3. 저장 경로 설정 및 폴더 생성
# ------------------------------------------------------------

save_dir <- "E:/R simulation/simu_260227/simu_result"

if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
  cat("폴더 생성:", save_dir, "\n")
}

# 실행 시작 시각 기록
run_start_time <- Sys.time()
run_timestamp  <- format(run_start_time, "%Y%m%d_%H%M%S")  # 예: 20260227_143022

# ------------------------------------------------------------
# 4. 시뮬레이션 파라미터 설정
# ------------------------------------------------------------

truer    <- rcontam <- 4
maxr     <- 8
sigL     <- sigL_ed <- 0.05
sigL_cvm <- 0.05
ntest    <- 3
model    <- rLOS <- "rgev11"
maxrun   <- 200

kpa  <- c(-0.35, -0.2, -1e-4, 0.2, 0.35)
nsam <- c(30, 50, 80)

# 파라미터 요약 문자열 (파일명에 사용)
param_tag <- paste0(
  "model-", model,
  "_truer-", truer,
  "_maxr-", maxr,
  "_run-", maxrun,
  "_n-", paste(nsam, collapse="+"),
  "_kpa-", paste(kpa, collapse="+")
)

# ------------------------------------------------------------
# 5. 실험 설정 메타데이터 저장 (txt)
# ------------------------------------------------------------

meta_filename <- file.path(save_dir,
                           paste0(run_timestamp, "_", model, "_meta.txt"))

meta_text <- paste0(
  "============================================================\n",
  "  시뮬레이션 실험 설정 기록\n",
  "============================================================\n",
  "실행 일시     : ", format(run_start_time, "%Y년 %m월 %d일 %H:%M:%S"), "\n",
  "GitHub 패키지 : Pjihong/select-r-ccdf\n",
  "모델          : ", model, "\n",
  "------------------------------------------------------------\n",
  "파라미터 설정\n",
  "------------------------------------------------------------\n",
  "truer (실제 r): ", truer,    "\n",
  "rcontam       : ", rcontam,  "\n",
  "maxr          : ", maxr,     "\n",
  "sigL (ED)     : ", sigL_ed,  "\n",
  "sigL (cvm)    : ", sigL_cvm, "\n",
  "maxrun        : ", maxrun,   "\n",
  "nsam          : ", paste(nsam, collapse=", "), "\n",
  "kpa           : ", paste(kpa,  collapse=", "), "\n",
  "검정 방법     : ed, ccdf, spacing\n",
  "오염 모델     : contam_rlosr (mixp=0.5)\n",
  "============================================================\n"
)

writeLines(meta_text, meta_filename)
cat("메타데이터 저장:", meta_filename, "\n")

# ------------------------------------------------------------
# 6. 시뮬레이션 메인 루프
# ------------------------------------------------------------

all_results <- list()

for (insam in 1:length(nsam)) {
  
  n <- nsam[insam]
  cat("\n========================================\n")
  cat("  샘플 크기 n =", n, "\n")
  cat("========================================\n")
  
  for (ipa in 1:length(kpa)) {
    
    para <- c(0, 0.1, 1, 0.02, kpa[ipa])
    cat("\n--- kappa =", kpa[ipa], "---\n")
    
    r_est <- matrix(NA, maxrun, ntest)
    colnames(r_est) <- c("ed", "ccdf", "spacing")
    Bader <- list()
    
    loop_start <- Sys.time()
    
    for (irun in 1:maxrun) {
      
      Bader[[irun]] <- contam_rlosr(n       = n,
                                    maxr    = maxr,
                                    model   = rLOS,
                                    para    = para,
                                    rcontam = rcontam,
                                    mixp    = 0.5)$cdata
      
      xdat <- as.matrix(Bader[[irun]])
      
      r_est[irun, 1] <- tryCatch(
        rsel.rgev11(xdat, model='rgev11', method="ed",
                    seq.cut=TRUE, qqplot=FALSE, sigL=sigL_cvm)$r.sel,
        error = function(e) NA)
      
      r_est[irun, 2] <- tryCatch(
        rsel.rgev11(xdat, model='rgev11', method="ccdf",
                    seq.cut=TRUE, qqplot=FALSE, sigL=sigL_cvm)$r.sel,
        error = function(e) NA)
      
      r_est[irun, 3] <- tryCatch(
        rsel.rgev11(xdat, model='rgev11', method="spacing",
                    seq.cut=TRUE, qqplot=FALSE, sigL=sigL_cvm)$r.sel,
        error = function(e) NA)
      
      cat("irun, truer, r(ed/ccdf/spacing)=",
          irun, truer, r_est[irun, ], "\n")
      
    } # end for irun
    
    loop_end  <- Sys.time()
    loop_time <- round(as.numeric(difftime(loop_end, loop_start, units="mins")), 2)
    
    # ---- 성능 지표: Bias, SE, RMSE ----
    bias <- rmse <- se <- rep(NA, ntest)
    for (k in 1:ntest) {
      bias[k] <- mean(r_est[, k] - truer, na.rm=TRUE)
      rmse[k] <- sqrt(bias[k]^2 + var(r_est[, k], na.rm=TRUE))
    }
    se <- sqrt(rmse^2 - bias^2)
    
    cat("\n[결과] n =", n, "/ kappa =", kpa[ipa], "\n")
    cat("       방법:        ed      ccdf  spacing\n")
    cat("       bias  =", round(bias, 4), "\n")
    cat("       se    =", round(se,   4), "\n")
    cat("       rmse  =", round(rmse, 4), "\n")
    cat("       소요시간:", loop_time, "분\n\n")
    
    # 결과 저장
    key <- paste0("n", n, "_kpa", kpa[ipa])
    all_results[[key]] <- list(
      n         = n,
      kpa       = kpa[ipa],
      para      = para,
      r_est     = r_est,
      bias      = bias,
      se        = se,
      rmse      = rmse,
      loop_time_min = loop_time
    )
    
    # ---- 조건별 r_est CSV 즉시 저장 ----
    # 파일명: 20260227_143022_rgev11_n30_kpa-0.35_r_est.csv
    kpa_str <- gsub("-", "neg", gsub("\\.", "p", as.character(kpa[ipa])))
    csv_filename <- file.path(save_dir,
                              paste0(run_timestamp, "_", model,
                                     "_n", n, "_kpa", kpa_str, "_r_est.csv"))
    
    r_est_df <- as.data.frame(r_est)
    r_est_df$irun <- 1:maxrun
    r_est_df <- r_est_df[, c("irun", "ed", "ccdf", "spacing")]
    write.csv(r_est_df, file=csv_filename, row.names=FALSE)
    cat("  CSV 저장:", basename(csv_filename), "\n")
    
  } # end for ipa
  
} # end for insam

# ------------------------------------------------------------
# 7. 전체 결과 RDS 저장 (재분석용)
# ------------------------------------------------------------

rds_filename <- file.path(save_dir,
                          paste0(run_timestamp, "_", model,
                                 "_run", maxrun, "_allresults.rds"))

saveRDS(all_results, file=rds_filename)
cat("\nRDS 저장:", rds_filename, "\n")

# ------------------------------------------------------------
# 8. 요약 테이블 CSV 저장
# ------------------------------------------------------------

summary_rows <- list()

for (key in names(all_results)) {
  res <- all_results[[key]]
  summary_rows[[key]] <- data.frame(
    조건     = key,
    n        = res$n,
    kpa      = res$kpa,
    bias_ed      = round(res$bias[1], 4),
    bias_ccdf    = round(res$bias[2], 4),
    bias_spacing = round(res$bias[3], 4),
    se_ed        = round(res$se[1],   4),
    se_ccdf      = round(res$se[2],   4),
    se_spacing   = round(res$se[3],   4),
    rmse_ed      = round(res$rmse[1], 4),
    rmse_ccdf    = round(res$rmse[2], 4),
    rmse_spacing = round(res$rmse[3], 4),
    소요시간_분  = res$loop_time_min
  )
}

summary_df <- do.call(rbind, summary_rows)
rownames(summary_df) <- NULL

summary_filename <- file.path(save_dir,
                              paste0(run_timestamp, "_", model,
                                     "_run", maxrun, "_summary.csv"))

write.csv(summary_df, file=summary_filename, row.names=FALSE)
cat("요약 CSV 저장:", summary_filename, "\n")

# ------------------------------------------------------------
# 9. 전체 실행 완료 메시지
# ------------------------------------------------------------

run_end_time   <- Sys.time()
total_time_min <- round(as.numeric(difftime(run_end_time, run_start_time, units="mins")), 2)

cat("\n============================================================\n")
cat("  시뮬레이션 완료\n")
cat("  시작:", format(run_start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("  종료:", format(run_end_time,   "%Y-%m-%d %H:%M:%S"), "\n")
cat("  총 소요시간:", total_time_min, "분\n")
cat("  저장 위치:", save_dir, "\n")
cat("============================================================\n")

# 저장된 파일 목록 출력
cat("\n저장된 파일 목록:\n")
saved_files <- list.files(save_dir,
                          pattern = run_timestamp,
                          full.names = FALSE)
for (f in saved_files) cat("  -", f, "\n")