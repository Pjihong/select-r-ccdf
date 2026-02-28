# test_install.R
# 패키지 설치 후 간단 테스트 스크립트
# 실행: Rscript test_install.R

cat("=== selectrgev 설치 및 기본 테스트 ===\n")

# 1. GitHub 에서 설치 (처음 한 번만)
# install.packages("remotes")
# remotes::install_github("Pjihong/select-r-ccdf")

# 로컬 개발 중일 때
# devtools::load_all(".")

library(selectrgev)

cat("패키지 로드 성공\n")

# ---- 테스트 1: RGEV (stationary) ----
cat("\n--- 테스트 1: rsel.rgev (stationary RGEV) ---\n")
library(eva)
set.seed(42)
xdat1 <- gevrSim(n = 50, r = 8, gumbel = TRUE)   # Gumbel (kappa=0) data

res1 <- rsel.rgev(xdat1, method = "ed", sigL = 0.05)
cat("Selected r =", res1$r.sel, "\n")
cat("MLE (mu, sigma, kappa) =", round(res1$mle, 4), "\n")
cat("p-values =", round(res1$pval, 4), "\n")

# ---- 테스트 2: RGEV11 (nonstationary) ----
cat("\n--- 테스트 2: rsel.rgev11 (nonstationary RGEV11) ---\n")
set.seed(123)
para <- c(0, 0.1, 1, 0.02, -0.1)   # (mu0, mu1, sigma0, sigma1, kappa)
sim  <- ran.rgev11(para, nsample = 40, rlarg = 8)
xdat2 <- sim$xdat

res2 <- rsel.rgev11(xdat2, model = "rgev11",
                    method = "ed", sigL = 0.05)
cat("Selected r =", res2$r.sel, "\n")
cat("MLE =", round(res2$mle, 4), "\n")
cat("p-values =", round(res2$pval, 4), "\n")

# ---- 테스트 3: ccdf method ----
cat("\n--- 테스트 3: rsel.rgev11 (ccdf method) ---\n")
res3 <- rsel.rgev11(xdat2, model = "rgev11",
                    method = "ccdf", sigL = 0.05)
cat("Selected r (ccdf) =", res3$r.sel, "\n")

cat("\n=== 모든 테스트 통과 ===\n")

# ---- 테스트 4: 실제 데이터 (산청) ----
cat("\n--- 테스트 4: sancheong 실제 데이터 ---\n")
data(sancheong)
cat("데이터 크기:", nrow(sancheong), "행 x", ncol(sancheong), "열\n")
cat("연도 범위:", min(sancheong$year), "-", max(sancheong$year), "\n")
cat("연최대값(X1) 범위:", min(sancheong$X1), "-", max(sancheong$X1), "mm\n")

xdat_sc <- as.matrix(sancheong[, 2:11])   # X1 ~ X10, r 최대 10

res_sc <- rsel.rgev(xdat_sc, method = "ed", sigL = 0.05)
cat("rsel.rgev  선택 r =", res_sc$r.sel, "\n")

res_sc2 <- rsel.rgev11(xdat_sc, model = "rgev11",
                        method = "ed", sigL = 0.05)
cat("rsel.rgev11 선택 r =", res_sc2$r.sel, "\n")
cat("MLE (mu0,mu1,sig0,sig1,kappa) =", round(res_sc2$mle, 4), "\n")
