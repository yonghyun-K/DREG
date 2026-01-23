# ============================================================
# Supplementary figures (extended): WLS, K-sensitivity, and
# nonparametric predictors (spline / random forest)
#
# This script extends `supp_simulation_figures_fixed.R` so you can
# directly generate *additional* Supplementary Material (SM) figures:
#   (1) WLS versions of the main grid figures (vary p / vary r)
#   (2) Sensitivity to the number of folds K
#   (3) A nonlinear DGP + nonparametric working models (spline, RF)
#
# It is written to be PSOCK/SLURM friendly (workers start clean):
#   - all custom functions/objects are explicitly exported in foreach
#   - package functions use explicit namespaces whenever feasible
#
# ---------------------------
# Usage (examples)
# ---------------------------
# 0) Baseline (same as before; OLS; produces fig/plot_p.png, plot_r.png, ...)
#    Rscript supp_simulation_figures_SM.R --task=grid --simnum=500
#
# 1) WLS grid figures (adds *_wls.png variants)
#    Rscript supp_simulation_figures_SM.R --task=grid --fit=wls --simnum=500
#
# 2) K-sensitivity (fixed p=90, r=-0.75 by default)
#    Rscript supp_simulation_figures_SM.R --task=K --K_grid=2,5,10,20 --simnum=500
#    Rscript supp_simulation_figures_SM.R --task=K --fit=wls --K_grid=2,5,10,20 --simnum=500
#
# 3) Nonlinear DGP + spline/RF working models (keep SIMNUM modest)
#    Rscript supp_simulation_figures_SM.R --task=nonlinear --simnum=200 --K=5 --np_dim=5 --spline_df=4 --rf_trees=200
#
# ---------------------------
# Outputs (default outdir=fig)
# ---------------------------
# Existing (grid / OLS):
#   fig/plot_p.png,  fig/plot_r.png
#   fig/plot_p2.png, fig/plot_r2.png
# Grid / WLS:
#   fig/plot_p_wls.png,  fig/plot_r_wls.png
#   fig/plot_p2_wls.png, fig/plot_r2_wls.png
# K-sensitivity:
#   fig/supp_K_stratified.png,  fig/supp_K_rejective.png
#   fig/supp_K_stratified_wls.png, fig/supp_K_rejective_wls.png
# Nonlinear + nonparam:
#   fig/supp_nonlinear_stratified.png, fig/supp_nonlinear_rejective.png
#   (and *_wls variants if --fit=wls)
# ============================================================

## ---------------------------
## Small utilities
## ---------------------------
`%||%` <- function(x, y) if (!is.null(x) && length(x) && !is.na(x[1])) x else y

parse_kv_args <- function(args) {
  # Parses arguments of the form "--key=value".
  out <- list()
  for (a in args) {
    if (!startsWith(a, "--") || !grepl("=", a, fixed = TRUE)) next
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    out[[kv[1]]] <- kv[2]
  }
  out
}

as_logical01 <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  x <- tolower(as.character(x))
  x %in% c("1", "true", "t", "yes", "y")
}

parse_csv_int <- function(x, default) {
  if (is.null(x) || !nzchar(x)) return(default)
  v <- suppressWarnings(as.integer(strsplit(gsub("\\s+", "", x), ",", fixed = TRUE)[[1]]))
  v <- v[is.finite(v) & v > 1]
  if (!length(v)) return(default)
  unique(v)
}

parse_csv_chr <- function(x, default) {
  if (is.null(x) || !nzchar(x)) return(default)
  v <- tolower(strsplit(gsub("\\s+", "", x), ",", fixed = TRUE)[[1]])
  v <- v[nzchar(v)]
  if (!length(v)) return(default)
  unique(v)
}

## ---------------------------
## Parameters (defaults)
## ---------------------------
DEFAULTS <- list(
  task     = "grid",     # grid | K | nonlinear | all
  fit      = "ols",      # ols | wls | both
  simnum   = 500,
  design   = "all",      # all | stratified | rejective
  vary     = "all",      # all | p | r  (only for task=grid)
  outdir   = "fig",
  base_size = 16,
  K        = 10,          # default K used for task=grid and task=nonlinear
  K_grid   = "2,5,10,20",# only for task=K
  workers  = NA,
  workers_cap = 80,     # cap default workers to avoid R connection exhaustion (R default max connections is typically 128)
  cluster  = "psock",   # psock | fork | auto
  cluster_outfile = "", # where PSOCK worker output is sent ("" -> stdout; use a file path to capture)
  seed_pop = 2,
  seed_sim = 11,
  # Nonparametric / nonlinear options
  np_methods = "spline,rf",  # spline, rf
  np_dim     = 5,             # number of (non-intercept) covariates used for spline/RF
  spline_df  = 4,
  rf_trees   = 200,
  rf_mtry    = NA,
  rf_min_node = 5
)

## ---------------------------
## Packages
## ---------------------------
suppressPackageStartupMessages({
  library(mvtnorm)
  library(glmnet)
  library(ggplot2)
  library(foreach)
  library(doParallel)
  library(doRNG)
  library(sampling)
})

## ---------------------------
## Simulation grid (manuscript values)
## ---------------------------
N <- 1000
n <- 300
s <- 5

seq_p <- c(10, 50, 90, 130, 170, 210, 250)
seq_r <- c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)

p_fixed <- 90
r_fixed <- -0.75

## ---------------------------
## Core generators
## ---------------------------
ar1_cor <- function(p, rho) {
  rho ^ abs(outer(seq_len(p), seq_len(p), "-"))
}

## =========================================================
## Sampling designs
## =========================================================

sample_design_stratified <- function(N, n, len, z_sorted_idx, y, X) {
  n_h <- as.integer(c(15, 20, 30, 35) / 100 * n)
  stopifnot(sum(n_h) == n)

  Index <- integer(n)
  d1 <- numeric(N)

  cumn_h <- cumsum(n_h)
  st <- c(1L, head(cumn_h, -1L) + 1L)

  for (h in 1:4) {
    Idx_z <- (len * (h - 1) + 1):(len * h)
    Idx <- z_sorted_idx[Idx_z]

    d1[Idx] <- len / n_h[h]
    Index[st[h]:cumn_h[h]] <- sample(Idx, size = n_h[h], replace = FALSE)
  }

  list(
    Index = Index,
    d1 = d1,
    d1_s = d1[Index],
    y_s = y[Index],
    X_s = X[Index, , drop = FALSE],
    n_h = n_h,
    st = st,
    cumn_h = cumn_h,
    len = len
  )
}

quad_Omega_strata <- function(e, n_h, len, st, cumn_h) {
  out <- 0
  for (h in 1:4) {
    idx <- st[h]:cumn_h[h]
    eh <- e[idx]
    nh <- n_h[h]

    c_off <- len^2/nh^2 - len*(len-1)/nh/(nh-1)
    d_diag <- (1 - nh/len) * len^2/nh^2

    s1 <- sum(eh)
    s2 <- sum(eh*eh)
    out <- out + c_off * (s1*s1) + (d_diag - c_off) * s2
  }
  out
}

make_pi_from_z <- function(z, n, gamma = 1) {
  m <- gamma * 1 / (1 + exp(-z))
  pi <- n * m / sum(m)

  fixed <- rep(FALSE, length(pi))
  while (any(pi > 1)) {
    idx1 <- which(pi > 1 & !fixed)
    if (!length(idx1)) break
    pi[idx1] <- 1
    fixed[idx1] <- TRUE

    n_rem <- n - sum(pi[fixed])
    if (n_rem <= 0) break

    idx0 <- which(!fixed)
    m0 <- m[idx0]
    pi[idx0] <- n_rem * m0 / sum(m0)
  }

  if (abs(sum(pi) - n) > 1e-8) {
    idx0 <- which(pi < 1)
    pi[idx0] <- pi[idx0] * (n - sum(pi[pi == 1])) / sum(pi[idx0])
  }

  pmax(pmin(pi, 1), 1e-12)
}

sample_design_rejective <- function(N, n, z, y, X, gamma = 1) {
  pik <- make_pi_from_z(z, n, gamma = gamma)
  s01 <- sampling::UPmaxentropy(pik)
  Index <- which(s01 == 1)
  d1 <- 1 / pik

  list(
    Index = Index,
    pi = pik,
    d1 = d1,
    d1_s = d1[Index],
    y_s = y[Index],
    X_s = X[Index, , drop = FALSE]
  )
}

var_poisson_approx <- function(e_s, pi_s) {
  sum((1 - pi_s) / (pi_s^2) * (e_s^2))
}

## =========================================================
## Estimation / inference (linear + lasso)
## =========================================================

run_estimators_oneK <- function(
  sampling,
  N, n, s, y, X, mu,
  Index, d1, d1_s, y_s, X_s,
  pi = NULL,
  n_h = NULL, st = NULL, cumn_h = NULL, len = NULL,
  K = 10,
  wls = FALSE,
  lasso_nfolds = 10
) {
  # ----- Base regression (GREG)
  if (wls) {
    lm_obj <- lm(y_s ~ 0 + X_s, weights = d1_s)
  } else {
    lm_obj <- lm(y_s ~ 0 + X_s)
  }
  beta_hat <- lm_obj$coefficients

  # Oracle regression uses intercept + first s signals
  X_s_oracle <- X[Index, 1:(1 + s), drop = FALSE]
  if (wls) {
    lm_obj_oracle <- lm(y_s ~ 0 + X_s_oracle, weights = d1_s)
  } else {
    lm_obj_oracle <- lm(y_s ~ 0 + X_s_oracle)
  }
  beta_hat_oracle <- lm_obj_oracle$coefficients

  # Core totals
  y_HT   <- sum(y_s * d1_s)
  y_diff <- drop(sum(mu) + sum((y_s - mu[Index]) * d1_s))
  y_GREG <- drop(colSums(X) %*% beta_hat + sum((y_s - drop(X_s %*% beta_hat)) * d1_s))

  y_GREG_oracle <- drop(
    colSums(X[, 1:(1 + s), drop = FALSE]) %*% beta_hat_oracle +
      sum((y_s - drop(X_s_oracle %*% beta_hat_oracle)) * d1_s)
  )

  # Lasso (CV)
  if (wls) {
    cv_model <- glmnet::cv.glmnet(X_s, y_s, weights = d1_s, nfolds = lasso_nfolds)
    best_model <- glmnet::glmnet(X_s, y_s, weights = d1_s, lambda = cv_model$lambda.min)
  } else {
    cv_model <- glmnet::cv.glmnet(X_s, y_s, nfolds = lasso_nfolds)
    best_model <- glmnet::glmnet(X_s, y_s, lambda = cv_model$lambda.min)
  }
  beta_hat_Lasso <- as.vector(coef(best_model))[-1]
  y_Lasso <- drop(colSums(X) %*% beta_hat_Lasso + sum((y_s - drop(X_s %*% beta_hat_Lasso)) * d1_s))

  # ----- Variances (design-specific)
  if (sampling == "stratified") {
    sigma_HT   <- sqrt(quad_Omega_strata(y_s, n_h, len, st, cumn_h))
    sigma_diff <- sqrt(quad_Omega_strata(y_s - mu[Index], n_h, len, st, cumn_h))
    sigma_GREG <- sqrt(quad_Omega_strata(y_s - drop(X_s %*% beta_hat), n_h, len, st, cumn_h))
    sigma_GREG_oracle <- sqrt(quad_Omega_strata(y_s - drop(X_s_oracle %*% beta_hat_oracle), n_h, len, st, cumn_h))
    sigma_Lasso <- sqrt(quad_Omega_strata(y_s - drop(X_s %*% beta_hat_Lasso), n_h, len, st, cumn_h))
  } else {
    sigma_HT   <- sqrt(var_poisson_approx(y_s,               pi[Index]))
    sigma_diff <- sqrt(var_poisson_approx(y_s - mu[Index],   pi[Index]))
    sigma_GREG <- sqrt(var_poisson_approx(y_s - drop(X_s %*% beta_hat),       pi[Index]))
    sigma_GREG_oracle <- sqrt(var_poisson_approx(y_s - drop(X_s_oracle %*% beta_hat_oracle), pi[Index]))
    sigma_Lasso<- sqrt(var_poisson_approx(y_s - drop(X_s %*% beta_hat_Lasso), pi[Index]))
  }

  # ----- Bias diagnostic proxy
  tmpIndex <- numeric(length(y)); tmpIndex[Index] <- d1_s
  bias_proxy <- c(
    HT = 0,
    Diff = 0,
    GREG.oracle = cov(tmpIndex, drop(y - X[, 1:(1 + s), drop = FALSE] %*% beta_hat_oracle)) * N,
    GREG = cov(tmpIndex, drop(y - X %*% beta_hat)) * N,
    GREG.Lasso = cov(tmpIndex, drop(y - X %*% beta_hat_Lasso)) * N
  )

  # FDR/FNR for Lasso selection (X includes explicit intercept as col 1)
  signal_idx <- 2:(1 + s)
  FDR <- sum(beta_hat_Lasso[-signal_idx] != 0) / (length(beta_hat_Lasso) - s)
  FNR <- sum(beta_hat_Lasso[signal_idx] == 0) / s

  # ----- Sample-splitting (SREG) with K folds
  perm <- sample.int(N)
  groups <- cut(seq_along(perm), breaks = K, labels = FALSE)
  SubIndex_list <- split(perm, groups)

  y_sreg_parts <- numeric(K)
  y_sreg_lasso_parts <- numeric(K)

  e_oof <- numeric(length(Index)); names(e_oof) <- as.character(Index)
  e_oof_lasso <- numeric(length(Index)); names(e_oof_lasso) <- as.character(Index)

  bias_parts <- matrix(0, nrow = K, ncol = 2)

  for (k in 1:K) {
    SubIndex <- SubIndex_list[[k]]

    Index_sub1 <- Index[Index %in% SubIndex]
    Index_sub2 <- Index[!(Index %in% SubIndex)]

    y_s1 <- y[Index_sub1]
    X_s1 <- X[Index_sub1, , drop = FALSE]
    d1_s1 <- d1[Index_sub1]

    y_s2 <- y[Index_sub2]
    X_s2 <- X[Index_sub2, , drop = FALSE]
    d1_s2 <- d1[Index_sub2]

    # out-of-fold linear fit
    if (wls) {
      lm2 <- lm(y_s2 ~ 0 + X_s2, weights = d1_s2)
    } else {
      lm2 <- lm(y_s2 ~ 0 + X_s2)
    }
    beta2 <- lm2$coefficients

    y_sreg_parts[k] <- drop(colSums(X[SubIndex, , drop = FALSE]) %*% beta2 +
                              sum((y_s1 - drop(X_s1 %*% beta2)) * d1_s1))

    e_tmp <- drop(y_s1 - X_s1 %*% beta2)
    e_oof[as.character(Index_sub1)] <- e_tmp

    bias_parts[k, 1] <- cov(tmpIndex[SubIndex], drop(y[SubIndex] - X[SubIndex, , drop = FALSE] %*% beta2)) * length(SubIndex)

    # out-of-fold lasso fit
    if (wls) {
      cv2 <- glmnet::cv.glmnet(X_s2, y_s2, weights = d1_s2, nfolds = lasso_nfolds)
      fit2 <- glmnet::glmnet(X_s2, y_s2, weights = d1_s2, lambda = cv2$lambda.min)
    } else {
      cv2 <- glmnet::cv.glmnet(X_s2, y_s2, nfolds = lasso_nfolds)
      fit2 <- glmnet::glmnet(X_s2, y_s2, lambda = cv2$lambda.min)
    }
    beta2_L <- as.vector(coef(fit2))[-1]

    y_sreg_lasso_parts[k] <- drop(colSums(X[SubIndex, , drop = FALSE]) %*% beta2_L +
                                    sum((y_s1 - drop(X_s1 %*% beta2_L)) * d1_s1))

    e_tmp_L <- drop(y_s1 - X_s1 %*% beta2_L)
    e_oof_lasso[as.character(Index_sub1)] <- e_tmp_L

    bias_parts[k, 2] <- cov(tmpIndex[SubIndex], drop(y[SubIndex] - X[SubIndex, , drop = FALSE] %*% beta2_L)) * length(SubIndex)
  }

  y_SREG <- sum(y_sreg_parts)
  y_SREG_Lasso <- sum(y_sreg_lasso_parts)

  e_s_ss  <- e_oof[as.character(Index)]
  e_s_ssl <- e_oof_lasso[as.character(Index)]

  if (sampling == "stratified") {
    sigma_SREG   <- sqrt(quad_Omega_strata(e_s_ss,  n_h, len, st, cumn_h))
    sigma_SREG_L <- sqrt(quad_Omega_strata(e_s_ssl, n_h, len, st, cumn_h))
  } else {
    sigma_SREG   <- sqrt(var_poisson_approx(e_s_ss,  pi[Index]))
    sigma_SREG_L <- sqrt(var_poisson_approx(e_s_ssl, pi[Index]))
  }

  y_res <- c(
    HT = y_HT,
    Diff = y_diff,
    GREG.oracle = y_GREG_oracle,
    GREG = y_GREG,
    SREG = y_SREG,
    GREG.Lasso = y_Lasso,
    SREG.Lasso = y_SREG_Lasso
  )

  sigma_res <- c(
    HT = sigma_HT,
    Diff = sigma_diff,
    GREG.oracle = sigma_GREG_oracle,
    GREG = sigma_GREG,
    SREG = sigma_SREG,
    GREG.Lasso = sigma_Lasso,
    SREG.Lasso = sigma_SREG_L
  )

  bias_proxy <- c(
    bias_proxy,
    SREG = sum(bias_parts[, 1]),
    SREG.Lasso = sum(bias_parts[, 2])
  )

  list(y_res = y_res, sigma_res = sigma_res, bias_proxy = bias_proxy, FDR = FDR, FNR = FNR)
}

## =========================================================
## Parallel helper
## =========================================================
get_workers <- function(simnum, cap = 80L) {
  # Using very large PSOCK clusters can hit R's max connection limit (commonly 128)
  # and/or be inefficient due to overhead. We therefore cap the *default*.
  nc <- parallel::detectCores(logical = TRUE)
  cap <- as.integer(cap)
  max(1L, min(simnum, max(1L, nc - 1L), cap))
}

make_cluster_safe <- function(workers, type = c("PSOCK", "FORK", "AUTO"), outfile = "") {
  type <- toupper(type[1])
  if (type %in% c("AUTO")) {
    type <- if (.Platform$OS.type == "windows") "PSOCK" else "FORK"
  }

  workers <- as.integer(workers)
  workers <- max(1L, workers)

  try_make <- function(tt) {
    if (tt == "FORK") {
      if (.Platform$OS.type == "windows") {
        stop("FORK cluster is not supported on Windows; use PSOCK.")
      }
      parallel::makeForkCluster(workers, outfile = outfile)
    } else {
      parallel::makePSOCKcluster(workers, outfile = outfile)
    }
  }

  # Try once; if we hit the classic 'all connections are in use', fall back or shrink.
  cl <- tryCatch(try_make(type), error = function(e) e)
  if (!inherits(cl, "error")) return(cl)

  msg <- conditionMessage(cl)
  if (grepl("all connections are in use", msg, fixed = TRUE)) {
    message("[cluster] ERROR: 'all connections are in use' while creating a PSOCK cluster.")
    message("[cluster] Fix options: (i) reduce --workers, (ii) run with --cluster=fork (Linux/macOS), (iii) start R with a larger --max-connections.")

    # Try to recover by closing stray connections and shrinking workers.
    try(closeAllConnections(), silent = TRUE)

    # If on Unix and not already fork, try fork.
    if (type != "FORK" && .Platform$OS.type != "windows") {
      message("[cluster] Retrying with FORK cluster (avoids socket connections).")
      cl2 <- tryCatch(try_make("FORK"), error = function(e) e)
      if (!inherits(cl2, "error")) return(cl2)
    }

    # Last resort: shrink workers and retry PSOCK.
    w2 <- max(1L, min(workers, 60L))
    if (w2 < workers) {
      message(sprintf("[cluster] Retrying with fewer workers: %d -> %d", workers, w2))
      workers <- w2
      cl3 <- tryCatch(try_make("PSOCK"), error = function(e) e)
      if (!inherits(cl3, "error")) return(cl3)
      cl <- cl3
    }
  }

  stop(cl)
}

## =========================================================
## Run a grid (vary p or r) for one sampling design
## =========================================================
run_grid <- function(
  SIMNUM,
  sampling,
  vary_over,
  K,
  wls,
  workers_override = NA,
  workers_cap = 80L,
  cluster_type = "PSOCK",
  cluster_outfile = "",
  seed_pop,
  seed_sim
) {
  stopifnot(sampling %in% c("stratified", "rejective"))
  stopifnot(vary_over %in% c("p", "r"))

  # population generation (fixed across all scenarios)
  max_p <- max(seq_p)
  set.seed(seed_pop)
  X0 <- mvtnorm::rmvnorm(
    n = N,
    mean = rep(2, max_p),
    sigma = ar1_cor(max_p, 0.2)
  )
  e <- rnorm(N, 0, 1)

  beta_base <- c(rep(1, s), rep(0, max_p - s))
  mu_base <- X0 %*% beta_base
  y <- drop(mu_base + e)
  t_y <- sum(y)

  len <- round(N / 4)
  grid <- if (vary_over == "p") seq_p else seq_r

  est_names <- c("HT","Diff","GREG.oracle","GREG","SREG","GREG.Lasso","SREG.Lasso")
  BIAS_res <- matrix(NA_real_, nrow = length(est_names), ncol = length(grid), dimnames = list(est_names, as.character(grid)))
  SE_res   <- BIAS_res
  RMSE_res <- BIAS_res
  RB_res   <- BIAS_res
  CR_res   <- BIAS_res

  FDR_mean <- rep(NA_real_, length(grid)); names(FDR_mean) <- as.character(grid)
  FNR_mean <- rep(NA_real_, length(grid)); names(FNR_mean) <- as.character(grid)

  workers <- if (!is.na(workers_override)) {
    max(1L, as.integer(workers_override))
  } else {
    get_workers(SIMNUM, cap = workers_cap)
  }

  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

  cl <- make_cluster_safe(workers, type = cluster_type, outfile = cluster_outfile)
  doParallel::registerDoParallel(cl)
  on.exit({ try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)

  message(sprintf("[%s/%s] SIMNUM=%d, workers=%d, K=%d, fit=%s", sampling, vary_over, SIMNUM, workers, K, if (wls) "WLS" else "OLS"))

  for (gi in seq_along(grid)) {
    val <- grid[gi]
    p_use <- if (vary_over == "p") val else p_fixed
    r_use <- if (vary_over == "r") val else r_fixed

    X <- X0[, 1:p_use, drop = FALSE]
    beta <- c(rep(1, s), rep(0, p_use - s))
    mu <- X %*% beta
    X <- cbind(1, X)

    z <- r_use * e + rnorm(N, 0, sqrt(1 - r_use^2))
    z_sorted_idx <- order(z)

    doRNG::registerDoRNG(seed = seed_sim + gi)

    export_vars <- c(
      # sampling
      "sample_design_stratified", "quad_Omega_strata",
      "sample_design_rejective", "make_pi_from_z", "var_poisson_approx",
      # estimation
      "run_estimators_oneK",
      # objects
      "sampling", "N", "n", "s", "len", "z_sorted_idx", "z", "y", "X", "mu", "t_y", "K", "wls"
    )

    final_res <- foreach::foreach(
      sim = 1:SIMNUM,
      .export = export_vars,
      .packages = c("glmnet"),
      .errorhandling = "pass"
    ) %dopar% {
      if (sampling == "stratified") {
        samp <- sample_design_stratified(N = N, n = n, len = len, z_sorted_idx = z_sorted_idx, y = y, X = X)
        out <- run_estimators_oneK(
          sampling = sampling,
          N = N, n = n, s = s, y = y, X = X, mu = mu,
          Index = samp$Index, d1 = samp$d1, d1_s = samp$d1_s, y_s = samp$y_s, X_s = samp$X_s,
          n_h = samp$n_h, st = samp$st, cumn_h = samp$cumn_h, len = samp$len,
          K = K,
          wls = wls
        )
      } else {
        samp <- sample_design_rejective(N = N, n = n, z = z, y = y, X = X, gamma = 1)
        out <- run_estimators_oneK(
          sampling = sampling,
          N = N, n = n, s = s, y = y, X = X, mu = mu,
          Index = samp$Index, d1 = samp$d1, d1_s = samp$d1_s, y_s = samp$y_s, X_s = samp$X_s,
          pi = samp$pi,
          K = K,
          wls = wls
        )
      }

      list(
        y = out$y_res,
        sigma = out$sigma_res,
        cover = ifelse(abs(out$y_res - t_y) > 1.96 * out$sigma_res, 0, 1),
        FDR = out$FDR,
        FNR = out$FNR
      )
    }

    ok <- vapply(final_res, function(x) is.list(x) && is.numeric(unlist(x$y)), logical(1))
    n_fail <- sum(!ok)
    if (n_fail > 0) message(sprintf("  grid[%d]=%s: failures=%d/%d", gi, as.character(val), n_fail, SIMNUM))

    if (sum(ok) == 0) {
      fail_msgs <- vapply(final_res[!ok], function(x) {
        if (inherits(x, "error") || inherits(x, "condition")) return(conditionMessage(x))
        if (is.character(x) && length(x) == 1) return(x)
        paste0("<", paste(class(x), collapse = ","), ">")
      }, character(1))
      tab <- sort(table(fail_msgs), decreasing = TRUE)
      top <- head(tab, 5)
      message("  Top failure messages:")
      for (j in seq_along(top)) message(sprintf("    - %s  [%d]", names(top)[j], as.integer(top[[j]])))
      stop(sprintf("All %d simulations failed for %s design at %s=%s.", SIMNUM, sampling, vary_over, as.character(val)))
    }

    ys <- do.call(rbind, lapply(final_res[ok], `[[`, "y"))
    sig <- do.call(rbind, lapply(final_res[ok], `[[`, "sigma"))
    covr <- do.call(rbind, lapply(final_res[ok], `[[`, "cover"))

    bias <- colMeans(ys - t_y)
    se   <- apply(ys, 2, function(x) sqrt(var(x) * (length(x) - 1) / length(x)))
    rmse <- apply(ys - t_y, 2, function(x) sqrt(mean(x^2)))

    bias_var <- colMeans(sig^2) - se^2
    rb <- bias_var / se^2
    cr <- colMeans(covr)

    BIAS_res[, gi] <- bias[est_names]
    SE_res[, gi]   <- se[est_names]
    RMSE_res[, gi] <- rmse[est_names]
    RB_res[, gi]   <- rb[est_names]
    CR_res[, gi]   <- cr[est_names]

    FDR_mean[gi] <- mean(vapply(final_res[ok], `[[`, numeric(1), "FDR"))
    FNR_mean[gi] <- mean(vapply(final_res[ok], `[[`, numeric(1), "FNR"))
  }

  list(
    sampling = sampling,
    vary_over = vary_over,
    wls = wls,
    grid = grid,
    BIAS = BIAS_res,
    RMSE = RMSE_res,
    RB = RB_res,
    CR = CR_res,
    FDR = FDR_mean,
    FNR = FNR_mean
  )
}

## =========================================================
## Plotting (grid)
## =========================================================

make_bias_rmse_plot <- function(res, show_legend = TRUE, base_size = 16, legend_position = "right") {
  xlab <- if (res$vary_over == "p") "p" else "r"
  xvals <- res$grid

  keep <- c("GREG", "SREG", "GREG.Lasso", "SREG.Lasso")

  df <- rbind(
    data.frame(
      estimator = rep(rownames(res$RMSE), times = ncol(res$RMSE)),
      x = rep(xvals, each = nrow(res$RMSE)),
      value = as.vector(res$RMSE),
      metric = "RMSE"
    ),
    data.frame(
      estimator = rep(rownames(res$BIAS), times = ncol(res$BIAS)),
      x = rep(xvals, each = nrow(res$BIAS)),
      value = as.vector(res$BIAS),
      metric = "Bias"
    )
  )
  df <- subset(df, estimator %in% keep)
  df$estimator <- factor(df$estimator, levels = keep)
  df$x <- factor(df$x, levels = xvals)

  gg <- ggplot(df, aes(x = x, y = value, group = interaction(estimator, metric), color = estimator, linetype = metric)) +
    geom_line(linewidth = 1) +
    geom_point(data = subset(df, metric == "RMSE"), size = 2) +
    scale_linetype_manual(values = c(RMSE = "solid", Bias = "dashed")) +
    labs(
      x = xlab,
      y = "RMSE (solid) / Bias (dashed)",
      color = "Estimator",
      linetype = NULL
    ) +
    theme_bw(base_size = base_size) +
    theme(
      legend.title = element_text(size = base_size * 0.9),
      legend.text  = element_text(size = base_size * 0.85),
      axis.title   = element_text(size = base_size * 0.95),
      strip.text   = element_text(size = base_size * 0.95)
    )

  if (!show_legend) {
    gg <- gg + theme(legend.position = "none")
  } else {
    gg <- gg + theme(legend.position = legend_position)
  }
  gg
}

save_png <- function(plot_obj, path, width = 7.2, height = 4.2, dpi = 300) {
  ggplot2::ggsave(filename = path, plot = plot_obj, width = width, height = height, dpi = dpi)
}


# Save two ggplot objects side-by-side (left and right) to a single PNG.
# To ensure BOTH plot panels have equal size even when a legend is present,
# we extract the legend (from the right plot) and draw it in a third column.
# This uses only ggplot2 + grid (no gtable/gridExtra helpers required).
extract_legend_grob <- function(p) {
  g <- ggplot2::ggplotGrob(p)
  idx <- which(grepl("guide-box", g$layout$name))
  if (!length(idx)) return(NULL)
  g$grobs[[idx[1]]]
}

save_two_panel_png <- function(p_left, p_right, path, width = 12.0, height = 4.4, dpi = 300) {
  # Pull legend from the right plot (if any), then remove legends from both plots.
  leg <- extract_legend_grob(p_right)

  p_left_noleg  <- p_left  + theme(legend.position = "none")
  p_right_noleg <- p_right + theme(legend.position = "none")

  gl <- ggplot2::ggplotGrob(p_left_noleg)
  gr <- ggplot2::ggplotGrob(p_right_noleg)

  # Align widths/heights so the panel regions match.
  max_heights <- grid::unit.pmax(gl$heights, gr$heights)
  max_widths  <- grid::unit.pmax(gl$widths,  gr$widths)
  gl$heights <- max_heights; gr$heights <- max_heights
  gl$widths  <- max_widths;  gr$widths  <- max_widths

  grDevices::png(filename = path, width = width, height = height, units = "in", res = dpi)
  on.exit(grDevices::dev.off(), add = TRUE)

  grid::grid.newpage()

  if (!is.null(leg)) {
    # Give legend its natural width (+ a touch of padding).
    leg_w <- grid::grobWidth(leg) + grid::unit(0.8, "lines")
    lay <- grid::grid.layout(
      nrow = 1, ncol = 3,
      widths = grid::unit.c(grid::unit(1, "null"), grid::unit(1, "null"), leg_w)
    )
    grid::pushViewport(grid::viewport(layout = lay))

    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid::grid.draw(gl)
    grid::popViewport()

    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    grid::grid.draw(gr)
    grid::popViewport()

    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 3))
    grid::grid.draw(leg)
    grid::popViewport()

    grid::popViewport()
  } else {
    # No legend case: just 2 equal columns.
    lay <- grid::grid.layout(nrow = 1, ncol = 2, widths = grid::unit(c(1, 1), "null"))
    grid::pushViewport(grid::viewport(layout = lay))

    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid::grid.draw(gl)
    grid::popViewport()

    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    grid::grid.draw(gr)
    grid::popViewport()

    grid::popViewport()
  }
}

# Combine two ggplots side-by-side while keeping legend from the right plot.
# Uses only ggplot2 + gtable (dependency of ggplot2), no extra packages.
combine_side_by_side <- function(p_left, p_right) {
  gl <- ggplot2::ggplotGrob(p_left)
  gr <- ggplot2::ggplotGrob(p_right)

  # Align panel heights
  max_heights <- grid::unit.pmax(gl$heights, gr$heights)
  gl$heights <- max_heights
  gr$heights <- max_heights

  # NOTE: In some gtable versions (common on clusters), gtable_cbind() exists
  # in the namespace but is not exported. We therefore look it up from the
  # namespace (works whether or not it is exported). If it doesn't exist (rare),
  # fall back to gridExtra::arrangeGrob if available.
  if (exists("gtable_cbind", envir = asNamespace("gtable"), inherits = FALSE)) {
    f <- get("gtable_cbind", envir = asNamespace("gtable"))
    return(f(gl, gr, size = "first"))
  }
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    return(gridExtra::arrangeGrob(gl, gr, ncol = 2))
  }
  stop("Could not combine plots: gtable_cbind not found in gtable namespace and gridExtra is not installed.")
}

## =========================================================
## Task: K-sensitivity (linear + lasso; fixed p/r)
## =========================================================

run_estimators_multiK <- function(
  sampling,
  N, n, s, y, X, mu,
  Index, d1, d1_s, y_s, X_s,
  pi = NULL,
  n_h = NULL, st = NULL, cumn_h = NULL, len = NULL,
  K_vec,
  wls = FALSE,
  lasso_nfolds = 10
) {
  # Base fits (do once)
  if (wls) {
    lm_obj <- lm(y_s ~ 0 + X_s, weights = d1_s)
  } else {
    lm_obj <- lm(y_s ~ 0 + X_s)
  }
  beta_hat <- lm_obj$coefficients

  if (wls) {
    cv_model <- glmnet::cv.glmnet(X_s, y_s, weights = d1_s, nfolds = lasso_nfolds)
    best_model <- glmnet::glmnet(X_s, y_s, weights = d1_s, lambda = cv_model$lambda.min)
  } else {
    cv_model <- glmnet::cv.glmnet(X_s, y_s, nfolds = lasso_nfolds)
    best_model <- glmnet::glmnet(X_s, y_s, lambda = cv_model$lambda.min)
  }
  beta_hat_Lasso <- as.vector(coef(best_model))[-1]

  y_HT   <- sum(y_s * d1_s)
  y_diff <- drop(sum(mu) + sum((y_s - mu[Index]) * d1_s))
  y_GREG <- drop(colSums(X) %*% beta_hat + sum((y_s - drop(X_s %*% beta_hat)) * d1_s))
  y_Lasso <- drop(colSums(X) %*% beta_hat_Lasso + sum((y_s - drop(X_s %*% beta_hat_Lasso)) * d1_s))

  # Variances for base estimators
  if (sampling == "stratified") {
    sigma_HT   <- sqrt(quad_Omega_strata(y_s, n_h, len, st, cumn_h))
    sigma_diff <- sqrt(quad_Omega_strata(y_s - mu[Index], n_h, len, st, cumn_h))
    sigma_GREG <- sqrt(quad_Omega_strata(y_s - drop(X_s %*% beta_hat), n_h, len, st, cumn_h))
    sigma_Lasso <- sqrt(quad_Omega_strata(y_s - drop(X_s %*% beta_hat_Lasso), n_h, len, st, cumn_h))
  } else {
    sigma_HT   <- sqrt(var_poisson_approx(y_s,               pi[Index]))
    sigma_diff <- sqrt(var_poisson_approx(y_s - mu[Index],   pi[Index]))
    sigma_GREG <- sqrt(var_poisson_approx(y_s - drop(X_s %*% beta_hat),       pi[Index]))
    sigma_Lasso<- sqrt(var_poisson_approx(y_s - drop(X_s %*% beta_hat_Lasso), pi[Index]))
  }

  # Shared permutation for all K (reduces Monte Carlo noise across K)
  perm <- sample.int(N)

  y_SREG_byK <- numeric(length(K_vec)); names(y_SREG_byK) <- paste0("SREG.K", K_vec)
  y_SREG_L_byK <- numeric(length(K_vec)); names(y_SREG_L_byK) <- paste0("SREG.Lasso.K", K_vec)

  sig_SREG_byK <- numeric(length(K_vec)); names(sig_SREG_byK) <- paste0("SREG.K", K_vec)
  sig_SREG_L_byK <- numeric(length(K_vec)); names(sig_SREG_L_byK) <- paste0("SREG.Lasso.K", K_vec)

  tmpIndex <- numeric(length(y)); tmpIndex[Index] <- d1_s

  for (ii in seq_along(K_vec)) {
    K <- K_vec[ii]
    groups <- cut(seq_along(perm), breaks = K, labels = FALSE)
    SubIndex_list <- split(perm, groups)

    y_parts <- numeric(K)
    y_parts_L <- numeric(K)
    e_oof <- numeric(length(Index)); names(e_oof) <- as.character(Index)
    e_oof_L <- numeric(length(Index)); names(e_oof_L) <- as.character(Index)

    for (k in 1:K) {
      SubIndex <- SubIndex_list[[k]]
      Index_sub1 <- Index[Index %in% SubIndex]
      Index_sub2 <- Index[!(Index %in% SubIndex)]

      y_s1 <- y[Index_sub1]
      X_s1 <- X[Index_sub1, , drop = FALSE]
      d1_s1 <- d1[Index_sub1]

      y_s2 <- y[Index_sub2]
      X_s2 <- X[Index_sub2, , drop = FALSE]
      d1_s2 <- d1[Index_sub2]

      # Linear
      if (wls) {
        lm2 <- lm(y_s2 ~ 0 + X_s2, weights = d1_s2)
      } else {
        lm2 <- lm(y_s2 ~ 0 + X_s2)
      }
      beta2 <- lm2$coefficients
      y_parts[k] <- drop(colSums(X[SubIndex, , drop = FALSE]) %*% beta2 +
                           sum((y_s1 - drop(X_s1 %*% beta2)) * d1_s1))
      e_oof[as.character(Index_sub1)] <- drop(y_s1 - X_s1 %*% beta2)

      # Lasso
      if (wls) {
        cv2 <- glmnet::cv.glmnet(X_s2, y_s2, weights = d1_s2, nfolds = lasso_nfolds)
        fit2 <- glmnet::glmnet(X_s2, y_s2, weights = d1_s2, lambda = cv2$lambda.min)
      } else {
        cv2 <- glmnet::cv.glmnet(X_s2, y_s2, nfolds = lasso_nfolds)
        fit2 <- glmnet::glmnet(X_s2, y_s2, lambda = cv2$lambda.min)
      }
      beta2_L <- as.vector(coef(fit2))[-1]
      y_parts_L[k] <- drop(colSums(X[SubIndex, , drop = FALSE]) %*% beta2_L +
                             sum((y_s1 - drop(X_s1 %*% beta2_L)) * d1_s1))
      e_oof_L[as.character(Index_sub1)] <- drop(y_s1 - X_s1 %*% beta2_L)
    }

    y_SREG_byK[ii] <- sum(y_parts)
    y_SREG_L_byK[ii] <- sum(y_parts_L)

    e_s_ss  <- e_oof[as.character(Index)]
    e_s_ssl <- e_oof_L[as.character(Index)]

    if (sampling == "stratified") {
      sig_SREG_byK[ii]   <- sqrt(quad_Omega_strata(e_s_ss,  n_h, len, st, cumn_h))
      sig_SREG_L_byK[ii] <- sqrt(quad_Omega_strata(e_s_ssl, n_h, len, st, cumn_h))
    } else {
      sig_SREG_byK[ii]   <- sqrt(var_poisson_approx(e_s_ss,  pi[Index]))
      sig_SREG_L_byK[ii] <- sqrt(var_poisson_approx(e_s_ssl, pi[Index]))
    }
  }

  y_res <- c(HT = y_HT, Diff = y_diff, GREG = y_GREG, GREG.Lasso = y_Lasso, y_SREG_byK, y_SREG_L_byK)
  sigma_res <- c(HT = sigma_HT, Diff = sigma_diff, GREG = sigma_GREG, GREG.Lasso = sigma_Lasso, sig_SREG_byK, sig_SREG_L_byK)

  list(y_res = y_res, sigma_res = sigma_res)
}

run_K_grid <- function(
  SIMNUM,
  sampling,
  K_vec,
  p_use,
  r_use,
  wls,
  workers_override = NA,
  workers_cap = 80L,
  cluster_type = "PSOCK",
  cluster_outfile = "",
  seed_pop,
  seed_sim
) {
  stopifnot(sampling %in% c("stratified", "rejective"))
  K_vec <- sort(unique(as.integer(K_vec)))

  set.seed(seed_pop)
  X0 <- mvtnorm::rmvnorm(
    n = N,
    mean = rep(2, p_use),
    sigma = ar1_cor(p_use, 0.2)
  )
  e <- rnorm(N, 0, 1)
  beta <- c(rep(1, s), rep(0, p_use - s))
  mu <- X0 %*% beta
  y <- drop(mu + e)
  t_y <- sum(y)
  X <- cbind(1, X0)

  len <- round(N / 4)
  z <- r_use * e + rnorm(N, 0, sqrt(1 - r_use^2))
  z_sorted_idx <- order(z)

  workers <- if (!is.na(workers_override)) {
    max(1L, as.integer(workers_override))
  } else {
    get_workers(SIMNUM, cap = workers_cap)
  }
  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

  cl <- make_cluster_safe(workers, type = cluster_type, outfile = cluster_outfile)
  doParallel::registerDoParallel(cl)
  on.exit({ try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)

  message(sprintf("[%s/K] SIMNUM=%d, workers=%d, fit=%s, K_grid={%s}",
                  sampling, SIMNUM, workers, if (wls) "WLS" else "OLS", paste(K_vec, collapse = ",")))

  doRNG::registerDoRNG(seed = seed_sim + 1000L)

  export_vars <- c(
    "sample_design_stratified", "quad_Omega_strata",
    "sample_design_rejective", "make_pi_from_z", "var_poisson_approx",
    "run_estimators_multiK",
    "sampling", "N", "n", "s", "len", "z_sorted_idx", "z", "y", "X", "mu", "t_y", "wls", "K_vec"
  )

  final_res <- foreach::foreach(
    sim = 1:SIMNUM,
    .export = export_vars,
    .packages = c("glmnet"),
    .errorhandling = "pass"
  ) %dopar% {
    if (sampling == "stratified") {
      samp <- sample_design_stratified(N = N, n = n, len = len, z_sorted_idx = z_sorted_idx, y = y, X = X)
      out <- run_estimators_multiK(
        sampling = sampling,
        N = N, n = n, s = s, y = y, X = X, mu = mu,
        Index = samp$Index, d1 = samp$d1, d1_s = samp$d1_s, y_s = samp$y_s, X_s = samp$X_s,
        n_h = samp$n_h, st = samp$st, cumn_h = samp$cumn_h, len = samp$len,
        K_vec = K_vec,
        wls = wls
      )
    } else {
      samp <- sample_design_rejective(N = N, n = n, z = z, y = y, X = X, gamma = 1)
      out <- run_estimators_multiK(
        sampling = sampling,
        N = N, n = n, s = s, y = y, X = X, mu = mu,
        Index = samp$Index, d1 = samp$d1, d1_s = samp$d1_s, y_s = samp$y_s, X_s = samp$X_s,
        pi = samp$pi,
        K_vec = K_vec,
        wls = wls
      )
    }

    list(
      y = out$y_res,
      sigma = out$sigma_res,
      cover = ifelse(abs(out$y_res - t_y) > 1.96 * out$sigma_res, 0, 1)
    )
  }

  ok <- vapply(final_res, function(x) is.list(x) && is.numeric(unlist(x$y)), logical(1))
  if (sum(ok) == 0) stop("All simulations failed in K-grid run.")

  ys <- do.call(rbind, lapply(final_res[ok], `[[`, "y"))

  bias <- colMeans(ys - t_y)
  se   <- apply(ys, 2, function(x) sqrt(var(x) * (length(x) - 1) / length(x)))
  rmse <- apply(ys - t_y, 2, function(x) sqrt(mean(x^2)))

  list(
    sampling = sampling,
    wls = wls,
    K_vec = K_vec,
    bias = bias,
    rmse = rmse,
    se = se
  )
}


make_K_plot <- function(resK, show_legend = TRUE, base_size = 16, legend_position = "bottom") {
  # Plot RMSE (solid) and Bias (dashed) together against K.
  # Baseline estimators (GREG, GREG.Lasso) are constant in K but are shown as horizontal lines.
  K_vec <- sort(unique(as.integer(resK$K_vec)))

  get_split_vals <- function(prefix, K_vec, vec) {
    nm <- paste0(prefix, K_vec)
    as.numeric(vec[nm])
  }

  df <- rbind(
    data.frame(K = K_vec, estimator = "GREG",       metric = "RMSE", value = rep(unname(resK$rmse["GREG"]),       length(K_vec))),
    data.frame(K = K_vec, estimator = "GREG.Lasso", metric = "RMSE", value = rep(unname(resK$rmse["GREG.Lasso"]), length(K_vec))),
    data.frame(K = K_vec, estimator = "SREG",       metric = "RMSE", value = get_split_vals("SREG.K",       K_vec, resK$rmse)),
    data.frame(K = K_vec, estimator = "SREG.Lasso", metric = "RMSE", value = get_split_vals("SREG.Lasso.K", K_vec, resK$rmse)),
    data.frame(K = K_vec, estimator = "GREG",       metric = "Bias", value = rep(unname(resK$bias["GREG"]),       length(K_vec))),
    data.frame(K = K_vec, estimator = "GREG.Lasso", metric = "Bias", value = rep(unname(resK$bias["GREG.Lasso"]), length(K_vec))),
    data.frame(K = K_vec, estimator = "SREG",       metric = "Bias", value = get_split_vals("SREG.K",       K_vec, resK$bias)),
    data.frame(K = K_vec, estimator = "SREG.Lasso", metric = "Bias", value = get_split_vals("SREG.Lasso.K", K_vec, resK$bias))
  )

  df$type <- ifelse(df$estimator %in% c("GREG", "GREG.Lasso"), "Baseline", "SREG")
  df$metric <- factor(df$metric, levels = c("RMSE", "Bias"))
  df$estimator <- factor(df$estimator, levels = c("GREG", "GREG.Lasso", "SREG", "SREG.Lasso"))

  gg <- ggplot(df, aes(
    x = K, y = value,
    color = estimator,
    linetype = metric,
    group = interaction(estimator, metric)
  )) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey35") +
    geom_line(aes(alpha = type), linewidth = 1) +
    geom_point(data = subset(df, metric == "RMSE" & type == "SREG"), size = 2) +
    scale_alpha_manual(values = c(Baseline = 0.55, SREG = 1), guide = "none") +
    scale_linetype_manual(values = c(RMSE = "solid", Bias = "dashed")) +
    scale_x_continuous(breaks = K_vec) +
    labs(
      x = "K",
      y = "RMSE (solid) / Bias (dashed)",
      color = "Estimator",
      linetype = NULL
    ) +
    theme_bw(base_size = base_size) +
    theme(
      legend.position = if (show_legend) legend_position else "none",
      legend.box = "vertical",
      legend.title = element_text(size = base_size * 0.95),
      legend.text  = element_text(size = base_size * 0.9),
      axis.title   = element_text(size = base_size * 1.0),
      strip.text   = element_text(size = base_size * 1.0)
    )

  gg
}

make_K_plot_combined <- function(resK_list, base_size = 14) {
  designs <- names(resK_list)
  if (is.null(designs) || any(!nzchar(designs))) designs <- paste0("design", seq_along(resK_list))
  
  build_df <- function(resK, design) {
    K_vec <- sort(unique(as.integer(resK$K_vec)))
    
    get_split_vals <- function(prefix, K_vec, vec) {
      nm <- paste0(prefix, K_vec)
      as.numeric(vec[nm])
    }
    
    # RMSE values
    df_rmse <- rbind(
      data.frame(design=design, K=K_vec, estimator="GREG",       metric="RMSE", value=rep(unname(resK$rmse["GREG"]),       length(K_vec))),
      data.frame(design=design, K=K_vec, estimator="GREG.Lasso", metric="RMSE", value=rep(unname(resK$rmse["GREG.Lasso"]), length(K_vec))),
      data.frame(design=design, K=K_vec, estimator="SREG",       metric="RMSE", value=get_split_vals("SREG.K",       K_vec, resK$rmse)),
      data.frame(design=design, K=K_vec, estimator="SREG.Lasso", metric="RMSE", value=get_split_vals("SREG.Lasso.K", K_vec, resK$rmse))
    )
    
    # NEGATIVE bias values (this is the key change)
    df_nbias <- rbind(
      data.frame(design=design, K=K_vec, estimator="GREG",       metric="-Bias", value=rep(-unname(resK$bias["GREG"]),       length(K_vec))),
      data.frame(design=design, K=K_vec, estimator="GREG.Lasso", metric="-Bias", value=rep(-unname(resK$bias["GREG.Lasso"]), length(K_vec))),
      data.frame(design=design, K=K_vec, estimator="SREG",       metric="-Bias", value=-get_split_vals("SREG.K",       K_vec, resK$bias)),
      data.frame(design=design, K=K_vec, estimator="SREG.Lasso", metric="-Bias", value=-get_split_vals("SREG.Lasso.K", K_vec, resK$bias))
    )
    
    df <- rbind(df_rmse, df_nbias)
    df$type <- ifelse(df$estimator %in% c("GREG","GREG.Lasso"), "Baseline", "SREG")
    df
  }
  
  df_all <- do.call(rbind, Map(build_df, resK_list, designs))
  
  df_all$design <- factor(df_all$design, levels = c("stratified","rejective"))
  df_all$metric <- factor(df_all$metric, levels = c("RMSE","-Bias"))
  df_all$estimator <- factor(df_all$estimator, levels = c("GREG","GREG.Lasso","SREG","SREG.Lasso"))
  K_breaks <- sort(unique(as.integer(df_all$K)))
  
  gg <- ggplot2::ggplot(df_all, ggplot2::aes(
    x = K, y = value,
    color = estimator,
    linetype = metric,
    group = interaction(design, estimator, metric)
  )) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "grey35") +
    ggplot2::geom_line(ggplot2::aes(alpha = type), linewidth = 1.0) +
    ggplot2::geom_point(data = subset(df_all, metric == "RMSE" & type == "SREG"), size = 2.2) +
    ggplot2::scale_alpha_manual(values = c(Baseline = 0.55, SREG = 1), guide = "none") +
    ggplot2::scale_linetype_manual(values = c(RMSE = "solid", `-Bias` = "dashed")) +
    ggplot2::scale_x_log10(breaks = K_breaks, labels = K_breaks) +
    ggplot2::facet_wrap(~ design, nrow = 1, scales = "free_y") +
    ggplot2::labs(
      x = "K",
      y = "RMSE (solid) / -Bias (dashed)",
      color = "Estimator",
      linetype = NULL
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      # ylab: smaller + add margin so it doesn't get clipped
      axis.title.y = ggplot2::element_text(size = base_size * 0.95,
                                           margin = ggplot2::margin(r = 8)),
      axis.title.x = ggplot2::element_text(size = base_size * 1.00),
      strip.text   = ggplot2::element_text(size = base_size * 1.00),
      
      # overall plot margins: more left space helps prevent clipping
      plot.margin = ggplot2::margin(t=6, r=6, b=6, l=12),
      
      # legend: make it compact
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = ggplot2::element_text(size = base_size * 0.90),
      legend.text  = ggplot2::element_text(size = base_size * 0.85),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.key.width  = grid::unit(0.80, "cm"),
      legend.spacing.x  = grid::unit(0.15, "cm"),
      legend.margin     = ggplot2::margin(t=2, r=2, b=2, l=2)
    ) +
    ggplot2::guides(
      # force one-row legend (much shorter)
      color = ggplot2::guide_legend(nrow = 1, byrow = TRUE),
      linetype = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
    )
  
  gg
}
## =========================================================
## Task: Nonlinear DGP + nonparametric working models
## =========================================================

mu_nonlinear <- function(X0, s) {
  # X0 is N x p (no intercept). We build a smooth nonlinear mean using first s covariates.
  # Center around the population mean (2) so scales stay reasonable.
  xc <- X0[, 1:s, drop = FALSE] - 2
  # A modest nonlinearity; tweakable.
  drop(1 +
         1.0 * xc[, 1] +
         0.5 * (xc[, 2]^2) +
         sin(xc[, 3]) +
         0.75 * (xc[, 4] * xc[, 5]))
}

fit_spline <- function(X_train, y_train, w_train = NULL, df = 4) {
  # X_train: matrix n x d (no intercept)
  d <- ncol(X_train)
  basis_list <- vector("list", d)
  B_train <- NULL
  for (j in 1:d) {
    bj <- splines::ns(X_train[, j], df = df)
    basis_list[[j]] <- bj
    B_train <- if (is.null(B_train)) bj else cbind(B_train, bj)
  }
  Z <- cbind(1, B_train)
  if (!is.null(w_train)) {
    fit <- lm(y_train ~ 0 + Z, weights = w_train)
  } else {
    fit <- lm(y_train ~ 0 + Z)
  }
  list(coef = fit$coefficients, basis = basis_list)
}

pred_spline <- function(model, X_new) {
  d <- ncol(X_new)
  B_new <- NULL
  for (j in 1:d) {
    bj <- predict(model$basis[[j]], newx = X_new[, j])
    B_new <- if (is.null(B_new)) bj else cbind(B_new, bj)
  }
  Z <- cbind(1, B_new)
  drop(Z %*% model$coef)
}

fit_rf <- function(X_train, y_train, w_train = NULL, num.trees = 200, mtry = NULL, min.node.size = 5) {
  # Prefer ranger (fast; supports case.weights). Fall back to randomForest (no weights).
  d <- ncol(X_train)
  if (is.null(mtry) || !is.finite(mtry)) mtry <- max(1L, floor(sqrt(d)))
  df_train <- as.data.frame(X_train)
  df_train$y <- y_train

  if (requireNamespace("ranger", quietly = TRUE)) {
    fit <- ranger::ranger(
      y ~ .,
      data = df_train,
      num.trees = num.trees,
      mtry = mtry,
      min.node.size = min.node.size,
      case.weights = w_train,
      respect.unordered.factors = "order",
      seed = 1
    )
    list(engine = "ranger", fit = fit)
  } else if (requireNamespace("randomForest", quietly = TRUE)) {
    if (!is.null(w_train)) {
      warning("randomForest does not support case weights for regression; fitting unweighted RF.")
    }
    fit <- randomForest::randomForest(
      x = df_train[, setdiff(names(df_train), "y"), drop = FALSE],
      y = df_train$y,
      ntree = num.trees,
      mtry = mtry
    )
    list(engine = "randomForest", fit = fit)
  } else {
    stop("Random forest requested, but neither 'ranger' nor 'randomForest' is installed.")
  }
}

pred_rf <- function(model, X_new) {
  df_new <- as.data.frame(X_new)
  if (identical(model$engine, "ranger")) {
    drop(predict(model$fit, data = df_new)$predictions)
  } else {
    drop(stats::predict(model$fit, newdata = df_new))
  }
}

greg_from_preds <- function(m_pop, m_s, y_s, d1_s) {
  sum(m_pop) + sum((y_s - m_s) * d1_s)
}

sreg_from_fit <- function(
  fit_fun,
  pred_fun,
  X_cov,
  y,
  d1,
  Index,
  K,
  wls,
  sampling,
  pi = NULL,
  n_h = NULL, st = NULL, cumn_h = NULL, len = NULL,
  fit_args = list(),
  pred_args = list()
) {
  perm <- sample.int(nrow(X_cov))
  groups <- cut(seq_along(perm), breaks = K, labels = FALSE)
  SubIndex_list <- split(perm, groups)

  y_parts <- numeric(K)
  e_oof <- numeric(length(Index)); names(e_oof) <- as.character(Index)

  for (k in 1:K) {
    SubIndex <- SubIndex_list[[k]]
    Index_sub1 <- Index[Index %in% SubIndex]
    Index_sub2 <- Index[!(Index %in% SubIndex)]

    X_tr <- X_cov[Index_sub2, , drop = FALSE]
    y_tr <- y[Index_sub2]
    w_tr <- if (wls) d1[Index_sub2] else NULL

    model <- do.call(fit_fun, c(list(X_train = X_tr, y_train = y_tr, w_train = w_tr), fit_args))

    # Fold population prediction sum
    m_fold_pop <- do.call(pred_fun, c(list(model = model, X_new = X_cov[SubIndex, , drop = FALSE]), pred_args))
    # Fold sampled prediction for residual correction
    m_fold_s <- do.call(pred_fun, c(list(model = model, X_new = X_cov[Index_sub1, , drop = FALSE]), pred_args))

    y_parts[k] <- greg_from_preds(m_fold_pop, m_fold_s, y[Index_sub1], d1[Index_sub1])
    e_oof[as.character(Index_sub1)] <- y[Index_sub1] - m_fold_s
  }

  y_sreg <- sum(y_parts)
  e_s_ss <- e_oof[as.character(Index)]

  if (sampling == "stratified") {
    sigma <- sqrt(quad_Omega_strata(e_s_ss, n_h, len, st, cumn_h))
  } else {
    sigma <- sqrt(var_poisson_approx(e_s_ss, pi[Index]))
  }

  list(est = y_sreg, sigma = sigma)
}

run_nonlinear <- function(
  SIMNUM,
  sampling,
  p_use,
  r_use,
  K,
  wls,
  np_methods,
  np_dim,
  spline_df,
  rf_trees,
  rf_mtry,
  rf_min_node,
  workers_override = NA,
  workers_cap = 80L,
  cluster_type = "PSOCK",
  cluster_outfile = "",
  seed_pop,
  seed_sim
) {
  stopifnot(sampling %in% c("stratified", "rejective"))

  set.seed(seed_pop)
  X0 <- mvtnorm::rmvnorm(
    n = N,
    mean = rep(2, p_use),
    sigma = ar1_cor(p_use, 0.2)
  )
  e <- rnorm(N, 0, 1)
  mu <- mu_nonlinear(X0, s)
  y <- drop(mu + e)
  t_y <- sum(y)

  len <- round(N / 4)
  z <- r_use * e + rnorm(N, 0, sqrt(1 - r_use^2))
  z_sorted_idx <- order(z)

  # Feature matrices
  X_lin <- cbind(1, X0)                           # linear uses intercept + all p
  X_np  <- X0[, 1:np_dim, drop = FALSE]           # spline/RF use first np_dim covariates (no intercept)

  # Which nonparam methods to run
  np_methods <- intersect(np_methods, c("spline", "rf"))

  workers <- if (!is.na(workers_override)) {
    max(1L, as.integer(workers_override))
  } else {
    get_workers(SIMNUM, cap = workers_cap)
  }
  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

  cl <- make_cluster_safe(workers, type = cluster_type, outfile = cluster_outfile)
  doParallel::registerDoParallel(cl)
  on.exit({ try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)

  message(sprintf("[%s/nonlinear] SIMNUM=%d, workers=%d, K=%d, fit=%s, np={%s}",
                  sampling, SIMNUM, workers, K, if (wls) "WLS" else "OLS", paste(np_methods, collapse = ",")))

  doRNG::registerDoRNG(seed = seed_sim + 2000L)

  export_vars <- c(
    # sampling
    "sample_design_stratified", "quad_Omega_strata",
    "sample_design_rejective", "make_pi_from_z", "var_poisson_approx",
    # nonlinear helpers
    "mu_nonlinear", "fit_spline", "pred_spline", "fit_rf", "pred_rf",
    "greg_from_preds", "sreg_from_fit",
    # objects
    "sampling", "N", "n", "s", "len", "z_sorted_idx", "z", "y", "mu", "t_y", "K", "wls",
    "X_lin", "X_np", "np_methods", "np_dim", "spline_df", "rf_trees", "rf_mtry", "rf_min_node"
  )

  final_res <- foreach::foreach(
    sim = 1:SIMNUM,
    .export = export_vars,
    .packages = character(0),
    .errorhandling = "pass"
  ) %dopar% {
    # Draw sample
    if (sampling == "stratified") {
      samp <- sample_design_stratified(N = N, n = n, len = len, z_sorted_idx = z_sorted_idx, y = y, X = X_lin)
      Index <- samp$Index
      d1 <- samp$d1
      d1_s <- samp$d1_s
      y_s <- samp$y_s
      # for variance
      n_h <- samp$n_h; st <- samp$st; cumn_h <- samp$cumn_h
      pi <- NULL
    } else {
      samp <- sample_design_rejective(N = N, n = n, z = z, y = y, X = X_lin, gamma = 1)
      Index <- samp$Index
      d1 <- samp$d1
      d1_s <- samp$d1_s
      y_s <- samp$y_s
      n_h <- st <- cumn_h <- NULL
      pi <- samp$pi
    }

    # Oracle difference
    y_diff <- drop(sum(mu) + sum((y_s - mu[Index]) * d1_s))
    if (sampling == "stratified") {
      sig_diff <- sqrt(quad_Omega_strata(y_s - mu[Index], n_h, len, st, cumn_h))
    } else {
      sig_diff <- sqrt(var_poisson_approx(y_s - mu[Index], pi[Index]))
    }

    out <- list(Diff = y_diff)
    sig <- list(Diff = sig_diff)

    # --- Linear working model (lm)
    Xs_lin <- X_lin[Index, , drop = FALSE]
    if (wls) {
      lm_fit <- lm(y_s ~ 0 + Xs_lin, weights = d1_s)
    } else {
      lm_fit <- lm(y_s ~ 0 + Xs_lin)
    }
    beta_hat <- lm_fit$coefficients
    m_pop <- drop(X_lin %*% beta_hat)
    m_s   <- drop(Xs_lin %*% beta_hat)
    out$GREG.linear <- greg_from_preds(m_pop, m_s, y_s, d1_s)

    e_s <- y_s - m_s
    if (sampling == "stratified") {
      sig$GREG.linear <- sqrt(quad_Omega_strata(e_s, n_h, len, st, cumn_h))
    } else {
      sig$GREG.linear <- sqrt(var_poisson_approx(e_s, pi[Index]))
    }

    # SREG.linear
    fit_lm <- function(X_train, y_train, w_train = NULL) {
      if (!is.null(w_train)) {
        fit <- lm(y_train ~ 0 + X_train, weights = w_train)
      } else {
        fit <- lm(y_train ~ 0 + X_train)
      }
      list(beta = fit$coefficients)
    }
    pred_lm <- function(model, X_new) drop(X_new %*% model$beta)

    sreg_lin <- sreg_from_fit(
      fit_fun = fit_lm,
      pred_fun = pred_lm,
      X_cov = X_lin,
      y = y,
      d1 = d1,
      Index = Index,
      K = K,
      wls = wls,
      sampling = sampling,
      pi = pi,
      n_h = n_h, st = st, cumn_h = cumn_h, len = len
    )
    out$SREG.linear <- sreg_lin$est
    sig$SREG.linear <- sreg_lin$sigma

    # --- Spline working model
    if ("spline" %in% np_methods) {
      Xs_np <- X_np[Index, , drop = FALSE]
      w_fit <- if (wls) d1_s else NULL
      sp_fit <- fit_spline(X_train = Xs_np, y_train = y_s, w_train = w_fit, df = spline_df)
      m_pop_sp <- pred_spline(sp_fit, X_np)
      m_s_sp   <- pred_spline(sp_fit, Xs_np)
      out$GREG.spline <- greg_from_preds(m_pop_sp, m_s_sp, y_s, d1_s)

      e_s_sp <- y_s - m_s_sp
      if (sampling == "stratified") {
        sig$GREG.spline <- sqrt(quad_Omega_strata(e_s_sp, n_h, len, st, cumn_h))
      } else {
        sig$GREG.spline <- sqrt(var_poisson_approx(e_s_sp, pi[Index]))
      }

      sreg_sp <- sreg_from_fit(
        fit_fun = fit_spline,
        pred_fun = pred_spline,
        X_cov = X_np,
        y = y,
        d1 = d1,
        Index = Index,
        K = K,
        wls = wls,
        sampling = sampling,
        pi = pi,
        n_h = n_h, st = st, cumn_h = cumn_h, len = len,
        fit_args = list(df = spline_df)
      )
      out$SREG.spline <- sreg_sp$est
      sig$SREG.spline <- sreg_sp$sigma
    }

    # --- Random forest working model
    if ("rf" %in% np_methods) {
      Xs_np <- X_np[Index, , drop = FALSE]
      w_fit <- if (wls) d1_s else NULL

      rf_fit <- fit_rf(X_train = Xs_np, y_train = y_s, w_train = w_fit, num.trees = rf_trees, mtry = rf_mtry, min.node.size = rf_min_node)
      m_pop_rf <- pred_rf(rf_fit, X_np)
      m_s_rf   <- pred_rf(rf_fit, Xs_np)
      out$GREG.rf <- greg_from_preds(m_pop_rf, m_s_rf, y_s, d1_s)

      e_s_rf <- y_s - m_s_rf
      if (sampling == "stratified") {
        sig$GREG.rf <- sqrt(quad_Omega_strata(e_s_rf, n_h, len, st, cumn_h))
      } else {
        sig$GREG.rf <- sqrt(var_poisson_approx(e_s_rf, pi[Index]))
      }

      sreg_rf <- sreg_from_fit(
        fit_fun = fit_rf,
        pred_fun = pred_rf,
        X_cov = X_np,
        y = y,
        d1 = d1,
        Index = Index,
        K = K,
        wls = wls,
        sampling = sampling,
        pi = pi,
        n_h = n_h, st = st, cumn_h = cumn_h, len = len,
        fit_args = list(num.trees = rf_trees, mtry = rf_mtry, min.node.size = rf_min_node)
      )
      out$SREG.rf <- sreg_rf$est
      sig$SREG.rf <- sreg_rf$sigma
    }

    out <- unlist(out)
    sig <- unlist(sig)

    list(
      y = out,
      sigma = sig,
      cover = ifelse(abs(out - t_y) > 1.96 * sig, 0, 1)
    )
  }

  ok <- vapply(final_res, function(x) is.list(x) && is.numeric(unlist(x$y)), logical(1))
  if (sum(ok) == 0) stop("All simulations failed in nonlinear run.")

  ys <- do.call(rbind, lapply(final_res[ok], `[[`, "y"))
  bias <- colMeans(ys - t_y)
  rmse <- apply(ys - t_y, 2, function(x) sqrt(mean(x^2)))

  list(sampling = sampling, wls = wls, bias = bias, rmse = rmse, err = ys - t_y)
}

make_nonlinear_summary_plot <- function(resNL, show_legend = TRUE) {
  df <- rbind(
    data.frame(estimator = names(resNL$rmse), metric = "RMSE", value = unname(resNL$rmse)),
    data.frame(estimator = names(resNL$bias), metric = "Bias", value = unname(resNL$bias))
  )
  # Put oracle Diff first
  ord <- unique(c("Diff", sort(setdiff(unique(df$estimator), "Diff"))))
  df$estimator <- factor(df$estimator, levels = ord)

  gg <- ggplot(df, aes(x = estimator, y = value, group = metric)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_point(size = 2) +
    facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  if (!show_legend) gg <- gg + theme(legend.position = "none")
  gg
}




make_nonlinear_boxplot <- function(res_list, design_labels = NULL, base_size = 16, facet_scales = "fixed") {
  # res_list: named list of run_nonlinear() results; each must contain $err (SIMNUM x est)
  if (is.null(design_labels)) {
    design_labels <- names(res_list)
  }
  if (is.null(design_labels) || !length(design_labels)) {
    design_labels <- rep("design", length(res_list))
  }

  build_df <- function(res, design) {
    err <- res$err
    if (is.null(dim(err))) {
      err <- matrix(err, ncol = 1)
      colnames(err) <- names(res$bias)
    }
    data.frame(
      design = design,
      estimator = rep(colnames(err), each = nrow(err)),
      error = as.vector(err),
      stringsAsFactors = FALSE
    )
  }

  df <- do.call(rbind, Map(build_df, res_list, design_labels))

  # Ordering / nicer labels
  est_order <- c(
    "Diff",
    "GREG.linear", "SREG.linear",
    "GREG.spline", "SREG.spline",
    "GREG.rf", "SREG.rf"
  )
  est_order <- est_order[est_order %in% unique(df$estimator)]

  pretty <- c(
    Diff = "Diff",
    `GREG.linear` = "GREG (linear)",
    `SREG.linear` = "SREG (linear)",
    `GREG.spline` = "GREG (spline)",
    `SREG.spline` = "SREG (spline)",
    `GREG.rf`     = "GREG (RF)",
    `SREG.rf`     = "SREG (RF)"
  )

  df$design <- factor(df$design, levels = c("stratified", "rejective"))
  df$estimator <- factor(df$estimator, levels = est_order, labels = unname(pretty[est_order]))

  ggplot(df, aes(x = estimator, y = error)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey35") +
    geom_boxplot(outlier.size = 0.7, linewidth = 0.45) +
    facet_wrap(~ design, nrow = 1, scales = facet_scales) +
    labs(x = NULL, y = expression(hat(T) - T)) +
    theme_bw(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1),
      strip.text  = element_text(size = base_size * 1.0)
    )
}
## =========================================================
## Main
## =========================================================

args <- commandArgs(trailingOnly = TRUE)
opt <- parse_kv_args(args)

TASK    <- tolower(opt$task   %||% DEFAULTS$task)
FIT     <- tolower(opt$fit    %||% DEFAULTS$fit)
SIMNUM  <- as.integer(opt$simnum   %||% DEFAULTS$simnum)
DESIGN  <- tolower(opt$design %||% DEFAULTS$design)
VARY    <- tolower(opt$vary   %||% DEFAULTS$vary)
OUTDIR  <- opt$outdir %||% DEFAULTS$outdir
BASE_SIZE <- as.numeric(opt$base_size %||% DEFAULTS$base_size)
K       <- as.integer(opt$K %||% DEFAULTS$K)
K_VEC   <- parse_csv_int(opt$K_grid %||% DEFAULTS$K_grid, default = parse_csv_int(DEFAULTS$K_grid, default = c(2,5,10,20)))
WORKERS <- suppressWarnings(as.integer(opt$workers %||% DEFAULTS$workers))
WORKERS_CAP <- as.integer(opt$workers_cap %||% DEFAULTS$workers_cap)
CLUSTER <- tolower(opt$cluster %||% DEFAULTS$cluster)
CLUSTER_OUTFILE <- opt$cluster_outfile %||% DEFAULTS$cluster_outfile
CLUSTER_TYPE <- switch(CLUSTER,
  "psock" = "PSOCK",
  "fork"  = "FORK",
  "auto"  = "AUTO",
  stop("Unknown --cluster=", CLUSTER)
)
SEED_POP <- as.integer(opt$seed_pop %||% DEFAULTS$seed_pop)
SEED_SIM <- as.integer(opt$seed_sim %||% DEFAULTS$seed_sim)

NP_METHODS <- parse_csv_chr(opt$np_methods %||% DEFAULTS$np_methods, default = c("spline","rf"))
NP_DIM     <- as.integer(opt$np_dim %||% DEFAULTS$np_dim)
SPLINE_DF  <- as.integer(opt$spline_df %||% DEFAULTS$spline_df)
RF_TREES   <- as.integer(opt$rf_trees %||% DEFAULTS$rf_trees)
RF_MTRY    <- suppressWarnings(as.integer(opt$rf_mtry %||% DEFAULTS$rf_mtry))
RF_MINNODE <- as.integer(opt$rf_min_node %||% DEFAULTS$rf_min_node)

if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

designs <- switch(DESIGN,
  "all" = c("stratified", "rejective"),
  "stratified" = "stratified",
  "rejective"  = "rejective",
  stop("Unknown --design=", DESIGN)
)

vary_list <- switch(VARY,
  "all" = c("p", "r"),
  "p" = "p",
  "r" = "r",
  stop("Unknown --vary=", VARY)
)

wls_flags <- switch(FIT,
  "ols" = c(FALSE),
  "wls" = c(TRUE),
  "both" = c(FALSE, TRUE),
  stop("Unknown --fit=", FIT)
)

run_grid_task <- function() {
  # We run p- and r-grids and additionally create a paired figure with
  # varying p (left) and varying r (right) for each (design, fit) case.
  for (sampling in designs) {
    for (wls in wls_flags) {
      res_by_vary <- list()

      for (vary_over in vary_list) {
        res <- run_grid(
          SIMNUM = SIMNUM,
          sampling = sampling,
          vary_over = vary_over,
          K = K,
          wls = wls,
          workers_override = WORKERS,
          workers_cap = WORKERS_CAP,
          cluster_type = CLUSTER_TYPE,
          cluster_outfile = CLUSTER_OUTFILE,
          seed_pop = SEED_POP,
          seed_sim = SEED_SIM
        )

        res_by_vary[[vary_over]] <- res

        # Save results
        out_rds <- file.path(OUTDIR, sprintf("results_%s_%s%s.rds", sampling, vary_over, if (wls) "_wls" else ""))
        saveRDS(res, out_rds)

        # Save the single-panel version (kept for backward compatibility)
        suffix_design <- if (sampling == "stratified") "" else "2"
        suffix_wls <- if (wls) "_wls" else ""
        fname <- sprintf("plot_%s%s%s.png", vary_over, suffix_design, suffix_wls)

        p1 <- make_bias_rmse_plot(res, show_legend = FALSE, base_size = BASE_SIZE)
        save_png(p1, file.path(OUTDIR, fname), width = 7.6, height = 4.1)
        message(sprintf("Wrote %s", file.path(OUTDIR, fname)))
      }

      # Paired (p | r) figure: requested SM layout
      if (all(c("p", "r") %in% names(res_by_vary))) {
        p_left  <- make_bias_rmse_plot(res_by_vary[["p"]], show_legend = FALSE, base_size = BASE_SIZE)
        p_right <- make_bias_rmse_plot(res_by_vary[["r"]], show_legend = TRUE,  base_size = BASE_SIZE, legend_position = "right")

        # Slight spacing tweak so two panels don't collide
        p_left  <- p_left  + theme(plot.margin = margin(5.5, 3, 5.5, 5.5))
        p_right <- p_right + theme(plot.margin = margin(5.5, 5.5, 5.5, 3))

        fit_tag <- if (wls) "wls" else "ols"
        fname_pr <- sprintf("plot_pr_%s_%s.png", sampling, fit_tag)
        save_two_panel_png(p_left, p_right, file.path(OUTDIR, fname_pr), width = 12.0, height = 4.4, dpi = 300)
        message(sprintf("Wrote %s", file.path(OUTDIR, fname_pr)))
      }
    }
  }
}

run_K_task <- function() {
  for (wls in wls_flags) {
    res_list <- list()

    for (sampling in designs) {
      resK <- run_K_grid(
        SIMNUM = SIMNUM,
        sampling = sampling,
        K_vec = K_VEC,
        p_use = p_fixed,
        r_use = r_fixed,
        wls = wls,
        workers_override = WORKERS,
        workers_cap = WORKERS_CAP,
        cluster_type = CLUSTER_TYPE,
        cluster_outfile = CLUSTER_OUTFILE,
        seed_pop = SEED_POP,
        seed_sim = SEED_SIM
      )
      res_list[[sampling]] <- resK

      suffix_wls <- if (wls) "_wls" else ""
      fname <- sprintf("supp_K_%s%s.png", sampling, suffix_wls)
      pK <- make_K_plot(resK, show_legend = TRUE, base_size = BASE_SIZE)
      save_png(pK, file.path(OUTDIR, fname), width = 7.4, height = 4.2)
      message(sprintf("Wrote %s", file.path(OUTDIR, fname)))
    }

    # Combined design figure (when both designs are requested)
    if (length(res_list) >= 2) {
      suffix_wls <- if (wls) "_wls" else ""
      fname2 <- sprintf("supp_K_both%s.png", suffix_wls)
      pK2 <- make_K_plot_combined(res_list, base_size = BASE_SIZE)
      save_png(pK2, file.path(OUTDIR, fname2), width = 11.0, height = 4.2)
      message(sprintf("Wrote %s", file.path(OUTDIR, fname2)))
    }
  }
}

run_nonlinear_task <- function() {
  for (wls in wls_flags) {
    res_list <- list()

    for (sampling in designs) {
      resNL <- run_nonlinear(
        SIMNUM = SIMNUM,
        sampling = sampling,
        p_use = p_fixed,
        r_use = r_fixed,
        K = K,
        wls = wls,
        np_methods = NP_METHODS,
        np_dim = NP_DIM,
        spline_df = SPLINE_DF,
        rf_trees = RF_TREES,
        rf_mtry = RF_MTRY,
        rf_min_node = RF_MINNODE,
        workers_override = WORKERS,
        workers_cap = WORKERS_CAP,
        cluster_type = CLUSTER_TYPE,
        cluster_outfile = CLUSTER_OUTFILE,
        seed_pop = SEED_POP,
        seed_sim = SEED_SIM
      )
      res_list[[sampling]] <- resNL

      suffix_wls <- if (wls) "_wls" else ""
      fname <- sprintf("supp_nonlinear_%s%s.png", sampling, suffix_wls)
      pNL <- make_nonlinear_boxplot(list(resNL), design_labels = sampling, base_size = BASE_SIZE)
      save_png(pNL, file.path(OUTDIR, fname), width = 7.8, height = 5.2)
      message(sprintf("Wrote %s", file.path(OUTDIR, fname)))
    }

    # Combined design figure (when both designs are requested)
    if (length(res_list) >= 2) {
      suffix_wls <- if (wls) "_wls" else ""
      fname2 <- sprintf("supp_nonlinear_both%s.png", suffix_wls)
      pNL2 <- make_nonlinear_boxplot(res_list, base_size = BASE_SIZE)
      save_png(pNL2, file.path(OUTDIR, fname2), width = 11.0, height = 5.2)
      message(sprintf("Wrote %s", file.path(OUTDIR, fname2)))
    }
  }
}

if (TASK %in% c("grid", "all")) run_grid_task()
if (TASK %in% c("k", "all")) run_K_task()
if (TASK %in% c("nonlinear", "all")) run_nonlinear_task()

message("Done.")
