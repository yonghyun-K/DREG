# Usage (batch):
#   Rscript main.R --run --B <SIMNUM> [--N <pop>] [--n <samp>] \
#                  [--sampling stratified|rejective] [--vary_over p|r] \
#                  [--p <num|list>] [--r <num|list>] \
#                  [--wls] [--an true|false] [--lasso_only] [--K <int>] \
#                  [--cores <int>] [--out <dir>]
#
# Notes:
#   --run         acts as a "do it" flag (any value accepted; presence is enough).
#   --B           number of Monte Carlo replicates (required in batch mode).
#   --N, --n      override population and sample size (defaults: 1000, 300).
#   --p, --r      override sweep grid / pinned value; accepts scalar or "a,b,c".
#                 When --vary_over=p, --p is the sweep grid and --r is the pinned
#                 value (must be scalar); vice versa for --vary_over=r.
#   --wls         use weighted least squares (weights = 1/pi) for the regression
#                 fits. Bare flag means TRUE; you can also write --wls true/false.
#                 Default is OLS (unweighted fits).
#   --an          also compute the conditional cross-fitted estimator of
#                 An, Dagdoug, Haziza and Tille (2026). Default TRUE.
#   --lasso_only  skip all least-squares fits on the full covariate set (for p >= n).
#   --K           number of cross-fitting folds (default 10).
#   --cores       number of parallel workers; overrides SLURM / detectCores().
#   --out         output directory name; defaults to a timestamped folder.
#
# Revision (JSSAM-2026-121):
#   * Lasso: glmnet's own intercept is added to the coefficient of the column of
#     ones, so lasso predictions include the intercept.
#   * Rejective sampling: exact second-order inclusion probabilities
#     (sampling::UPmaxentropypi2) are used in the Horvitz-Thompson variance
#     estimator (4.11) for every estimator.
#   * New estimators: Diff.star (oracle built on the target of the fitted
#     regression) and AnCF / AnCF.Lasso (An et al., 2026).
#   * Output: BIAS, SD, RMSE, RB, CR, Monte Carlo SE of CR and mean CI length per
#     grid value (summary_<value>.csv) and all replicate draws (replicates.rds).
#   * Seeds: population set.seed(2); replicates doRNG seed 11 (reset per grid value).

## -------- CLI parser (flag-style) --------
.parse_flags <- function(args) {
  out <- list()
  i <- 1L
  # Flags that may appear without an argument (presence == TRUE).
  bool_flags <- c("run", "wls", "an", "lasso_only")
  while (i <= length(args)) {
    a <- args[i]
    if (!startsWith(a, "--")) stop("Expected --flag, got: ", a)
    key <- sub("^--", "", a)
    if (key %in% bool_flags) {
      nxt <- if (i + 1L <= length(args)) args[i + 1L] else NA_character_
      if (is.na(nxt) || startsWith(nxt, "--")) {
        out[[key]] <- "true"; i <- i + 1L
      } else {
        out[[key]] <- nxt; i <- i + 2L
      }
    } else {
      if (i + 1L > length(args)) stop("Missing value for --", key)
      out[[key]] <- args[i + 1L]; i <- i + 2L
    }
  }
  out
}

.parse_num_list <- function(s) {
  if (is.null(s)) return(NULL)
  if (!is.character(s)) {
    stop(".parse_num_list expected character, got ", class(s)[1],
         " (value: ", paste(s, collapse = ", "), ").")
  }
  toks <- strsplit(s, "[,[:space:]]+", perl = TRUE)[[1]]
  toks <- toks[nzchar(toks)]
  if (length(toks) == 0L) return(NULL)
  v <- suppressWarnings(as.numeric(toks))
  if (any(is.na(v))) stop("Could not parse numeric list from '", s, "'.")
  v
}

.parse_bool <- function(x, default) {
  if (is.null(x)) return(default)
  v <- tolower(as.character(x))
  if (v %in% c("true", "t", "yes", "y", "1", "")) return(TRUE)
  if (v %in% c("false", "f", "no", "n", "0")) return(FALSE)
  stop("Could not parse boolean value '", x, "'.")
}

if (!interactive()) {
  opts <- .parse_flags(commandArgs(trailingOnly = TRUE))
  # NOTE: always use opts[["x"]], NEVER opts$x (partial matching).
  if (is.null(opts[["B"]])) stop("Must supply --B <SIMNUM>.")
  SIMNUM <- as.numeric(opts[["B"]])
  if (!is.finite(SIMNUM) || SIMNUM <= 0) stop("--B must be a positive number.")

  sampling  <- if (!is.null(opts[["sampling"]]))  tolower(opts[["sampling"]])  else "stratified"
  vary_over <- if (!is.null(opts[["vary_over"]])) tolower(opts[["vary_over"]]) else "p"

  N_user    <- if (!is.null(opts[["N"]]))     as.integer(opts[["N"]])     else 1000L
  n_user    <- if (!is.null(opts[["n"]]))     as.integer(opts[["n"]])     else 300L
  cores_req <- if (!is.null(opts[["cores"]])) as.integer(opts[["cores"]]) else NA_integer_
  K_user    <- if (!is.null(opts[["K"]]))     as.integer(opts[["K"]])     else 10L
  out_dir   <- opts[["out"]]

  p_user <- .parse_num_list(opts[["p"]])
  r_user <- .parse_num_list(opts[["r"]])

  wls_use    <- .parse_bool(opts[["wls"]], FALSE)
  an_use     <- .parse_bool(opts[["an"]], TRUE)
  lasso_only <- .parse_bool(opts[["lasso_only"]], FALSE)
} else {
  SIMNUM    <- 20
  sampling  <- "stratified"      # or "rejective"
  vary_over <- "r"               # "p" or "r"
  N_user    <- 1000L
  n_user    <- 300L
  cores_req <- NA_integer_
  K_user    <- 10L
  out_dir   <- NULL
  p_user    <- NULL
  r_user    <- -0.75
  wls_use   <- FALSE
  an_use    <- TRUE
  lasso_only <- FALSE
}

if (!sampling %in% c("stratified", "rejective")) {
  stop("Unknown sampling='", sampling, "'. Use 'stratified' or 'rejective'.")
}
if (!vary_over %in% c("p", "r")) {
  stop("Unknown vary_over='", vary_over, "'. Use 'p' or 'r'.")
}
message("Sampling: ", sampling, " | varying over: ", vary_over)
message("Regression fit: ", if (wls_use) "WLS (weighted by 1/pi)" else "OLS",
        " | An et al. estimator: ", an_use, " | lasso only: ", lasso_only, " | K = ", K_user)

timenow1 <- Sys.time()
timenow0 <- gsub(" ", "_", gsub("[-:]", "", timenow1))
timenow0 <- paste0(c(sampling, vary_over, timenow0), collapse = "_")
if (!is.null(out_dir)) timenow0 <- out_dir
timenow  <- paste0(basename(timenow0), ".txt")

## -------- Packages --------
suppressPackageStartupMessages({
  library(xtable)
  library(mvtnorm)
  library(glmnet)
  library(ggplot2)
  library(foreach)
  library(doParallel)
  library(doRNG)
  library(sampling)
})

## -------- Parameters --------
N <- N_user
n <- n_user
s <- 5
K_use <- K_user

default_seq_p <- c(10, 50, 90, 130, 170, 210, 250)
default_seq_r <- c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)
default_p     <- 90
default_r     <- -0.75

if (vary_over == "p") {
  seq_p   <- if (!is.null(p_user)) p_user else default_seq_p
  if (!is.null(r_user) && length(r_user) > 1L) stop("--r must be a single value when --vary_over p")
  r_fixed <- if (!is.null(r_user)) r_user else default_r
  p_fixed <- NA_real_
  seq     <- seq_p
} else {
  if (!is.null(p_user) && length(p_user) > 1L) stop("--p must be a single value when --vary_over r")
  p_fixed <- if (!is.null(p_user)) as.integer(p_user) else default_p
  seq_r   <- if (!is.null(r_user)) r_user else default_seq_r
  r_fixed <- NA_real_
  seq     <- seq_r
}

message(sprintf("Sweep dimension: %s = {%s}", vary_over, paste(seq, collapse = ", ")))

## -------- Utilities --------
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
  rho ^ exponent
}

## =========================================================
## Sampling designs
## =========================================================

## ---- Stratified: four strata of equal size defined by sorted z ----
sample_design_stratified <- function(N, n, len, z_sorted_idx) {
  n_h <- as.integer(c(15, 20, 30, 35) / 100 * n)
  stopifnot(sum(n_h) == n)
  Index <- integer(n)
  d1 <- numeric(N)
  cumn_h <- cumsum(n_h)
  st <- c(1L, head(cumn_h, -1L) + 1L)
  for (i in 1:4) {
    Idx <- z_sorted_idx[(len * (i - 1) + 1):(len * i)]
    d1[Idx] <- len / n_h[i]
    Index[st[i]:cumn_h[i]] <- sample(Idx, size = n_h[i], replace = FALSE)
  }
  list(Index = Index, d1 = d1, n_h = n_h, st = st, cumn_h = cumn_h, len = len)
}

# Exact stratified-SRSWOR HT variance estimator (4.11); e ordered by strata blocks.
quad_Omega_strata <- function(e, n_h, len, st, cumn_h) {
  out <- 0
  for (i in seq_along(n_h)) {
    eh <- e[st[i]:cumn_h[i]]
    nh <- n_h[i]
    c_off <- len^2/nh^2 - len*(len-1)/nh/(nh-1)
    d_diag <- (1 - nh/len) * len^2/nh^2
    s1 <- sum(eh); s2 <- sum(eh*eh)
    out <- out + c_off * (s1*s1) + (d_diag - c_off) * s2
  }
  out
}

## ---- Rejective (conditional Poisson / maximum entropy) ----
make_pi_from_z <- function(z, n, gamma = 1) {
  m  <- gamma * 1 / (1 + exp(-z))
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
    pi[idx0] <- n_rem * m[idx0] / sum(m[idx0])
  }
  if (abs(sum(pi) - n) > 1e-8) {
    idx0 <- which(pi < 1)
    pi[idx0] <- pi[idx0] * (n - sum(pi[pi == 1])) / sum(pi[idx0])
  }
  pmax(pmin(pi, 1), 1e-12)
}

# Maximum-entropy (conditional Poisson) design with prescribed first-order
# probabilities pik (Chen, Dempster and Liu, 1994). The working probabilities and the
# matrix qmat = UPMEqfromw(w, n) are computed once per population; drawing with
# UPMEsfromq(qmat) is then identical to sampling::UPmaxentropy(pik) when all pik < 1.
sample_design_rejective <- function(pik, qmat) {
  Index <- which(sampling::UPMEsfromq(qmat) == 1)
  list(Index = Index, d1 = 1 / pik)
}

# HT variance estimator (4.11) with exact pi_ij: Ds = (pi_ij - pi_i pi_j) / pi_ij on the sample.
var_ht_exact <- function(e_s, pi_s, Ds) {
  a <- e_s / pi_s
  drop(crossprod(a, Ds %*% a))
}

## =========================================================
## Fitting helpers
## =========================================================
fit_ls <- function(Xs, ys, w = NULL) {
  b <- if (is.null(w)) stats::lm.fit(Xs, ys)$coefficients else stats::lm.wfit(Xs, ys, w)$coefficients
  b[is.na(b)] <- 0
  b
}

# Lasso with 10-fold CV (lambda.min, standardize = TRUE). X has a column of ones in
# position 1; glmnet fits its own (unpenalized) intercept, which is added to that
# column's coefficient so that X %*% beta reproduces the lasso predictions.
fit_lasso_beta <- function(Xs, ys, w = NULL) {
  cvm <- if (is.null(w)) glmnet::cv.glmnet(Xs, ys) else glmnet::cv.glmnet(Xs, ys, weights = w)
  b <- as.vector(stats::coef(cvm, s = "lambda.min"))
  beta <- b[-1]
  beta[1] <- beta[1] + b[1]
  beta
}

# Out-of-fold predictions: mhat[i] = prediction for unit i from the fit on A \ A_{kappa_i}.
cf_predict <- function(kappa, K, Index, X, y, d1, fit_fun, wls) {
  mhat <- numeric(nrow(X))
  for (k in 1:K) {
    tr <- Index[kappa[Index] != k]
    te <- which(kappa == k)
    b <- fit_fun(X[tr, , drop = FALSE], y[tr], if (wls) d1[tr] else NULL)
    mhat[te] <- drop(X[te, , drop = FALSE] %*% b)
  }
  mhat
}

# Balanced random partition of U into K folds, independent of the design.
make_pop_folds <- function(N, K) {
  kappa <- integer(N)
  kappa[sample.int(N)] <- cut(seq_len(N), breaks = K, labels = FALSE)
  kappa
}

## ---- An et al. (2026): conditional cross-fitted estimator ----
# Stratified SRSWOR: design-aware folds within strata (their Algorithm 3 with the
# split of Lu, Shi, Liu and Ding, 2025): in each stratum the sampled and the
# nonsampled units are split separately into K folds of (nearly) equal size.
make_folds_strat_aware <- function(strat_id, Index, K) {
  kappa <- integer(length(strat_id))
  inS <- logical(length(strat_id)); inS[Index] <- TRUE
  for (h in sort(unique(strat_id))) {
    Uh <- which(strat_id == h)
    Ah <- Uh[inS[Uh]]; Bh <- Uh[!inS[Uh]]
    Nh <- length(Uh); nh <- length(Ah)
    o <- sample.int(K)
    Nhk <- (rep(Nh %/% K, K) + (seq_len(K) <= Nh %% K))[o]
    nhk <- (rep(nh %/% K, K) + (seq_len(K) <= nh %% K))[o]
    kappa[Ah[sample.int(length(Ah))]] <- rep(seq_len(K), times = nhk)
    kappa[Bh[sample.int(length(Bh))]] <- rep(seq_len(K), times = Nhk - nhk)
  }
  kappa
}

# Estimator and variance with weights 1/pi_{i|A} = N_hk / n_hk; given the fold sample
# sizes, each stratum-by-fold cell is an SRSWOR of size n_hk from N_hk units.
an_strat <- function(mhat, kappa, strat_id, Index, y) {
  cell_U <- paste(strat_id, kappa)
  cell_A <- cell_U[Index]
  N_c <- table(cell_U); n_c <- table(cell_A)
  e <- y[Index] - mhat[Index]
  Nk <- as.numeric(N_c[cell_A]); nk <- as.numeric(n_c[cell_A])
  total <- sum(mhat) + sum(e * Nk / nk)
  V <- 0
  for (cc in names(n_c)) {
    ec <- e[cell_A == cc]
    Ncc <- as.numeric(N_c[cc]); ncc <- length(ec)
    if (ncc > 1) V <- V + Ncc^2 * (1 - ncc / Ncc) * stats::var(ec) / ncc
  }
  c(total, sqrt(max(V, 0)))
}

# Rejective sampling: same population folds as DREG; given the fold sample sizes n_k,
# the sample in fold k is conditional Poisson of size n_k on U_k with the same
# canonical weights wcan, so pi_{i|A} and pi_{ij|A} follow from the CPS recursions.
an_rej_cond <- function(kappa, K, Index, wcan) {
  lapply(1:K, function(k) {
    Uk <- which(kappa == k)
    Ak <- Index[kappa[Index] == k]
    pc <- sampling::UPMEpikfromq(sampling::UPMEqfromw(wcan[Uk], length(Ak)))
    p2 <- sampling::UPMEpik2frompikw(pc, wcan[Uk])
    pos <- match(Ak, Uk)
    pcA <- pc[pos]
    D <- (p2[pos, pos, drop = FALSE] - tcrossprod(pcA)) / p2[pos, pos, drop = FALSE]
    list(Ak = Ak, pcA = pcA, D = D)
  })
}

an_rej <- function(mhat, cond, y) {
  total <- sum(mhat); V <- 0
  for (cc in cond) {
    a <- (y[cc$Ak] - mhat[cc$Ak]) / cc$pcA
    total <- total + sum(a)
    V <- V + drop(crossprod(a, cc$D %*% a))
  }
  c(total, sqrt(max(V, 0)))
}

## =========================================================
## Estimation / inference for one sample
## =========================================================
run_estimators <- function(sampling, N, n, s, y, X, mu, mu_star, samp, K,
                           wls, an, lasso_only, pik = NULL, Dmat = NULL,
                           wcan = NULL, strat_id = NULL) {
  Index <- samp$Index
  d1 <- samp$d1
  d1_s <- d1[Index]; y_s <- y[Index]; X_s <- X[Index, , drop = FALSE]
  w_fit <- if (wls) d1_s else NULL

  if (sampling == "stratified") {
    vfun <- function(e_s) quad_Omega_strata(e_s, samp$n_h, samp$len, samp$st, samp$cumn_h)
  } else {
    Ds <- Dmat[Index, Index]
    pis <- pik[Index]
    vfun <- function(e_s) var_ht_exact(e_s, pis, Ds)
  }

  est <- c(); sig <- c()
  add <- function(name, mpop) {
    e_s <- y_s - mpop[Index]
    est[name] <<- sum(mpop) + sum(e_s * d1_s)
    sig[name] <<- sqrt(max(vfun(e_s), 0))
  }

  ## HT and oracles
  est["HT"] <- sum(y_s * d1_s); sig["HT"] <- sqrt(max(vfun(y_s), 0))
  add("Diff", mu)
  add("Diff.star", mu_star)

  ## GREG on the true support (intercept + first s covariates)
  b_or <- fit_ls(X_s[, 1:(1 + s), drop = FALSE], y_s, w_fit)
  add("GREG.oracle", drop(X[, 1:(1 + s), drop = FALSE] %*% b_or))

  ## Full-sample fits
  if (!lasso_only) {
    b <- fit_ls(X_s, y_s, w_fit)
    add("GREG", drop(X %*% b))
  }
  bL <- fit_lasso_beta(X_s, y_s, w_fit)
  add("GREG.Lasso", drop(X %*% bL))

  ## DREG: population folds independent of the design
  kappa <- make_pop_folds(N, K)
  if (!lasso_only) {
    mhat <- cf_predict(kappa, K, Index, X, y, d1, fit_ls, wls)
    add("DREG", mhat)
  }
  mhatL <- cf_predict(kappa, K, Index, X, y, d1, fit_lasso_beta, wls)
  add("DREG.Lasso", mhatL)

  ## An et al. (2026)
  if (an) {
    if (sampling == "stratified") {
      kappa_aw <- make_folds_strat_aware(strat_id, Index, K)
      if (!lasso_only) {
        r <- an_strat(cf_predict(kappa_aw, K, Index, X, y, d1, fit_ls, wls), kappa_aw, strat_id, Index, y)
        est["AnCF"] <- r[1]; sig["AnCF"] <- r[2]
      }
      r <- an_strat(cf_predict(kappa_aw, K, Index, X, y, d1, fit_lasso_beta, wls), kappa_aw, strat_id, Index, y)
      est["AnCF.Lasso"] <- r[1]; sig["AnCF.Lasso"] <- r[2]
    } else {
      cond <- an_rej_cond(kappa, K, Index, wcan)
      if (!lasso_only) {
        r <- an_rej(mhat, cond, y); est["AnCF"] <- r[1]; sig["AnCF"] <- r[2]
      }
      r <- an_rej(mhatL, cond, y); est["AnCF.Lasso"] <- r[1]; sig["AnCF.Lasso"] <- r[2]
    }
  }

  ## Lasso selection (columns 2..(1+s) are the signals; column 1 is the intercept)
  FDR <- sum(bL[-(1:(1 + s))] != 0) / (ncol(X) - 1 - s)
  FNR <- sum(bL[2:(1 + s)] == 0) / s

  ord <- c("HT", "Diff", "Diff.star", "GREG.oracle", "GREG", "DREG", "AnCF",
           "GREG.Lasso", "DREG.Lasso", "AnCF.Lasso")
  ord <- ord[ord %in% names(est)]
  list(y_res = est[ord], sigma_res = sig[ord], FDR = FDR, FNR = FNR)
}

## =========================================================
## Parallel setup
## =========================================================
get_allocated_cores <- function() {
  candidates <- suppressWarnings(as.integer(c(Sys.getenv("SLURM_CPUS_PER_TASK"),
                                              Sys.getenv("PBS_NP"), Sys.getenv("LSB_DJOB_NUMPROC"))))
  candidates <- candidates[is.finite(candidates) & candidates > 0]
  if (length(candidates)) return(max(candidates))
  parallel::detectCores(logical = TRUE)
}

if (!interactive()) {
  dir.create(timenow0, showWarnings = FALSE, recursive = TRUE)
  setwd(timenow0)
  sink(timenow, append = TRUE, split = TRUE)
  nc <- if (is.finite(cores_req) && cores_req > 0L) cores_req else get_allocated_cores()
  reserve <- if (is.finite(cores_req) && cores_req > 0L) 0L else 1L
  workers <- max(1L, min(SIMNUM, nc - reserve, 100L))
  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
  cl <- parallel::makeCluster(workers, type = "PSOCK", outfile = "cluster_worker.log")
  doParallel::registerDoParallel(cl)
  on.exit({ try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)
  message(sprintf("Batch backend: allocated=%d, workers=%d", nc, foreach::getDoParWorkers()))
} else {
  workers <- min(8L, parallel::detectCores() - 1L)
  cl <- parallel::makeCluster(workers, type = "PSOCK")
  doParallel::registerDoParallel(cl)
  on.exit({ try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)
}

## =========================================================
## Driver
## =========================================================
BIAS_res <- NULL; SE_res <- NULL; RMSE_res <- NULL; RB_res <- NULL; CR_res <- NULL
FDR_res <- NULL; FNR_res <- NULL
replicates <- list(); summaries <- list()

set.seed(2)
X0_dim <- if (vary_over == "p") max(seq) else as.integer(p_fixed)
X0 <- rmvnorm(n = N, mean = rep(0, X0_dim), sigma = ar1_cor(X0_dim, 0.2)) + 2
e <- rnorm(N, 0, 1)
len <- round(N / 4)

print(paste("SIMNUM =", SIMNUM, "| sampling =", sampling, "| fit =", if (wls_use) "WLS" else "OLS"))

for (val in seq) {
  if (vary_over == "p") { p <- val; r <- r_fixed } else { r <- val; p <- p_fixed }
  print(paste("p =", p, ", r =", r))

  z <- r * e + rnorm(n = N, 0, sqrt(1 - r^2))
  z_sorted_idx <- order(z)

  X <- X0[, 1:p, drop = FALSE]
  beta <- c(rep(1, s), rep(0, p - s))
  mu <- drop(X %*% beta)
  y <- drop(mu + e)
  t_y <- sum(y)
  X <- cbind(1, X)  # add intercept

  ## Design quantities that do not change across replicates
  pik <- NULL; Dmat <- NULL; wcan <- NULL; qmat <- NULL
  strat_id <- integer(N)
  for (h in 1:4) strat_id[z_sorted_idx[(len * (h - 1) + 1):(len * h)]] <- h
  if (sampling == "stratified") {
    piv <- (as.integer(c(15, 20, 30, 35) / 100 * n) / len)[strat_id]
  } else {
    pik  <- make_pi_from_z(z, n, gamma = 1)
    stopifnot(all(pik < 1))
    pi2  <- sampling::UPmaxentropypi2(pik)
    Dmat <- (pi2 - tcrossprod(pik)) / pi2
    rm(pi2)
    pikt <- sampling::UPMEpiktildefrompik(pik)
    wcan <- pikt / (1 - pikt)
    qmat <- sampling::UPMEqfromw(wcan, n)
    piv  <- pik
  }

  ## Target of the fitted working regression: intercept shift c_N of the census fit
  ## (pi-weighted for unweighted fits, unweighted for design-weighted fits).
  c_N <- if (wls_use) mean(e) else sum(piv * e) / sum(piv)
  mu_star <- mu + c_N
  print(paste("c_N =", signif(c_N, 4)))

  registerDoRNG(seed = 11)

  final_res <- foreach(
    simnum = 1:SIMNUM,
    .export = c("sample_design_stratified", "quad_Omega_strata", "sample_design_rejective",
                "var_ht_exact", "fit_ls", "fit_lasso_beta", "cf_predict", "make_pop_folds",
                "make_folds_strat_aware", "an_strat", "an_rej_cond", "an_rej", "run_estimators",
                "sampling", "N", "n", "s", "len", "z_sorted_idx", "y", "X", "mu", "mu_star",
                "K_use", "wls_use", "an_use", "lasso_only", "pik", "Dmat", "wcan", "qmat", "strat_id"),
    .packages = c("glmnet", "sampling"),
    .errorhandling = "pass"
  ) %dopar% {
    samp <- if (sampling == "stratified") {
      sample_design_stratified(N = N, n = n, len = len, z_sorted_idx = z_sorted_idx)
    } else {
      sample_design_rejective(pik, qmat)
    }
    out_w <- run_estimators(sampling = sampling, N = N, n = n, s = s, y = y, X = X,
                            mu = mu, mu_star = mu_star, samp = samp, K = K_use,
                            wls = wls_use, an = an_use, lasso_only = lasso_only,
                            pik = pik, Dmat = Dmat, wcan = wcan, strat_id = strat_id)
    list(res = out_w$y_res, res2 = out_w$sigma_res,
         res3 = as.numeric(abs(out_w$y_res - t_y) <= 1.96 * out_w$sigma_res),
         res4 = out_w$FDR, res5 = out_w$FNR)
  }

  ok <- sapply(final_res, function(x) is.list(x) && is.numeric(unlist(x$res)))
  print(paste("# of failure:", sum(!ok)))
  if (any(!ok)) print(final_res[!ok][[1]])

  res  <- do.call("rbind", lapply(final_res[ok], `[[`, "res"))
  res2 <- do.call("rbind", lapply(final_res[ok], `[[`, "res2"))
  res3 <- do.call("rbind", lapply(final_res[ok], `[[`, "res3"))
  colnames(res3) <- colnames(res)
  res4 <- unlist(lapply(final_res[ok], `[[`, "res4"))
  res5 <- unlist(lapply(final_res[ok], `[[`, "res5"))
  Bok <- nrow(res)

  BIAS <- colMeans(res - t_y)
  SD   <- apply(res, 2, function(x) sqrt(var(x) * (length(x) - 1) / length(x)))
  RMSE <- apply(res - t_y, 2, function(x) sqrt(mean(x^2)))
  RB   <- (colMeans(res2^2) - SD^2) / SD^2
  CR   <- colMeans(res3)
  CR_MCSE <- sqrt(CR * (1 - CR) / Bok)
  CIL  <- colMeans(2 * qnorm(0.975) * res2)

  summ <- data.frame(estimator = colnames(res), BIAS = BIAS, SD = SD, RMSE = RMSE,
                     RB = RB, CR = CR, CR_MCSE = CR_MCSE, CIL = CIL, B = Bok,
                     p = p, r = r, c_N = c_N, row.names = NULL)
  write.csv(summ, sprintf("summary_%s_%s.csv", vary_over, format(val)), row.names = FALSE)
  print(xtable(summ[, c("estimator", "BIAS", "SD", "RMSE", "RB", "CR", "CIL")], digits = 2,
               caption = sprintf("p = %s, r = %s", p, r)), include.rownames = FALSE)

  summaries[[as.character(val)]] <- summ
  replicates[[as.character(val)]] <- list(est = res, sigma = res2, t_y = t_y, c_N = c_N)

  SE_res   <- cbind(SE_res, SD)
  RMSE_res <- cbind(RMSE_res, RMSE)
  BIAS_res <- cbind(BIAS_res, BIAS)
  RB_res   <- cbind(RB_res, RB)
  CR_res   <- cbind(CR_res, CR)
  FDR_res  <- cbind(FDR_res, res4)
  FNR_res  <- cbind(FNR_res, res5)
}

if (!interactive()) {
  saveRDS(list(summaries = summaries, replicates = replicates,
               settings = list(SIMNUM = SIMNUM, N = N, n = n, s = s, K = K_use,
                               sampling = sampling, vary_over = vary_over, seq = seq,
                               wls = wls_use, an = an_use, lasso_only = lasso_only,
                               seed_population = 2, seed_replicates = 11)),
          "replicates.rds")
}

timenow2 <- Sys.time()
print("Running time")
print(timenow2 - timenow1)

colnames(FDR_res) <- seq
print(colMeans(FDR_res))
print(colMeans(FNR_res))

## =========================================================
## Plot (single panel for this run)
## =========================================================
colnames(RMSE_res) <- as.character(seq)
colnames(BIAS_res) <- as.character(seq)

df_plot <- rbind(
  data.frame(estimator = rep(rownames(RMSE_res), times = ncol(RMSE_res)),
             x = rep(as.numeric(colnames(RMSE_res)), each = nrow(RMSE_res)),
             value = as.vector(RMSE_res), metric = "RMSE"),
  data.frame(estimator = rep(rownames(BIAS_res), times = ncol(BIAS_res)),
             x = rep(as.numeric(colnames(BIAS_res)), each = nrow(BIAS_res)),
             value = as.vector(BIAS_res), metric = "Bias")
)

# Keep only the four estimators shown in the paper figure.
keep_est <- c("GREG", "DREG", "GREG.Lasso", "DREG.Lasso")
keep_est <- keep_est[keep_est %in% df_plot$estimator]
df_plot  <- subset(df_plot, estimator %in% keep_est)
df_plot$estimator <- factor(df_plot$estimator, levels = keep_est)
df_plot$x         <- factor(df_plot$x, levels = seq)
df_plot$metric    <- factor(df_plot$metric, levels = c("Bias", "RMSE"))

if (!interactive()) png("linegraph.png", width = 960, height = 560, res = 110)
print(
  ggplot(df_plot, aes(x = x, y = value, group = interaction(estimator, metric),
                      color = estimator, linetype = metric)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_line(linewidth = 1) +
    geom_point(data = subset(df_plot, metric == "RMSE"), size = 2) +
    scale_linetype_manual(values = c(RMSE = "solid", Bias = "dashed")) +
    labs(x = vary_over, y = "RMSE (solid) / Bias (dashed)", color = "Estimator", linetype = NULL) +
    theme_bw(base_size = 13) +
    theme(legend.position = "right", legend.key.width = unit(1.2, "lines"))
)
if (!interactive()) dev.off()

## Save the plotting data frame for combine_pr_plot.R.
if (!interactive()) saveRDS(df_plot, file = "plot_df.rds")
