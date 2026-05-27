# Usage (batch):
#   Rscript main.R --run --B <SIMNUM> [--N <pop>] [--n <samp>] \
#                  [--sampling stratified|rejective] [--vary_over p|r] \
#                  [--p <num|list>] [--r <num|list>] \
#                  [--wls] [--cores <int>] [--out <dir>]
#
# Notes:
#   --run       acts as a "do it" flag (any value accepted; presence is enough).
#   --B         number of Monte Carlo replicates (required in batch mode).
#   --N, --n    override population and sample size (defaults: 1000, 300).
#   --p, --r    override sweep grid / pinned value; accepts scalar or "a,b,c".
#               When --vary_over=p, --p is the sweep grid and --r is the pinned
#               value (must be scalar); vice versa for --vary_over=r.
#   --wls       use weighted least squares (weights = 1/pi) for the regression
#               fits. Bare flag means TRUE; you can also write --wls true/false.
#               Default is OLS.
#   --cores     number of parallel workers; overrides SLURM / detectCores().
#   --out       output directory name; defaults to a timestamped folder.

## -------- CLI parser (flag-style) --------
.parse_flags <- function(args) {
  out <- list()
  i <- 1L
  # Flags that may appear without an argument (presence == TRUE).
  bool_flags <- c("run", "wls")
  while (i <= length(args)) {
    a <- args[i]
    if (!startsWith(a, "--")) stop("Expected --flag, got: ", a)
    key <- sub("^--", "", a)
    if (key %in% bool_flags) {
      # Bare boolean: if the next token is missing or another flag, treat as TRUE.
      # Otherwise consume it as an explicit value (true/false/yes/no/1/0).
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
  # Accepts NULL, a single number, or a comma/space-separated list of numbers.
  if (is.null(s)) return(NULL)
  if (!is.character(s)) {
    stop(".parse_num_list expected character, got ", class(s)[1],
         " (value: ", paste(s, collapse = ", "), "). ",
         "This usually means an opts$<flag> partial-match against another key.")
  }
  toks <- strsplit(s, "[,[:space:]]+", perl = TRUE)[[1]]
  toks <- toks[nzchar(toks)]
  if (length(toks) == 0L) return(NULL)
  v <- suppressWarnings(as.numeric(toks))
  if (any(is.na(v))) {
    stop("Could not parse numeric list from '", s, "'.")
  }
  v
}

if (!interactive()) {
  opts <- .parse_flags(commandArgs(trailingOnly = TRUE))

  # NOTE: always use opts[["x"]], NEVER opts$x. R's $ does partial matching
  # on list element names; e.g. opts$r would silently match opts$run if --r
  # is not supplied, returning TRUE instead of NULL.

  if (is.null(opts[["B"]])) stop("Must supply --B <SIMNUM>.")
  SIMNUM <- as.numeric(opts[["B"]])
  if (!is.finite(SIMNUM) || SIMNUM <= 0) stop("--B must be a positive number.")

  sampling  <- if (!is.null(opts[["sampling"]]))  tolower(opts[["sampling"]])  else "stratified"
  vary_over <- if (!is.null(opts[["vary_over"]])) tolower(opts[["vary_over"]]) else "p"

  N_user    <- if (!is.null(opts[["N"]]))     as.integer(opts[["N"]])     else 1000L
  n_user    <- if (!is.null(opts[["n"]]))     as.integer(opts[["n"]])     else 300L
  cores_req <- if (!is.null(opts[["cores"]])) as.integer(opts[["cores"]]) else NA_integer_
  out_dir   <- opts[["out"]]          # NULL → fall back to timestamped name

  p_user <- .parse_num_list(opts[["p"]])   # NULL, scalar, or vector
  r_user <- .parse_num_list(opts[["r"]])

  # --wls accepts the bare flag (TRUE) or an explicit value (TRUE/FALSE/yes/no/1/0).
  wls_use <- if (is.null(opts[["wls"]])) FALSE else {
    v <- tolower(as.character(opts[["wls"]]))
    if (v %in% c("true", "t", "yes", "y", "1", "")) TRUE
    else if (v %in% c("false", "f", "no", "n", "0")) FALSE
    else stop("Could not parse --wls value '", opts[["wls"]],
              "'. Use true/false (or omit for OLS).")
  }

} else {
  SIMNUM    <- 100
  sampling  <- "stratified"      # or "rejective"
  vary_over <- "r"               # "p" or "r"
  N_user    <- 1000L
  n_user    <- 300L
  cores_req <- NA_integer_
  out_dir   <- NULL
  p_user    <- NULL
  r_user    <- NULL
  wls_use   <- FALSE
}

if (!sampling %in% c("stratified", "rejective")) {
  stop("Unknown sampling='", sampling, "'. Use 'stratified' or 'rejective'.")
}
if (!vary_over %in% c("p", "r")) {
  stop("Unknown vary_over='", vary_over, "'. Use 'p' or 'r'.")
}
if (sampling == "rejective") message("Rejective sampling")
message("Varying over: ", vary_over)
message("Regression fit: ", if (wls_use) "WLS (weighted by 1/pi)" else "OLS")

timenow1 <- Sys.time()
timenow0 <- gsub(" ", "_", gsub("[-:]", "", timenow1))
timenow0 <- paste0(c(sampling, vary_over, timenow0), collapse = "_")
if (!is.null(out_dir)) timenow0 <- out_dir   # honor --out
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
})

## -------- Parameters --------
N <- N_user
n <- n_user
s <- 5

# Defaults (used when the corresponding flag is not supplied).
default_seq_p <- c(10, 50, 90, 130, 170, 210, 250)   # p = 250 takes too long.
default_seq_r <- c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)
default_p     <- 90
default_r     <- -0.75

if (vary_over == "p") {
  # Sweep over p; r must be a scalar.
  seq_p   <- if (!is.null(p_user)) p_user else default_seq_p
  if (!is.null(r_user) && length(r_user) > 1L) {
    stop("--r must be a single value when --vary_over p (got ",
         length(r_user), " values: ", paste(r_user, collapse = ", "), ")")
  }
  r_fixed <- if (!is.null(r_user)) r_user else default_r
  p_fixed <- NA_real_                  # unused on this branch
  seq     <- seq_p
} else {
  # Sweep over r; p must be a scalar.
  if (!is.null(p_user) && length(p_user) > 1L) {
    stop("--p must be a single value when --vary_over r (got ",
         length(p_user), " values: ", paste(p_user, collapse = ", "), ")")
  }
  p_fixed <- if (!is.null(p_user)) as.integer(p_user) else default_p
  seq_r   <- if (!is.null(r_user)) r_user else default_seq_r
  r_fixed <- NA_real_                  # unused on this branch
  seq     <- seq_r
}

# Initialize p, r so they exist before the loop overwrites them.
p <- if (vary_over == "p") seq[1] else p_fixed
r <- if (vary_over == "r") seq[1] else r_fixed

message(sprintf("Sweep dimension: %s = {%s}", vary_over, paste(seq, collapse = ", ")))
if (vary_over == "p") message(sprintf("Fixed r = %s", r_fixed))
if (vary_over == "r") message(sprintf("Fixed p = %s", p_fixed))

seq_K <- round(c(10))  # K = N is jackknife

## -------- Utilities --------
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
  rho ^ exponent
}

## =========================================================
## Sampling designs
## =========================================================

## ---- Stratified (original main.R) ----
sample_design_stratified <- function(N, n, len, z_sorted_idx, y, X) {
  n_h <- as.integer(c(15, 20, 30, 35) / 100 * n)
  stopifnot(sum(n_h) == n)

  Index <- integer(n)
  d1 <- numeric(N)

  cumn_h <- cumsum(n_h)
  st <- c(1L, head(cumn_h, -1L) + 1L)

  for (i in 1:4) {
    Idx_z <- (len * (i - 1) + 1):(len * i)
    Idx <- z_sorted_idx[Idx_z]

    d1[Idx] <- len / n_h[i]
    Index[st[i]:cumn_h[i]] <- sample(Idx, size = n_h[i], replace = FALSE)
  }

  list(
    Index = Index,
    d1 = d1,
    d1_s = d1[Index],
    y_s = y[Index],
    X_s = X[Index, , drop = FALSE],
    # strata bookkeeping for exact HT-style quadratic form
    n_h = n_h,
    st = st,
    cumn_h = cumn_h,
    len = len
  )
}

quad_Omega_strata <- function(e, n_h, len, st, cumn_h) {
  out <- 0
  for (i in 1:4) {
    idx <- st[i]:cumn_h[i]
    eh <- e[idx]
    nh <- n_h[i]

    # constants from the Omega_tmp construction
    c_off <- len^2/nh^2 - len*(len-1)/nh/(nh-1)
    d_diag <- (1 - nh/len) * len^2/nh^2

    s1 <- sum(eh)
    s2 <- sum(eh*eh)
    out <- out + c_off * (s1*s1) + (d_diag - c_off) * s2
  }
  out
}

## ---- Rejective (original main_rejective.R) ----
make_pi_from_z <- function(z, n, gamma = 1) {
  mz <- as.numeric(scale(z))
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
    m0 <- m[idx0]
    pi[idx0] <- n_rem * m0 / sum(m0)
  }

  if (abs(sum(pi) - n) > 1e-8) {
    idx0 <- which(pi < 1)
    pi[idx0] <- pi[idx0] * (n - sum(pi[pi == 1])) / sum(pi[idx0])
  }
  pi <- pmax(pmin(pi, 1), 1e-12)
  pi
}

rejective_sample <- function(pi, n, max_tries = 1e6) {
  N <- length(pi)
  for (t in 1:max_tries) {
    s <- rbinom(N, 1, pi)
    if (sum(s) == n) return(s)
  }
  stop("Rejective sampler failed to hit fixed size; try adjusting pi or max_tries.")
}

library(sampling)

sample_design_rejective <- function(N, n, z, y, X, gamma = 1) {
  pik <- make_pi_from_z(z, n, gamma = gamma)   # sum(pik) == n, 0<pik<1
  s   <- sampling::UPmaxentropy(pik)           # 0/1 vector, fixed size
  
  Index <- which(s == 1)
  d1    <- 1 / pik
  
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
## Estimation / inference (shared core)
## =========================================================
run_estimators <- function(
  sampling,
  N, n, s, y, X, mu,
  # sample object fields (some may be NULL depending on sampling)
  Index, d1, d1_s, y_s, X_s, pi = NULL, n_h = NULL, st = NULL, cumn_h = NULL, len = NULL,
  seq_K, wls = FALSE
) {
  ## Base regression (GREG)
  if (wls) {
    lm_obj <- lm(y_s ~ 0 + X_s, weights = d1_s)
  } else {
    lm_obj <- lm(y_s ~ 0 + X_s)
  }
  beta_hat <- lm_obj$coefficients

  ## Oracle regression uses intercept + first s signals (X already has intercept col)
  X_s_oracle <- X[Index, 1:(1 + s), drop = FALSE]
  if (wls) {
    lm_obj_oracle <- lm(y_s ~ 0 + X_s_oracle, weights = d1_s)
  } else {
    lm_obj_oracle <- lm(y_s ~ 0 + X_s_oracle)
  }
  beta_hat_oracle <- lm_obj_oracle$coefficients

  ## Core totals
  y_HT   <- sum(y_s * d1_s)
  y_diff <- drop(sum(mu) + sum((y_s - mu[Index]) * d1_s))
  y_GREG <- drop(colSums(X) %*% beta_hat +
                   sum((y_s - drop(X_s %*% beta_hat)) * d1_s))
  # ## ---- Fast Tukey computation block (identical structure) ----
  # w <- if (wls) d1_s else rep(1, n)
  # 
  # sw <- sqrt(w)
  # Xw <- X_s * sw
  # yw <- y_s * sw
  # 
  # XtWX <- crossprod(Xw)
  # XtWy <- crossprod(Xw, yw)
  # 
  # R <- chol(XtWX)
  # beta_hat <- backsolve(R, forwardsolve(t(R), XtWy))
  # 
  # yhat <- drop(X_s %*% beta_hat)
  # e <- y_s - yhat
  # 
  # A <- chol2inv(R)
  # XA <- X_s %*% A
  # h <- w * rowSums(XA * X_s)
  # 
  # rk <- w * e / (1 - h)
  # 
  # X_tot  <- colSums(X)
  # Wy_sum <- sum(w * y_s)
  # Xw_sum <- colSums(X_s * w)
  # 
  # AtXtot  <- drop(A %*% X_tot)
  # v <- drop(X_s %*% AtXtot)
  # 
  # AtXw_sum <- drop(A %*% Xw_sum)
  # u <- drop(X_s %*% AtXw_sum)
  # 
  # term1 <- drop(crossprod(X_tot, beta_hat)) - rk * v
  # 
  # xbeta_mk_at_k <- yhat - h * (e / (1 - h))
  # sum_wx_beta_mk <- drop(crossprod(Xw_sum, beta_hat)) - rk * u
  # sum_wx_beta_excl_k <- sum_wx_beta_mk - w * xbeta_mk_at_k
  # ressum_excl_k <- (Wy_sum - w * y_s) - sum_wx_beta_excl_k
  # 
  # term2 <- ressum_excl_k * n / (n - 1)

  # y_GREG <- drop(crossprod(X_tot, beta_hat) + sum((y_s - yhat) * w))
  # y_Tukey <- n * y_GREG - (n - 1) / n * sum(term1 + term2)

  ## Oracle GREG
  y_GREG_oracle <- drop(
    colSums(X[,1:(1 + s), drop = FALSE]) %*% beta_hat_oracle +
      sum((y_s - drop(X_s_oracle %*% beta_hat_oracle)) * d1_s)
  )

  ## Lasso (CV per simulation)
  if (wls) {
    cv_model <- cv.glmnet(X_s, y_s, weights = d1_s)
    best_model <- glmnet(X_s, y_s, weights = d1_s, lambda = cv_model$lambda.min)
  } else {
    cv_model <- cv.glmnet(X_s, y_s)
    best_model <- glmnet(X_s, y_s, lambda = cv_model$lambda.min)
  }
  beta_hat_Lasso <- as.vector(coef(best_model))[-1]  # length == ncol(X)

  y_Lasso <- drop(
    colSums(X) %*% beta_hat_Lasso +
      sum((y_s - drop(X_s %*% beta_hat_Lasso)) * d1_s)
  )

  ## -------- Variances (design-specific) --------
  if (sampling == "stratified") {
    sigma_HT   <- sqrt(quad_Omega_strata(y_s, n_h, len, st, cumn_h))

    e_s <- y_s - mu[Index]
    sigma_diff <- sqrt(quad_Omega_strata(e_s, n_h, len, st, cumn_h))

    e_s <- y_s - X_s %*% beta_hat
    sigma_GREG <- sqrt(quad_Omega_strata(e_s, n_h, len, st, cumn_h))

    # # Tukey JK by strata (your stratified version)
    # T_i <- term1 + term2
    # V_jk <- 0
    # for (hh in seq(length(st))) {
    #   idx <- st[hh]:cumn_h[hh]
    #   Th_bar <- mean(T_i[idx])
    #   V_jk <- V_jk + (n_h[hh] - 1) / n_h[hh] * sum((T_i[idx] - Th_bar)^2)
    # }
    # sigma_Tukey <- sqrt(V_jk)

    e_s <- y_s - X_s_oracle %*% beta_hat_oracle
    sigma_GREG_oracle <- sqrt(quad_Omega_strata(e_s, n_h, len, st, cumn_h))

    e_s <- y_s - X_s %*% beta_hat_Lasso
    sigma_Lasso <- sqrt(quad_Omega_strata(e_s, n_h, len, st, cumn_h))

  } else { # rejective
    sigma_HT   <- sqrt(var_poisson_approx(y_s,               pi[Index]))
    sigma_diff <- sqrt(var_poisson_approx(y_s - mu[Index],   pi[Index]))
    sigma_GREG <- sqrt(var_poisson_approx(y_s - drop(X_s %*% beta_hat),       pi[Index]))
    sigma_Lasso<- sqrt(var_poisson_approx(y_s - drop(X_s %*% beta_hat_Lasso), pi[Index]))

    # sigma_Tukey <- sqrt((n - 1) / n * sum((term1 + term2 - mean(term1 + term2))^2))

    # keep same label set (oracle available; variance uses Poisson approx)
    sigma_GREG_oracle <- sqrt(var_poisson_approx(y_s - drop(X_s_oracle %*% beta_hat_oracle), pi[Index]))
  }

  ## -------- Assemble base outputs --------
  y_res <- c(
    HT = y_HT,
    Diff = y_diff,
    # Tukey = y_Tukey,
    GREG.oracle = y_GREG_oracle,    
    GREG = y_GREG,
    GREG.Lasso = y_Lasso
  )

  sigma_res <- c(
    HT = sigma_HT,
    Diff = sigma_diff,
    # Tukey = sigma_Tukey,
    GREG.oracle = sigma_GREG_oracle,    
    GREG = sigma_GREG,
    GREG.Lasso = sigma_Lasso
  )

  tmpIndex <- numeric(length(y))
  tmpIndex[Index] <- d1_s

  bias_res <- c(
    HT = 0,
    Diff = 0,
    GREG.oracle = cov(tmpIndex, drop(y - X[,1:(1 + s), drop=FALSE] %*% beta_hat_oracle)) * N,
    GREG = cov(tmpIndex, drop(y - X %*% beta_hat)) * N,
    GREG.Lasso = cov(tmpIndex, drop(y - X %*% beta_hat_Lasso)) * N
  )

  ## FDR/FNR (same definition as originals)
  FDR <- sum(beta_hat_Lasso[-c(2:(1 + s))] != 0) / (ncol(X) - s)
  FNR <- sum(beta_hat_Lasso[c(2:(1 + s))] == 0) / s

  ## -------- Sample-splitting debiasing (shared; design-specific sigma) --------
  for (K in seq_K) {
    y_debiased_vec <- NULL
    y_debiased_Lasso_vec <- NULL
    e_s_vec <- list()
    e_s_vec_Lasso <- list()
    bias_vec <- matrix(0, nrow = K, ncol = 2)

    shuffled <- sample(1:N)
    groups <- cut(seq_along(shuffled), breaks = K, labels = FALSE)
    SubIndex_list <- split(shuffled, groups)

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

      if (wls) {
        lm_obj2 <- lm(y_s2 ~ 0 + X_s2, weights = d1_s2)
      } else {
        lm_obj2 <- lm(y_s2 ~ 0 + X_s2)
      }
      beta_hat2 <- lm_obj2$coefficients

      y_debiased <- drop(colSums(X[SubIndex, , drop = FALSE]) %*% beta_hat2 +
                           sum((y_s1 - drop(X_s1 %*% beta_hat2)) * d1_s1))

      e_tmp <- drop(y_s1 - X_s1 %*% beta_hat2)
      names(e_tmp) <- Index_sub1
      e_s_vec <- append(e_s_vec, list(e_tmp))

      bias_vec[k, 1] <- cov(tmpIndex[SubIndex],
                            drop(y[SubIndex] - X[SubIndex, , drop=FALSE] %*% beta_hat2)) * length(SubIndex)

      y_debiased_vec <- c(y_debiased_vec, y_debiased)

      if (wls) {
        cv2 <- cv.glmnet(X_s2, y_s2, weights = d1_s2)
        fit2 <- glmnet(X_s2, y_s2, weights = d1_s2, lambda = cv2$lambda.min)
      } else {
        cv2 <- cv.glmnet(X_s2, y_s2)
        fit2 <- glmnet(X_s2, y_s2, lambda = cv2$lambda.min)
      }
      beta_hat2_L <- as.vector(coef(fit2))[-1]

      y_debiased_Lasso <- drop(
        colSums(X[SubIndex, , drop = FALSE]) %*% beta_hat2_L +
          sum((y_s1 - drop(X_s1 %*% beta_hat2_L)) * d1_s1)
      )

      e_tmp_L <- drop(y_s1 - X_s1 %*% beta_hat2_L)
      names(e_tmp_L) <- Index_sub1
      e_s_vec_Lasso <- append(e_s_vec_Lasso, list(e_tmp_L))

      bias_vec[k, 2] <- cov(tmpIndex[SubIndex],
                            drop(y[SubIndex] - X[SubIndex, , drop = FALSE] %*% beta_hat2_L)) * length(SubIndex)

      y_debiased_Lasso_vec <- c(y_debiased_Lasso_vec, y_debiased_Lasso)
    }

    y_debiased <- sum(y_debiased_vec)
    y_debiased_Lasso <- sum(y_debiased_Lasso_vec)

    e_s_ss  <- unlist(e_s_vec)[as.character(Index)]
    e_s_ssl <- unlist(e_s_vec_Lasso)[as.character(Index)]

    if (sampling == "stratified") {
      sigma_SSGREG  <- sqrt(quad_Omega_strata(e_s_ss,  n_h, len, st, cumn_h))
      sigma_SSLasso <- sqrt(quad_Omega_strata(e_s_ssl, n_h, len, st, cumn_h))
    } else {
      sigma_SSGREG  <- sqrt(var_poisson_approx(e_s_ss,  pi[Index]))
      sigma_SSLasso <- sqrt(var_poisson_approx(e_s_ssl, pi[Index]))
    }

    y_res_tmp <- setNames(c(y_debiased, y_debiased_Lasso), c("DREG", "DREG.Lasso"))
    y_res <- c(y_res[1:(length(y_res) - 1)], y_res_tmp[1], y_res[length(y_res)], y_res_tmp[2])

    sigma_res_tmp <- setNames(c(sigma_SSGREG, sigma_SSLasso), c("DREG", "DREG.Lasso"))
    sigma_res <- c(sigma_res[1:(length(sigma_res) - 1)], sigma_res_tmp[1], sigma_res[length(sigma_res)], sigma_res_tmp[2])

    bias_res_tmp <- setNames(colSums(bias_vec), c("DREG", "DREG.Lasso"))
    bias_res <- c(bias_res[1:(length(bias_res) - 1)], bias_res_tmp[1], bias_res[length(bias_res)], bias_res_tmp[2])
  }

  list(
    y_res = y_res,
    sigma_res = sigma_res,
    bias_res = bias_res,
    FDR = FDR,
    FNR = FNR
  )
}

## =========================================================
## Parallel setup (shared)
## =========================================================
get_allocated_cores <- function() {
  slurm <- Sys.getenv("SLURM_CPUS_PER_TASK")
  pbs   <- Sys.getenv("PBS_NP")
  lsf   <- Sys.getenv("LSB_DJOB_NUMPROC")

  candidates <- suppressWarnings(as.integer(c(slurm, pbs, lsf)))
  candidates <- candidates[is.finite(candidates) & candidates > 0]

  if (length(candidates)) return(max(candidates))
  parallel::detectCores(logical = TRUE)
}

if (!interactive()) {
  dir.create(timenow0, showWarnings = FALSE, recursive = TRUE)
  setwd(timenow0)
  sink(timenow, append = TRUE, split = TRUE)   # split=TRUE: keep echoing to stdout

  if (is.finite(cores_req) && cores_req > 0L) {
    nc <- cores_req                            # honor --cores verbatim
  } else {
    nc <- get_allocated_cores()
  }
  reserve <- if (is.finite(cores_req) && cores_req > 0L) 0L else 1L

  workers_cap <- 100L
  workers <- max(1L, min(SIMNUM, nc - reserve, workers_cap))

  Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")

  suppressMessages(library(doParallel))
  cl <- parallel::makeCluster(workers, type = "PSOCK", outfile = "cluster_worker.log")

  doParallel::registerDoParallel(cl)
  on.exit({ try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)

  message(sprintf("Batch backend: requested=%s, allocated=%d, workers=%d",
                  ifelse(is.finite(cores_req), as.character(cores_req), "auto"),
                  nc, foreach::getDoParWorkers()))

  .ok <- tryCatch({
    aa <- foreach::foreach(i = 1:foreach::getDoParWorkers(), .combine = c) %dopar% i
    length(aa) == foreach::getDoParWorkers()
  }, error = function(e) FALSE)
  if (!.ok) stop("Parallel backend smoke test failed")

} else {
  suppressMessages(library(doParallel))
  print(paste("detectCores =", parallel::detectCores()))
  workers <- min(8L, parallel::detectCores() - 1L)
  cl <- parallel::makeCluster(workers, type = "PSOCK")
  doParallel::registerDoParallel(cl)
  on.exit({ try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)
}

## =========================================================
## Driver (shared)
## =========================================================
BIAS_res <- NULL
SE_res <- NULL
RMSE_res <- NULL
RB_res <- NULL
CR_res <- NULL
FDR_res <- NULL
FNR_res <- NULL
BIASapprox_res <- NULL

set.seed(2)
# X0 must be wide enough to support the largest p actually used.
X0_dim <- if (vary_over == "p") max(seq) else as.integer(p_fixed)
X0 <- rmvnorm(n = N, mean = rep(0, X0_dim),
              sigma = ar1_cor(X0_dim, 0.2)) + 2
e <- rnorm(N, 0, 1)
len <- round(N / 4)

print(paste("SIMNUM =", SIMNUM))
print(paste("sampling =", sampling))

for (val in seq) {
  if (vary_over == "p") {
    p <- val
    r <- r_fixed
  } else {
    r <- val
    p <- p_fixed
  }
  
  print(paste("p =", p, ", r =", r))

  z <- r * e + rnorm(n = N, 0, sqrt(1 - r^2))
  z_sorted_idx <- order(z)

  X <- X0[, 1:p, drop = FALSE]

  beta <- c(rep(1, s), rep(0, p - s))
  mu <- X %*% beta
  y <- drop(mu + e)
  t_y <- sum(y)

  X <- cbind(1, X)  # add intercept

  registerDoRNG(seed = 11)

  final_res <- foreach(
    simnum = 1:SIMNUM,
    .export = c(
      # sampling
      "sample_design_stratified", "quad_Omega_strata",
      "sample_design_rejective", "make_pi_from_z", "rejective_sample", "var_poisson_approx",
      # estimation
      "run_estimators",
      # params
      "sampling", "N", "n", "s", "p", "len", "z_sorted_idx", "z", "y", "X", "mu", "seq_K", "wls_use"
    ),
    .packages = c("glmnet"),
    .errorhandling = "pass"
  ) %dopar% {

    if (sampling == "stratified") {
      samp <- sample_design_stratified(N = N, n = n, len = len, z_sorted_idx = z_sorted_idx, y = y, X = X)

      out_w <- run_estimators(
        sampling = sampling,
        N=N, n=n, s=s, y=y, X=X, mu=mu,
        Index=samp$Index, d1=samp$d1, d1_s=samp$d1_s, y_s=samp$y_s, X_s=samp$X_s,
        n_h=samp$n_h, st=samp$st, cumn_h=samp$cumn_h, len=samp$len,
        seq_K=seq_K, wls = wls_use
      )

    } else {
      samp <- sample_design_rejective(N = N, n = n, z = z, y = y, X = X, gamma = 1)

      out_w <- run_estimators(
        sampling = sampling,
        N=N, n=n, s=s, y=y, X=X, mu=mu,
        Index=samp$Index, pi=samp$pi, d1=samp$d1, d1_s=samp$d1_s, y_s=samp$y_s, X_s=samp$X_s,
        seq_K=seq_K, wls = wls_use
      )
    }

    list(
      res  = c(out_w$y_res),
      res2 = c(out_w$sigma_res),
      res3 = c(ifelse(abs(out_w$y_res - t_y) > 1.96 * out_w$sigma_res, 0, 1)),
      res4 = out_w$FDR,
      res5 = out_w$FNR,
      res6 = c(out_w$bias_res)
    )
  }

  final_res1 <- lapply(final_res, `[[`, "res")
  final_res2 <- lapply(final_res, `[[`, "res2")
  final_res3 <- lapply(final_res, `[[`, "res3")
  final_res4 <- lapply(final_res, `[[`, "res4")
  final_res5 <- lapply(final_res, `[[`, "res5")
  final_res6 <- lapply(final_res, `[[`, "res6")

  ok <- sapply(final_res1, function(x) is.numeric(unlist(x)))
  print(paste("# of failure:", sum(!ok)))

  res  <- do.call("rbind", final_res1[ok])
  res2 <- do.call("rbind", final_res2[ok])
  res3 <- do.call("rbind", final_res3[ok])
  res4 <- do.call("rbind", final_res4[ok])
  res5 <- do.call("rbind", final_res5[ok])
  res6 <- do.call("rbind", final_res6[ok])

  BIAS <- colMeans(res - t_y)
  SE <- apply(res, 2, function(x) sqrt(var(x) * (length(x) - 1) / length(x)))
  RMSE <- apply(res - t_y, 2, function(x) sqrt(mean(x^2)))

  tmpdf21 <- cbind(BIAS, SE, RMSE)
  print(xtable(tmpdf21, digits = 2, caption = "Summary of point estimation"))

  BIAS2 <- colMeans(res2^2) - SE^2
  tmpdf22 <- cbind(RB = BIAS2 / SE^2, CR = colMeans(res3))
  print(xtable(tmpdf22, digits = 2, caption = "Summary of variance estimation"))

  boxplot(res, main = paste("p =", p, ", r =", r, ", sampling =", sampling), cex.axis = 0.5)
  abline(h = t_y, lty = 1, col = 2)
  abline(v = c(2.5, 6.5), lty = 3)

  SE_res <- cbind(SE_res, SE)
  RMSE_res <- cbind(RMSE_res, RMSE)
  BIAS_res <- cbind(BIAS_res, BIAS)
  RB_res <- cbind(RB_res, BIAS2 / SE^2)
  CR_res <- cbind(CR_res, colMeans(res3))
  FDR_res <- cbind(FDR_res, res4)
  FNR_res <- cbind(FNR_res, res5)
  BIASapprox_res <- cbind(BIASapprox_res, colMeans(res6))
}

if (!interactive()) save.image(paste0(basename(timenow0), ".RData"))

timenow2 <- Sys.time()
print("Running time")
print(timenow2 - timenow1)

print(xtable(cbind(
  cbind(BIAS = BIAS_res[,1], SE = SE_res[,1], RMSE = RMSE_res[,1]),
  cbind(BIAS = BIAS_res[,ncol(BIAS_res)], SE = SE_res[,ncol(SE_res)], RMSE = RMSE_res[,ncol(RMSE_res)])
)))

print(xtable(cbind(
  RB = RB_res[,1], CR = CR_res[,1],
  RB = RB_res[,ncol(RB_res)], CR = CR_res[,ncol(CR_res)]
)))

colnames(FDR_res) <- seq
colMeans(FDR_res)
colMeans(FNR_res)

## =========================================================
## Plot (single panel for this run)
## =========================================================
colnames(RMSE_res) <- as.character(seq)
colnames(BIAS_res) <- as.character(seq)

df_rmse <- data.frame(
  estimator = rep(rownames(RMSE_res), times = ncol(RMSE_res)),
  x         = rep(as.numeric(colnames(RMSE_res)), each = nrow(RMSE_res)),
  value     = as.vector(RMSE_res),
  metric    = "RMSE"
)
df_bias <- data.frame(
  estimator = rep(rownames(BIAS_res), times = ncol(BIAS_res)),
  x         = rep(as.numeric(colnames(BIAS_res)), each = nrow(BIAS_res)),
  value     = as.vector(BIAS_res),
  metric    = "Bias"
)

df_plot <- rbind(df_rmse, df_bias)

# Keep only the four estimators shown in the paper figure.
keep_est <- c("GREG", "DREG", "GREG.Lasso", "DREG.Lasso")
df_plot  <- subset(df_plot, estimator %in% keep_est)
df_plot$estimator <- factor(df_plot$estimator, levels = keep_est)
df_plot$x         <- factor(df_plot$x, levels = seq)
df_plot$metric    <- factor(df_plot$metric, levels = c("Bias", "RMSE"))

x_label <- vary_over   # "p" or "r"

if (!interactive()) png("linegraph.png", width = 960, height = 560, res = 110)

print(
  ggplot(
    df_plot,
    aes(x = x, y = value,
        group = interaction(estimator, metric),
        color = estimator,
        linetype = metric)
  ) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_line(linewidth = 1) +
    geom_point(data = subset(df_plot, metric == "RMSE"), size = 2) +
    scale_linetype_manual(values = c(RMSE = "solid", Bias = "dashed")) +
    labs(
      x        = x_label,
      y        = "RMSE (solid) / Bias (dashed)",
      color    = "Estimator",
      linetype = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "right",
          legend.key.width = unit(1.2, "lines"))
)

if (!interactive()) dev.off()

## Also save the plotting data frame so the combined-panel script can pick it up.
if (!interactive()) saveRDS(df_plot, file = "plot_df.rds")

############################################################
## End unified single-file script                           ##
############################################################
