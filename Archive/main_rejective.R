############################################################
## Single-file version (main + sample + run), no sourcing  ##
############################################################

print("Rejective sampling")

if (!interactive()) {
  args <- as.numeric(commandArgs(trailingOnly = TRUE))
} else{
  args <- c(15)
}

timenow1 <- Sys.time()
timenow0 = gsub(' ', '_', gsub('[-:]', '', timenow1))
timenow = paste(timenow0, ".txt", sep = "")

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

## Parameters
N <- 500
p <- 90
n <- 300
s <- 5
SIMNUM <- args[1]
r <- 0.75        # Noninformative (your original default)

# seq_p <- c(10, 20)
# seq_p <- c(10, 30, 50, 70, 90, 110, 130)
seq_p <- c(10, 50, 90, 130, 170, 210, 250) # p = 250 takes too long time.
# seq_r <- c(0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75)
seq_r <- c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)
seq_K <- round(c(N))  # K = N is jackknife

## -------- Utilities --------
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
  rho ^ exponent
}

## -------- Rejective sampling design --------
make_pi_from_z <- function(z, n, gamma = 1) {
  # Positive size measure from z
  mz <- as.numeric(scale(z))
  m  <- exp(gamma * mz)
  
  # Initial scaling
  pi <- n * m / sum(m)
  
  # Enforce 0 < pi <= 1 and sum(pi) = n via "peeling"
  # (Iteratively fix pi=1 units and rescale remaining.)
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
  
  # Guardrails
  if (abs(sum(pi) - n) > 1e-8) {
    # final tiny renormalization on unfixed units (rare)
    idx0 <- which(pi < 1)
    pi[idx0] <- pi[idx0] * (n - sum(pi[pi == 1])) / sum(pi[idx0])
  }
  pi <- pmax(pmin(pi, 1), 1e-12)
  pi
}

rejective_sample <- function(pi, n, max_tries = 1e6) {
  # Simple rejective sampler: independent Bernoulli(pi) until sum==n
  # Works well when n is not extreme and pi are not too concentrated.
  N <- length(pi)
  for (t in 1:max_tries) {
    s <- rbinom(N, 1, pi)
    if (sum(s) == n) return(s)
  }
  stop("Rejective sampler failed to hit fixed size; try adjusting pi or max_tries.")
}

sample_design_rejective <- function(N, n, z, y, X, gamma = 1) {
  pi <- make_pi_from_z(z, n, gamma = gamma)
  s  <- rejective_sample(pi, n)
  
  Index <- which(s == 1)
  d1 <- 1 / pi
  
  list(
    Index = Index,
    pi = pi,
    d1 = d1,
    d1_s = d1[Index],
    y_s = y[Index],
    X_s = X[Index, , drop = FALSE]
  )
}



## -------- Approx variance for rejective (Poisson-style) --------
var_poisson_approx <- function(e_s, pi_s) {
  sum((1 - pi_s) / (pi_s^2) * (e_s^2))
}



## -------- Estimation/inference (was run.R) --------
run_estimators <- function(N, n, s, y, X, mu, Index, pi, d1, d1_s, y_s, X_s, seq_K, wls = TRUE) {
  ## Base regression (GREG)
  if (wls) {
    lm_obj <- lm(y_s ~ 0 + X_s, weights = d1_s)
  } else{
    lm_obj <- lm(y_s ~ 0 + X_s)
  }
  beta_hat <- lm_obj$coefficients
  
  ## Oracle regression uses intercept + first s signals (your X already has intercept col)
  X_s_oracle <- X[Index, 1:(1 + s), drop = FALSE]
  if (wls) {
    lm_obj_oracle <- lm(y_s ~ 0 + X_s_oracle, weights = d1_s)
  } else{
    lm_obj_oracle <- lm(y_s ~ 0 + X_s_oracle)
  }
  beta_hat_oracle <- lm_obj_oracle$coefficients
  
  ## Core totals
  y_HT   <- sum(y_s * d1_s)
  y_diff <- drop(sum(mu) + sum((y_s - mu[Index]) * d1_s))
  # y_GREG <- drop(colSums(X) %*% beta_hat +
  #                  sum((y_s - drop(X_s %*% beta_hat)) * d1_s))
  # 
  # ## Tukey / jackknife variants (kept identical)
  # y_Tukey <- length(y_s) * y_GREG - (length(y_s) - 1) / length(y_s) * sum(
  #   sapply(1:length(y_s), function(k){
  #     if (wls) {
  #       betahat_k <- drop(solve(
  #         t(X_s[-k,,drop=FALSE]) %*% diag(d1_s[-k]) %*% X_s[-k,,drop=FALSE],
  #         t(X_s[-k,,drop=FALSE]) %*% diag(d1_s[-k]) %*% y_s[-k]
  #       ))
  #     } else{
  #       betahat_k <- drop(solve(
  #         t(X_s[-k,,drop=FALSE]) %*% X_s[-k,,drop=FALSE],
  #         t(X_s[-k,,drop=FALSE]) %*% y_s[-k]
  #       ))
  #     }
  #     drop(colSums(X) %*% betahat_k) +
  #       sum((y_s[-k] - X_s[-k,,drop=FALSE] %*% betahat_k) * d1_s[-k]) * n / (n - 1)
  #   })
  # )
  
  w <- if (wls) d1_s else rep(1, n)
  
  # ----- single WLS / OLS fit -----
  sw <- sqrt(w)
  Xw <- X_s * sw
  yw <- y_s * sw
  
  XtWX <- crossprod(Xw)
  XtWy <- crossprod(Xw, yw)
  
  R <- chol(XtWX)
  beta_hat <- backsolve(R, forwardsolve(t(R), XtWy))
  
  yhat <- drop(X_s %*% beta_hat)
  e <- y_s - yhat
  
  # ----- hat diagonals -----
  A <- chol2inv(R)
  XA <- X_s %*% A
  h <- w * rowSums(XA * X_s)
  
  rk <- w * e / (1 - h)
  
  # ----- precomputations -----
  X_tot  <- colSums(X)
  Wy_sum <- sum(w * y_s)
  Xw_sum <- colSums(X_s * w)
  
  AtXtot  <- drop(A %*% X_tot)
  v <- drop(X_s %*% AtXtot)
  
  AtXw_sum <- drop(A %*% Xw_sum)
  u <- drop(X_s %*% AtXw_sum)
  
  # ----- term 1 -----
  term1 <- drop(crossprod(X_tot, beta_hat)) - rk * v
  
  # ----- term 2 -----
  xbeta_mk_at_k <- yhat - h * (e / (1 - h))
  sum_wx_beta_mk <- drop(crossprod(Xw_sum, beta_hat)) - rk * u
  sum_wx_beta_excl_k <- sum_wx_beta_mk - w * xbeta_mk_at_k
  ressum_excl_k <- (Wy_sum - w * y_s) - sum_wx_beta_excl_k
  
  term2 <- ressum_excl_k * n / (n - 1)
  
  # ----- GREG and Tukey -----
  y_GREG <- drop(crossprod(X_tot, beta_hat) + sum((y_s - yhat) * w))
  y_Tukey <- n * y_GREG - (n - 1) / n * sum(term1 + term2)
  
  
  # y_NaiveJ <- mean(sapply(1:length(y_s), function(k){
  #   if (wls) {
  #     betahat_k <- drop(solve(
  #       t(X_s[-k,,drop=FALSE]) %*% diag(d1_s[-k]) %*% X_s[-k,,drop=FALSE],
  #       t(X_s[-k,,drop=FALSE]) %*% diag(d1_s[-k]) %*% y_s[-k]
  #     ))
  #   } else{
  #     betahat_k <- drop(solve(
  #       t(X_s[-k,,drop=FALSE]) %*% X_s[-k,,drop=FALSE],
  #       t(X_s[-k,,drop=FALSE]) %*% y_s[-k]
  #     ))
  #   }
  #   drop(colSums(X) %*% betahat_k) + sum((y_s - X_s %*% betahat_k) * d1_s)
  # }))
  # 
  # y_SampleJ <- mean(sapply(1:length(y_s), function(k){
  #   if (wls) {
  #     betahat_k <- drop(solve(
  #       t(X_s[-k,,drop=FALSE]) %*% diag(d1_s[-k]) %*% X_s[-k,,drop=FALSE],
  #       t(X_s[-k,,drop=FALSE]) %*% diag(d1_s[-k]) %*% y_s[-k]
  #     ))
  #   } else{
  #     betahat_k <- drop(solve(
  #       t(X_s[-k,,drop=FALSE]) %*% X_s[-k,,drop=FALSE],
  #       t(X_s[-k,,drop=FALSE]) %*% y_s[-k]
  #     ))
  #   }
  #   drop(colSums(X) %*% betahat_k + (y_s[k] - X_s[k,,drop=FALSE] %*% betahat_k) * d1_s[k])
  # }))
  
  y_GREG_oracle <- drop(
    colSums(X[,1:(1 + s), drop = FALSE]) %*% beta_hat_oracle +
      sum((y_s - drop(X_s_oracle %*% beta_hat_oracle)) * d1_s)
  )
  
  ## Lasso (do CV per simulation; no simnum==0 “cache” needed)
  if (wls) {
    cv_model <- cv.glmnet(X_s, y_s, weights = d1_s)
    best_model <- glmnet(X_s, y_s, weights = d1_s, lambda = cv_model$lambda.min)
  } else{
    cv_model <- cv.glmnet(X_s, y_s)
    best_model <- glmnet(X_s, y_s, lambda = cv_model$lambda.min)
  }
  beta_hat_Lasso <- as.vector(coef(best_model))[-1]  # length == ncol(X)
  
  y_Lasso <- drop(
    colSums(X) %*% beta_hat_Lasso +
      sum((y_s - drop(X_s %*% beta_hat_Lasso)) * d1_s)
  )
  
  sigma_HT   <- sqrt(var_poisson_approx(y_s,       pi[Index]))
  sigma_diff <- sqrt(var_poisson_approx(y_s - mu[Index], pi[Index]))
  sigma_GREG <- sqrt(var_poisson_approx(y_s - drop(X_s %*% beta_hat), pi[Index]))
  sigma_Lasso<- sqrt(var_poisson_approx(y_s - drop(X_s %*% beta_hat_Lasso), pi[Index]))
  
  
  sigma_Tukey <- sqrt((n - 1) / n * sum((term1 + term2 - mean(term1 + term2))^2))
  
  ## Assemble base outputs
  y_res <- c(
    HT = y_HT,
    Diff = y_diff,
    Tukey = y_Tukey,    
    GREG = y_GREG,
    GREG_oracle = y_GREG_oracle,
    # NaiveJ = y_NaiveJ,
    # SampleJ = y_SampleJ,
    GREG.Lasso = y_Lasso
  )
  
  sigma_res <- c(
    HT = sigma_HT,
    Diff = sigma_diff,
    Tukey = sigma_Tukey,    
    GREG = sigma_GREG,
    GREG_oracle = sigma_GREG_oracle,
    # NaiveJ = NA,
    # SampleJ = NA,
    GREG.Lasso = sigma_Lasso
  )
  
  tmpIndex <- numeric(length(y))
  tmpIndex[Index] <- d1_s
  
  bias_res <- c(
    HT = 0,
    Diff = 0,
    GREG = cov(tmpIndex, drop(y - X %*% beta_hat)) * N,
    GREG_oracle = cov(tmpIndex, drop(y - X[,1:(1 + s), drop=FALSE] %*% beta_hat_oracle)) * N,
    GREG.Lasso = cov(tmpIndex, drop(y - X %*% beta_hat_Lasso)) * N
  )
  
  ## FDR/FNR (same definitions as your original driver uses)
  # Note: beta_hat_Lasso corresponds to all columns of X (including intercept)
  # but here beta_hat_Lasso excludes intercept because coef(best_model)[-1].
  # In your original, X had intercept as first column and beta_hat_Lasso matched that.
  # To keep identical behavior, define "p+1" as ncol(X) and treat first coefficient as intercept.
  # We'll reconstruct an intercept-included vector for FDR/FNR parity.
  # FDR/FNR consistent with your original:
  # true signals are columns 2:(1+s) in X (since col 1 is the intercept column you added)
  FDR <- sum(beta_hat_Lasso[-c(2:(1 + s))] != 0) / (ncol(X) - s)
  FNR <- sum(beta_hat_Lasso[c(2:(1 + s))] == 0) / s
  
  
  ## Sample-splitting debiasing over seq_K (kept structurally same)
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
      
      ## OLS/WLS on complement
      if (wls) {
        lm_obj2 <- lm(y_s2 ~ 0 + X_s2, weights = d1_s2)
      } else{
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
      
      ## Lasso on complement (do CV on complement each fold; matches your “cv_model2” behavior)
      if (wls) {
        cv2 <- cv.glmnet(X_s2, y_s2, weights = d1_s2)
        fit2 <- glmnet(X_s2, y_s2, weights = d1_s2, lambda = cv2$lambda.min)
      } else{
        cv2 <- cv.glmnet(X_s2, y_s2)
        fit2 <- glmnet(X_s2, y_s2, lambda = cv2$lambda.min)
      }
      beta_hat2_L <- as.vector(coef(fit2))[-1]  # length == ncol(X)
      
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
    
    e_s_ss <- unlist(e_s_vec)[as.character(Index)]
    sigma_SSGREG   <- sqrt(var_poisson_approx(e_s_ss,       pi[Index]))
    
    e_s_ssl <- unlist(e_s_vec_Lasso)[as.character(Index)]
    sigma_SSLasso   <- sqrt(var_poisson_approx(e_s_ssl,       pi[Index]))    
    
    y_res_tmp <- setNames(c(y_debiased, y_debiased_Lasso),
                          c("SREG", "SREG.Lasso"))
    y_res = c(y_res[1:(length(y_res) - 1)], y_res_tmp[1], y_res[length(y_res)], y_res_tmp[2])
    
    # y_res <- c(y_res, setNames(c(y_debiased, y_debiased_Lasso),
    #                            c("SREG", "SREG.Lasso")))
    
    sigma_res_tmp <- setNames(c(sigma_SSGREG, sigma_SSLasso),
                              c("SREG", "SREG.Lasso"))
    sigma_res = c(sigma_res[1:(length(sigma_res) - 1)], sigma_res_tmp[1], sigma_res[length(sigma_res)], sigma_res_tmp[2])
    
    # sigma_res <- c(sigma_res, setNames(c(sigma_SSGREG, sigma_SSLasso),
    #                                    c("SREG", "SREG.Lasso")))
    
    bias_res_tmp <- setNames(colSums(bias_vec),
                             c("SREG", "SREG.Lasso"))
    bias_res = c(bias_res[1:(length(bias_res) - 1)], bias_res_tmp[1], bias_res[length(bias_res)], bias_res_tmp[2])
    
    # bias_res <- c(bias_res, setNames(colSums(bias_vec),
    #                                  c("SREG", "SREG.Lasso")))
  }
  
  list(
    y_res = y_res,
    sigma_res = sigma_res,
    bias_res = bias_res,
    beta_hat2_L = beta_hat2_L,
    FDR = FDR,
    FNR = FNR
  )
}

## -------- Driver (was main.R) --------

## Parallel setup
# cores <- min(parallel::detectCores() - 1, 101)
# print(paste("cores =", cores))
# cl <- makeCluster(cores)
# registerDoParallel(cl)

## ---------- Parallel setup (FAST + robust) ----------
get_allocated_cores <- function() {
  # Common schedulers
  slurm <- Sys.getenv("SLURM_CPUS_PER_TASK")
  pbs   <- Sys.getenv("PBS_NP")
  lsf   <- Sys.getenv("LSB_DJOB_NUMPROC")
  
  candidates <- suppressWarnings(as.integer(c(slurm, pbs, lsf)))
  candidates <- candidates[is.finite(candidates) & candidates > 0]
  
  if (length(candidates)) return(max(candidates))
  parallel::detectCores(logical = TRUE)
}

if (!interactive()) {
  dir.create(timenow0, showWarnings = FALSE)
  setwd(timenow0)
  sink(timenow, append = TRUE)
  
  nc <- get_allocated_cores()
  reserve <- 1L
  
  # IMPORTANT: cap workers to avoid OOM kills when MC is large
  workers_cap <- 100L                 # try 4L if your server RAM is tight
  workers <- max(1L, min(SIMNUM, nc - reserve, workers_cap))
  
  # prevent BLAS oversubscription
  Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
  
  suppressMessages(library(doParallel))
  
  # IMPORTANT: prefer PSOCK for robustness (avoid makeForkCluster)
  cl <- parallel::makeCluster(workers, type = "PSOCK", outfile = "cluster_worker.log")
  
  doParallel::registerDoParallel(cl)
  on.exit({ try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)
  
  message(sprintf("Batch backend: allocated=%d, workers=%d", nc, foreach::getDoParWorkers()))
  
  .ok <- tryCatch({
    aa <- foreach::foreach(i = 1:foreach::getDoParWorkers(), .combine = c) %dopar% i
    length(aa) == foreach::getDoParWorkers()
  }, error = function(e) FALSE)
  if (!.ok) stop("Parallel backend smoke test failed")
  
} else {
  suppressMessages(library(doParallel))
  print(paste("detectCores =", parallel::detectCores()))
  
  # also cap here for consistency
  workers <- min(8L, parallel::detectCores() - 1L)
  cl <- parallel::makeCluster(workers, type = "PSOCK")
  
  doParallel::registerDoParallel(cl)
  on.exit({ try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)
}


## Results containers
BIAS_res <- NULL
SE_res <- NULL
RMSE_res <- NULL
RB_res <- NULL
CR_res <- NULL
FDR_res <- NULL
FNR_res <- NULL
BIASapprox_res <- NULL

## Data generation (fixed across p, as in your original)
set.seed(2)
X0 <- rmvnorm(n = N, mean = rep(0, seq_p[length(seq_p)]),
              sigma = ar1_cor(seq_p[length(seq_p)], 0.2)) + 2
e <- rnorm(N, 0, 1)

len <- round(N / 4)

print(paste("SIMNUM =", SIMNUM))

# seq <- seq_p # Change this 1
# for (p in seq_p) {
  
  seq <- seq_r # Change this 2
  for (r in seq_r) {
  print(paste("p =", p, ", r =", r))
  
  z = r * e + rnorm(n = N, 0, 1 - r^2)
  # if (r == 0) {
  #   z <- rnorm(n = N, 0, 1)
  # } else{
  #   z <- rnorm(n = N, 0, sqrt((1 - r) / r)) + e
  # }
  z_sorted_idx <- order(z)
  
  X <- X0[, 1:p, drop = FALSE]
  
  beta <- c(rep(1, s), rep(0, p - s))
  # mu <- cbind(X[, 1:s, drop = FALSE]^(1:s %% 3 + 1) / 20, # nonlinear
  #             X[, (s + 1):ncol(X), drop = FALSE]) %*% beta
  mu <- X %*% beta
  y <- drop(mu + e)
  t_y <- sum(y)
  
  ## add intercept (matches your original)
  X <- cbind(1, X)
  
  registerDoRNG(seed = 11)
  
  final_res <- foreach(
    simnum = 1:SIMNUM,
    .export = c("sample_design_rejective", "make_pi_from_z", "rejective_sample",
                "run_estimators", "N", "n", "s", "p", "y", "X", "mu", "seq_K", "z"),
    .packages = c("glmnet"),
    .errorhandling = "pass"
  ) %dopar% {
    
    ## sampling
    samp <- sample_design_rejective(N = N, n = n, z = z, y = y, X = X, gamma = 1)
    
    ## run twice (unweighted then weighted) like your original logic
    # out_unw <- run_estimators(
    #   N=N, n=n, s=s, y=y, X=X, mu=mu,
    #   Index=samp$Index, d1=samp$d1, d1_s=samp$d1_s, y_s=samp$y_s, X_s=samp$X_s, n_h=samp$n_h, st = samp$st, cumn_h = samp$cumn_h, len = samp$len,
    #   seq_K=seq_K, wls=FALSE
    # )
    
    out_w <- run_estimators(
      N=N, n=n, s=s, y=y, X=X, mu=mu,
      Index=samp$Index, pi=samp$pi, d1=samp$d1, d1_s=samp$d1_s,
      y_s=samp$y_s, X_s=samp$X_s,
      seq_K=seq_K, wls=TRUE
    )
    
    
    # list(
    #   res = c(out_unw$y_res[1:2], out_w$y_res),
    #   res2 = c(out_unw$sigma_res[1:2], out_w$sigma_res),
    #   res3 = c(
    #     ifelse(abs(out_unw$y_res[1:2] - t_y) > 1.96 * out_unw$sigma_res[1:2], 0, 1),
    #     ifelse(abs(out_w$y_res - t_y) > 1.96 * out_w$sigma_res, 0, 1)
    #   ),
    #   res4 = out_w$FDR,
    #   res5 = out_w$FNR,
    #   res6 = c(out_unw$bias_res[1:2], out_w$bias_res)
    # )
    
    list(
      res = c(out_w$y_res),
      res2 = c(out_w$sigma_res),
      res3 = c(
        ifelse(abs(out_w$y_res - t_y) > 1.96 * out_w$sigma_res, 0, 1)
      ),
      res4 = out_w$FDR,
      res5 = out_w$FNR,
      res6 = c(out_w$bias_res)
    )
  }
  
  ## Aggregate
  final_res1 <- lapply(final_res, `[[`, "res")
  final_res2 <- lapply(final_res, `[[`, "res2")
  final_res3 <- lapply(final_res, `[[`, "res3")
  final_res4 <- lapply(final_res, `[[`, "res4")
  final_res5 <- lapply(final_res, `[[`, "res5")
  final_res6 <- lapply(final_res, `[[`, "res6")
  
  ## Drop failures
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
  # print(cbind(tmpdf21, tmpdf22))
  # print(xtable(cbind(tmpdf21, tmpdf22), digits = c(2)))
  
  # print("FD summary")
  # print(summary(res4 * ( (p + 1) - s )))
  # print("FN summary")
  # print(summary(res5 * s))
  
  boxplot(res, main = paste("p =", p), cex.axis = 0.5)
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
if(!interactive()) save.image(paste(timenow0, ".RData", sep = ""))

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

############################################################
## End single-file script                                  ##
############################################################

library(ggplot2)

# Ensure columns correspond to seq
colnames(RMSE_res) <- as.character(seq)
colnames(BIAS_res) <- as.character(seq)

# RMSE long
df_rmse <- data.frame(
  estimator = rep(rownames(RMSE_res), times = ncol(RMSE_res)),
  p = rep(as.numeric(colnames(RMSE_res)), each = nrow(RMSE_res)),
  value = as.vector(RMSE_res),
  metric = "RMSE"
)

# Bias long (absolute bias is usually what you want on RMSE scale)
df_bias <- data.frame(
  estimator = rep(rownames(BIAS_res), times = ncol(BIAS_res)),
  p = rep(as.numeric(colnames(BIAS_res)), each = nrow(BIAS_res)),
  value = (as.vector(BIAS_res)),
  metric = "Bias"
)

# Combine
df_plot <- rbind(df_rmse, df_bias)

# Drop Diff estimator
df_plot <- subset(df_plot, estimator != "Diff")
df_plot <- subset(df_plot, estimator != "HT")

# Fix estimator order
df_plot$estimator <- factor(df_plot$estimator, levels = rownames(RMSE_res))

# Fix p order
df_plot$p <- factor(df_plot$p, levels = seq)



if (!interactive()) png("linegraph.png", width = 960, height = 560)

ggplot(
  df_plot,
  aes(
    x = p,
    y = value,
    group = interaction(estimator, metric),
    color = estimator,
    linetype = metric
  )
) +
  geom_line(linewidth = 1) +
  geom_point(data = subset(df_plot, metric == "RMSE"), size = 2) +
  scale_linetype_manual(
    values = c(RMSE = "solid", Bias = "dashed")
  ) +
  labs(
    x = NULL,
    y = "RMSE (solid) / Bias (dashed)",
    color = "Estimator",
    linetype = NULL
  ) +
  theme_bw()

if (!interactive()) dev.off()


