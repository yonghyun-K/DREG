# ============================================================================
# Kynetec calibration workflow (organized)
#
# Goals
#   1) Build/merge crop-level survey data with weighting-file metadata
#   2) Construct calibration matrix X and universe totals for aux variables
#   3) Produce multiple weight sets (WeightGen, calibration/GREG, proposed CVXR)
#   4) Compare point estimators for SALES totals:
#        - HT
#        - REG (OLS GREG)
#        - SREG (sample-split OLS GREG; K-fold cross-fit)
#        - REG.Spline / SREG.Spline (GAM spline over CRD, by crop)
#        - REG.RF / SREG.RF (random forest via ranger on Crop x CRD)
#
# Notes
#   - This script keeps your original logic, but reorganizes it into sections
#     and small helper functions.
#   - It does NOT rely on setwd(); instead set `project_dir` below.
# ============================================================================

## ---- USER CONFIG -----------------------------------------------------------
# Path to the folder that contains ftns.R and UniverseFiles_revised/
project_dir <- "G:/My Drive/Kynetec_Internship/RShiny"  # <-- EDIT

# State abbreviation used in StateCodes.xlsx (e.g., "IA")
stateNm <- "SC"  # <-- EDIT CO, SC

# Calibration controls
min_cell_n <- 100                # min sample count for Farms-by-Size cells; smaller cells are merged
bound_l <- 0.15                  # lower bound multiplier (for diagnostics)
bound_u <- 4                     # upper bound multiplier (for diagnostics)

# Sample-splitting controls (apply once)
ss_K <- 5                         # default: 2-fold cross-fit

ss_seed <- 20250120               # set seed so your split is reproducible
ss_wls <- FALSE                   # match main.R default (unweighted regression)

# Nonparametric model controls
gam_k <- 10                      # basis dimension upper bound for spline over CRD (mgcv)
rf_num_trees <- 200              # number of trees for ranger (increase for final; decrease for speed)
rf_mtry <- NULL                  # if NULL, ranger default
rf_min_node_size <- 20           # larger => shallower trees => faster
rf_num_threads <- max(1, parallel::detectCores() - 1)  # use multiple CPU cores
rf_verbose <- FALSE              # suppress 'Growing trees...' progress output

# Output
save_rdata <- FALSE
save_csv   <- FALSE

# Targets of interest for estimation + plotting
#   - By default, focus on pesticide-type totals only for these 3 categories.
#   - Set to NULL to evaluate all pesticide types.
#   - Set include_company_targets <- TRUE to add company totals too.
# targets_of_interest <- c("Herbicide", "Fungicide")
targets_of_interest <- c("Herbicide")
include_company_targets <- FALSE

## ---- PACKAGES --------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  library(tidyr)
})


# plotly is optional (only used if you uncomment interactive plots below)

## ---- HELPERS ---------------------------------------------------------------
stop_if_missing <- function(path, label = NULL) {
  if (!file.exists(path)) {
    msg <- if (!is.null(label)) sprintf("Missing %s: %s", label, path) else sprintf("Missing file: %s", path)
    stop(msg)
  }
}

substr_right <- function(x, n) substr(x, nchar(x) - n + 1, nchar(x))

collapse_dup_cols <- function(df, pattern, new_name) {
  cols <- grep(pattern, names(df), value = TRUE)
  if (!length(cols)) return(df)

  df[[new_name]] <- rowSums(df[, cols, drop = FALSE], na.rm = TRUE)

  # remove all duplicated/suffixed versions but keep the newly created column
  rm_cols <- setdiff(cols, new_name)
  df[, rm_cols] <- NULL
  df
}

make_folds <- function(n, K = 2, seed = 1) {
  stopifnot(K >= 2, n >= K)
  set.seed(seed)
  fold_id <- sample(rep(seq_len(K), length.out = n))
  fold_id
}

safe_lm_coef <- function(y, X, w = NULL) {
  X <- as.matrix(X)
  y <- as.numeric(y)

  if (is.null(w)) {
    fit <- stats::lm.fit(x = X, y = y)
    beta <- fit$coefficients
  } else {
    fit <- stats::lm.wfit(x = X, y = y, w = as.numeric(w))
    beta <- fit$coefficients
  }

  beta[is.na(beta)] <- 0
  beta
}

estimate_total_ht <- function(y, d) {
  sum(d * y)
}

estimate_total_hajek <- function(y, d, N_total) {
  # Hajek (ratio) total = N * (HT mean)
  #   - If N_total is unknown, you can pass N_total = sum(d)
  (sum(d * y) / sum(d)) * N_total
}

greg_fit <- function(y, X_s, d_s, wls = FALSE) {
  # Fits a (possibly weighted) linear model with NO implicit intercept.
  # If you want an intercept, include a column of 1s in X_s.
  beta_hat <- safe_lm_coef(y = y, X = X_s, w = if (wls) d_s else NULL)
  yhat <- drop(as.matrix(X_s) %*% beta_hat)
  list(beta = beta_hat, resid = as.numeric(y) - yhat)
}

estimate_total_greg <- function(y, X_s, d_s, X_total, wls = FALSE, return_resid = FALSE) {
  # Regression (GREG) total using known auxiliary totals X_total:
  #   t_hat = X_total' beta_hat + sum_s d_i (y_i - x_i' beta_hat)
  fit <- greg_fit(y = y, X_s = X_s, d_s = d_s, wls = wls)

  total_hat <- sum(X_total * fit$beta) + sum(d_s * fit$resid)

  if (return_resid) {
    return(list(total = total_hat, beta = fit$beta, resid = fit$resid))
  }
  total_hat
}

estimate_total_sreg <- function(y, X_s, d_s, X_total, fold_id, wls = FALSE, return_resid = FALSE) {
  # Cross-fitted (sample-split) regression total.
  #
  # This adapts your main.R sample-splitting idea to the case where we only know
  # overall auxiliary totals X_total (not unit-level X for the whole population).
  #
  # Returns:
  #   - total estimate
  #   - (optionally) cross-fitted residuals for variance estimation

  X_s <- as.matrix(X_s)
  y <- as.numeric(y)
  d_s <- as.numeric(d_s)

  K <- max(fold_id)
  p <- ncol(X_s)

  # HT estimated X totals overall and by fold
  Xhat_total <- drop(crossprod(d_s, X_s))   # length p

  Xhat_by_fold <- vector("list", K)
  for (k in seq_len(K)) {
    idx_k <- which(fold_id == k)
    Xhat_by_fold[[k]] <- drop(crossprod(d_s[idx_k], X_s[idx_k, , drop = FALSE]))
  }

  # Allocate true totals to folds proportional to HT fold totals, variable-by-variable
  # (component-wise scaling). This guarantees sum_k X_total_k == X_total.
  scale_factor <- rep(0, p)
  nz <- abs(Xhat_total) > 1e-12
  scale_factor[nz] <- X_total[nz] / Xhat_total[nz]

  X_total_by_fold <- lapply(Xhat_by_fold, function(v) v * scale_factor)

  # Fold contributions + cross-fitted residuals (main.R-style)
  total_hat <- 0
  resid_cf <- rep(NA_real_, length(y))

  for (k in seq_len(K)) {
    test  <- which(fold_id == k)
    train <- which(fold_id != k)

    beta_k <- safe_lm_coef(y = y[train], X = X_s[train, , drop = FALSE], w = if (wls) d_s[train] else NULL)

    pred_k <- sum(X_total_by_fold[[k]] * beta_k)
    res_k  <- sum(d_s[test] * (y[test] - drop(X_s[test, , drop = FALSE] %*% beta_k)))

    total_hat <- total_hat + pred_k + res_k

    # store cross-fitted residuals for the held-out observations
    resid_cf[test] <- y[test] - drop(X_s[test, , drop = FALSE] %*% beta_k)
  }

  if (return_resid) {
    return(list(total = total_hat, resid = resid_cf))
  }

  total_hat
}




# --- Nonparametric model-assisted estimators: spline (GAM) + random forest (ranger) ----

estimate_total_gam_greg <- function(y, df_s, d_s, pop_df,
                                   wls = FALSE, gam_k = 10, seed = 1,
                                   return_resid = FALSE) {
  # Model-assisted total: sum_U mhat(x) + sum_s d_i (y_i - mhat(x_i))
  # mhat is fitted using mgcv::gam (spline over CRD), with crop-specific smooth.

  if (!requireNamespace("mgcv", quietly = TRUE)) {
    warning("Package 'mgcv' is not installed. REG.Spline/SREG.Spline will be NA. Install via install.packages('mgcv').")
    resid <- rep(NA_real_, length(y))
    if (return_resid) return(list(total = NA_real_, resid = resid))
    return(NA_real_)
  }

  set.seed(seed)

  y <- as.numeric(y)
  d_s <- as.numeric(d_s)

  df <- df_s
  df$y <- y

  # Choose a safe basis dimension based on available CRD support
  n_crd <- length(unique(df$CRD_num[is.finite(df$CRD_num)]))
  k_use <- max(4, min(as.integer(gam_k), as.integer(n_crd)))

  w <- if (wls) d_s else NULL

  # Try GAM; fall back to two-way factor model if GAM fails
  fit <- tryCatch({
    mgcv::gam(
      y ~ CropNm + mgcv::s(CRD_num, CropNm, bs = "fs", k = k_use),
      data = df,
      weights = w,
      method = "REML"
    )
  }, error = function(e) {
    NULL
  })

  if (is.null(fit)) {
    # fallback: saturated by (Crop, CRD) via factors
    fit <- stats::lm(y ~ CropNm + CRD_fac, data = df, weights = w)
    pred_s <- stats::predict(fit, newdata = df)
    pred_U <- stats::predict(fit, newdata = pop_df)
  } else {
    pred_s <- stats::predict(fit, newdata = df, type = "response")
    pred_U <- stats::predict(fit, newdata = pop_df, type = "response")
  }

  pred_s[!is.finite(pred_s)] <- 0
  pred_U[!is.finite(pred_U)] <- 0

  resid <- y - as.numeric(pred_s)
  total_pred_U <- sum(as.numeric(pop_df$Farms) * as.numeric(pred_U))
  total_hat <- total_pred_U + sum(d_s * resid)

  if (return_resid) return(list(total = total_hat, resid = resid))
  total_hat
}

estimate_total_gam_sreg <- function(y, df_s, d_s, pop_df, fold_id,
                                   wls = FALSE, gam_k = 10, seed = 1,
                                   return_resid = FALSE) {
  # Cross-fitted version of GAM model-assisted estimator.
  # Uses out-of-fold residuals for main.R-style variance estimation.

  if (!requireNamespace("mgcv", quietly = TRUE)) {
    warning("Package 'mgcv' is not installed. REG.Spline/SREG.Spline will be NA. Install via install.packages('mgcv').")
    resid_cf <- rep(NA_real_, length(y))
    if (return_resid) return(list(total = NA_real_, resid = resid_cf))
    return(NA_real_)
  }

  y <- as.numeric(y)
  d_s <- as.numeric(d_s)
  df_full <- df_s
  df_full$y <- y

  K <- max(fold_id)
  pred_totals <- numeric(K)
  resid_cf <- rep(NA_real_, length(y))

  for (k in seq_len(K)) {
    test  <- which(fold_id == k)
    train <- which(fold_id != k)

    df_tr <- df_full[train, , drop = FALSE]
    df_te <- df_full[test,  , drop = FALSE]

    w_tr <- if (wls) d_s[train] else NULL

    n_crd <- length(unique(df_tr$CRD_num[is.finite(df_tr$CRD_num)]))
    k_use <- max(4, min(as.integer(gam_k), as.integer(n_crd)))

    set.seed(seed + 1000L * k)

    fit <- tryCatch({
      mgcv::gam(
        y ~ CropNm + mgcv::s(CRD_num, CropNm, bs = "fs", k = k_use),
        data = df_tr,
        weights = w_tr,
        method = "REML"
      )
    }, error = function(e) {
      NULL
    })

    if (is.null(fit)) {
      fit <- stats::lm(y ~ CropNm + CRD_fac, data = df_tr, weights = w_tr)
      pred_U <- stats::predict(fit, newdata = pop_df)
      pred_te <- stats::predict(fit, newdata = df_te)
    } else {
      pred_U <- stats::predict(fit, newdata = pop_df, type = "response")
      pred_te <- stats::predict(fit, newdata = df_te, type = "response")
    }

    pred_U[!is.finite(pred_U)] <- 0
    pred_te[!is.finite(pred_te)] <- 0

    pred_totals[k] <- sum(as.numeric(pop_df$Farms) * as.numeric(pred_U))

    resid_cf[test] <- y[test] - as.numeric(pred_te)
  }

  total_hat <- mean(pred_totals) + sum(d_s * resid_cf)

  if (return_resid) return(list(total = total_hat, resid = resid_cf))
  total_hat
}

estimate_total_rf_greg <- function(y, df_s, d_s, pop_df,
                                  wls = FALSE,
                                  num_trees = 200, mtry = NULL, min_node_size = 20,
                                  num_threads = 1, verbose = FALSE,
                                  seed = 1,
                                  return_resid = FALSE) {
  # Random-forest model-assisted total using ranger.

  if (!requireNamespace("ranger", quietly = TRUE)) {
    warning("Package 'ranger' is not installed. REG.RF/SREG.RF will be NA. Install via install.packages('ranger').")
    resid <- rep(NA_real_, length(y))
    if (return_resid) return(list(total = NA_real_, resid = resid))
    return(NA_real_)
  }

  set.seed(seed)

  y <- as.numeric(y)
  d_s <- as.numeric(d_s)

  df <- df_s
  df$y <- y

  # ranger uses case.weights for weighted fitting
  cw <- if (wls) d_s else NULL

  fit <- ranger::ranger(
    formula = y ~ CropNm + CRD_fac,
    data = df,
    num.trees = as.integer(num_trees),
    mtry = mtry,
    min.node.size = as.integer(min_node_size),
    case.weights = cw,
    respect.unordered.factors = "partition",
    num.threads = as.integer(num_threads),
    verbose = isTRUE(verbose),
    seed = as.integer(seed)
  )

  pred_s <- predict(fit, data = df)$predictions
  pred_U <- predict(fit, data = pop_df)$predictions

  pred_s[!is.finite(pred_s)] <- 0
  pred_U[!is.finite(pred_U)] <- 0

  resid <- y - as.numeric(pred_s)
  total_pred_U <- sum(as.numeric(pop_df$Farms) * as.numeric(pred_U))
  total_hat <- total_pred_U + sum(d_s * resid)

  if (return_resid) return(list(total = total_hat, resid = resid))
  total_hat
}

estimate_total_rf_sreg <- function(y, df_s, d_s, pop_df, fold_id,
                                  wls = FALSE,
                                  num_trees = 200, mtry = NULL, min_node_size = 20,
                                  num_threads = 1, verbose = FALSE,
                                  seed = 1,
                                  return_resid = FALSE) {
  # Cross-fitted ranger estimator + out-of-fold residuals.

  if (!requireNamespace("ranger", quietly = TRUE)) {
    warning("Package 'ranger' is not installed. REG.RF/SREG.RF will be NA. Install via install.packages('ranger').")
    resid_cf <- rep(NA_real_, length(y))
    if (return_resid) return(list(total = NA_real_, resid = resid_cf))
    return(NA_real_)
  }

  y <- as.numeric(y)
  d_s <- as.numeric(d_s)

  df_full <- df_s
  df_full$y <- y

  K <- max(fold_id)
  pred_totals <- numeric(K)
  resid_cf <- rep(NA_real_, length(y))

  for (k in seq_len(K)) {
    test  <- which(fold_id == k)
    train <- which(fold_id != k)

    df_tr <- df_full[train, , drop = FALSE]
    df_te <- df_full[test,  , drop = FALSE]

    cw <- if (wls) d_s[train] else NULL

    fit <- ranger::ranger(
      formula = y ~ CropNm + CRD_fac,
      data = df_tr,
      num.trees = as.integer(num_trees),
      mtry = mtry,
      min.node.size = as.integer(min_node_size),
      case.weights = cw,
      respect.unordered.factors = "partition",
      num.threads = as.integer(num_threads),
      verbose = isTRUE(verbose),
      seed = as.integer(seed + 1000L * k)
    )

    pred_U <- predict(fit, data = pop_df)$predictions
    pred_te <- predict(fit, data = df_te)$predictions

    pred_U[!is.finite(pred_U)] <- 0
    pred_te[!is.finite(pred_te)] <- 0

    pred_totals[k] <- sum(as.numeric(pop_df$Farms) * as.numeric(pred_U))
    resid_cf[test] <- y[test] - as.numeric(pred_te)
  }

  total_hat <- mean(pred_totals) + sum(d_s * resid_cf)

  if (return_resid) return(list(total = total_hat, resid = resid_cf))
  total_hat
}

preprocess_farmsize_for_variance <- function(data, stateNm) {
  # Matches ad-hoc harmonization rules used in Sampling_states.R
  data <- data %>%
    mutate(FarmSize = ifelse(Crop %in% c(87, 92), 1, FarmSize)) %>%
    mutate(FarmSize = ifelse(stateNm == "WI" & Crop == 91, 3, FarmSize)) %>%
    mutate(FarmSize = ifelse(stateNm == "OK" & Crop == 31, 1, FarmSize)) %>%
    mutate(FarmSize = ifelse(stateNm == "NY" & Crop == 2 & FarmSize == 7, 6, FarmSize))
  data
}

build_design_from_farms_by_size <- function(data, FarmsBySize, CropCodes, stateNm) {
  # Build stratified SRSWOR design info using (Crop, FarmSize) strata.
  # N_h comes from FarmsBySize.xlsx; n_h from sample counts.

  # 1) Universe stratum sizes
  farms_long <- FarmsBySize %>%
    filter(State == stateNm) %>%
    select(-2, -3, -ncol(FarmsBySize))

  farms_long <- reshape2::melt(farms_long)
  farms_long <- farms_long[!is.na(farms_long$value), ]
  names(farms_long) <- c("CropDesc", "FarmSize", "Farms")

  farms_long$CropDesc <- trimws(as.character(farms_long$CropDesc))
  crop_map <- CropCodes %>% select(CropCode, CropDescription)
  farms_long <- left_join(farms_long, crop_map, by = c("CropDesc" = "CropDescription"))

  if (any(is.na(farms_long$CropCode))) {
    missing <- unique(farms_long$CropDesc[is.na(farms_long$CropCode)])
    stop(
      "Some crop names in FarmsBySize.xlsx could not be mapped to CropCodes.xlsx. Examples: ",
      paste(head(missing, 10), collapse = ", ")
    )
  }

  farms_long <- farms_long %>%
    transmute(
      Crop = as.integer(CropCode),
      FarmSize = as.double(FarmSize),
      Farms = as.numeric(Farms)
    )

  # 2) Sample counts by stratum
  data2 <- data %>% mutate(Crop = as.integer(Crop), FarmSize = as.double(FarmSize))
  samp_n <- data2 %>% count(Crop, FarmSize, name = "n")

  # 3) Join N_h and n_h to each sampled unit
  strata_tbl <- left_join(samp_n, farms_long, by = c("Crop", "FarmSize"))

  if (any(is.na(strata_tbl$Farms))) {
    bad <- strata_tbl %>% filter(is.na(Farms))
    stop(
      "Design strata not found in FarmsBySize.xlsx for some sampled (Crop, FarmSize) cells. ",
      "First few missing rows:\n",
      paste(utils::capture.output(print(utils::head(bad, 10))), collapse = "\n")
    )
  }

  data_design <- left_join(data2, strata_tbl %>% select(Crop, FarmSize, Farms, n), by = c("Crop", "FarmSize"))

  N_h <- as.numeric(data_design$Farms)
  n_h <- as.integer(data_design$n)

  if (any(!is.finite(N_h)) || any(!is.finite(n_h))) stop("Non-finite N_h or n_h in design build")
  if (any(n_h <= 0)) stop("Found nonpositive sample size (n_h) in design build")

  # Guard: if a universe cell is <= sample count, bump by 1 (as in Sampling_states.R)
  N_h_adj <- ifelse(N_h <= n_h, n_h + 1, N_h)

  d <- N_h_adj / n_h

  stratum_id <- interaction(data_design$Crop, data_design$FarmSize, drop = TRUE)

  strata_params <- data.frame(
    stratum = as.character(stratum_id),
    N_h = N_h_adj,
    n_h = n_h
  ) %>% distinct()

  # Main.R-compatible stratum quadratic-form constants
  strata_params <- strata_params %>%
    mutate(
      c_off = ifelse(n_h > 1,
                    N_h^2 / n_h^2 - N_h * (N_h - 1) / n_h / (n_h - 1),
                    0),
      d_diag = ifelse(n_h > 1,
                     (1 - n_h / N_h) * N_h^2 / n_h^2,
                     0)
    )

  N_total <- sum(strata_params$N_h)

  list(
    d = d,
    stratum_id = stratum_id,
    strata_params = strata_params,
    N_total = N_total
  )
}

varhat_stratified <- function(e, stratum_id, strata_params) {
  # Computes sum_h [ c_off_h * (sum e_h)^2 + (d_diag_h - c_off_h) * sum(e_h^2) ]
  # which matches main.R's quad_Omega_strata, but for arbitrary strata sizes.

  e <- as.numeric(e)
  if (length(e) != length(stratum_id)) stop("e and stratum_id lengths do not match")

  strata_key <- as.character(stratum_id)

  sum1 <- tapply(e, strata_key, sum)
  sum2 <- tapply(e * e, strata_key, sum)

  # align to strata_params
  s1 <- sum1[strata_params$stratum]
  s2 <- sum2[strata_params$stratum]

  s1[is.na(s1)] <- 0
  s2[is.na(s2)] <- 0

  with(strata_params, sum(c_off * (s1 * s1) + (d_diag - c_off) * s2))
}


collapse_sparse_farmsize_cells <- function(Xs, FarmsTotal, total, cropNames, min_n = 100) {
  # Merge farm-size dummy columns with small sample counts (< min_n)
  # by moving them into an adjacent size class within the same crop.

  Xs <- as.data.frame(Xs)

  while (TRUE) {
    moved <- 0

    for (CropNm in cropNames) {
      namesvec <- names(FarmsTotal)[grepl(CropNm, names(FarmsTotal))]
      if (length(namesvec) <= 1) next

      for (j in seq_along(namesvec)) {
        col_j <- namesvec[j]

        # choose neighbor (previous if last, else next)
        if (j == length(namesvec)) {
          col_to <- namesvec[j - 1]
        } else {
          col_to <- namesvec[j + 1]
        }

        # sample count in that cell
        # (Columns can be removed earlier in the same pass; guard against that.)
        if (!col_j %in% names(Xs)) next
        if (!col_j %in% names(FarmsTotal)) next
        if (!col_j %in% names(total)) next

        if (!col_to %in% names(Xs)) next
        if (!col_to %in% names(FarmsTotal)) next
        if (!col_to %in% names(total)) next

        if (sum(Xs[[col_j]], na.rm = TRUE) < min_n) {
          moved <- moved + 1

          # move 1s into neighbor and drop small column
          Xs[[col_to]][Xs[[col_j]] == 1] <- 1
          Xs[[col_to]][is.na(Xs[[col_to]])] <- 0

          Xs[[col_j]] <- NULL

          # update totals (universe + combined)
          # NOTE: FarmsTotal/total are named numeric vectors, not lists.
          # Using `[[name]] <- NULL` errors for vectors, so we remove by subsetting.
          if (is.na(FarmsTotal[col_to])) FarmsTotal[col_to] <- 0
          FarmsTotal[col_to] <- FarmsTotal[col_to] + FarmsTotal[col_j]
          FarmsTotal <- FarmsTotal[names(FarmsTotal) != col_j]

          if (is.na(total[col_to])) total[col_to] <- 0
          total[col_to] <- total[col_to] + total[col_j]
          total <- total[names(total) != col_j]
        }
      }
    }

    message("# merged farm-size cells this pass: ", moved)
    if (moved == 0) break
  }

  list(Xs = Xs, FarmsTotal = FarmsTotal, total = total)
}

## ---- PATHS -----------------------------------------------------------------
universe_dir <- file.path(project_dir, "UniverseFiles_revised")

paths <- list(
  crop_codes  = file.path(universe_dir, "CropCodes.xlsx"),
  state_codes = file.path(universe_dir, "StateCodes.xlsx"),
  sales_data  = file.path(universe_dir, "CSA_WeightFile1.txt"),

  init_weight_dir  = file.path(universe_dir, "2020 FarmTrak weighting files", "Initial Weight Files 2020"),
  final_weight_dir = file.path(universe_dir, "2020 FarmTrak weighting files", "Final Weight Files 2020"),
  universe_excels  = file.path(universe_dir, "2020 FarmTrak weighting files", "Universe Excels 2020"),

  acres_by_crd_xlsx = file.path(universe_dir, "2020 FarmTrak weighting files", "Universe Excels 2020", "2020 acres by CRD.xlsx"),
  farms_by_size_xlsx= file.path(universe_dir, "2020 FarmTrak weighting files", "Universe Excels 2020", "2020 Farms by Size.xlsx"),
  addl_univ_xlsx    = file.path(universe_dir, "2020 FarmTrak weighting files", "Universe Excels 2020", "AditionalUniverses.xlsx"),

  sales_universe_xlsx = file.path(universe_dir, "FarmTrakCP_ClientSalesUniverses_2020_20220601.xlsx")
)


stop_if_missing(paths$crop_codes, "CropCodes.xlsx")
stop_if_missing(paths$state_codes, "StateCodes.xlsx")
stop_if_missing(paths$sales_data, "CSA_WeightFile1.txt")
stop_if_missing(paths$acres_by_crd_xlsx, "2020 acres by CRD.xlsx")
stop_if_missing(paths$farms_by_size_xlsx, "2020 Farms by Size.xlsx")
stop_if_missing(paths$sales_universe_xlsx, "FarmTrakCP_ClientSalesUniverses_2020_20220601.xlsx")

## ---- LOAD REFERENCE TABLES -------------------------------------------------
CropCodes  <- readxl::read_xlsx(paths$crop_codes)
StateCodes <- readxl::read_xlsx(paths$state_codes)

state <- StateCodes$StateCode[StateCodes$StateDescription == stateNm]
if (length(state) != 1) stop("State not found or ambiguous in StateCodes.xlsx: ", stateNm)

# Load universe support tables once (faster than reading inside crop loop)
AcresByCRD_all <- readxl::read_xlsx(paths$acres_by_crd_xlsx)
colnames(AcresByCRD_all)[3:4] <- c("CurrYr", "Acres")

FarmsBySize <- readxl::read_xlsx(paths$farms_by_size_xlsx)

SalesUniverse <- readxl::read_xlsx(paths$sales_universe_xlsx, sheet = "Data")

## ---- LOAD SALES DATA -------------------------------------------------------
SalesData <- read.csv(paths$sales_data, sep = " ", header = FALSE)
colnames(SalesData) <- c(
  "ID", "State", "Crop", "Acres",
  "Herbicide", "Insecticide", "Fungicide",
  "BASF", "BAYER", "CORTEVA AGRISCIENCE", "FMC", "SYNGENTA", "Others"
)

SalesData$FarmID <- SalesData$ID %/% 1000

data <- SalesData %>% filter(State == state)

cropIDs <- sort(unique(data$Crop))

CropCodesvec <- CropCodes$CropDescription
names(CropCodesvec) <- CropCodes$CropCode
cropNames <- CropCodesvec[paste(cropIDs)]
print(cropNames)

## ---- BUILD CROP-LEVEL DESIGN / AUX COLUMNS --------------------------------
AcresByCRD0 <- data.frame()  # collects AcresByCRD rows used per crop

for (cropcnt in seq_along(cropIDs)) {
  message("[", cropcnt, "/", length(cropIDs), "] crop = ", cropNames[cropcnt])

  CropNm <- cropNames[cropcnt]
  CropID <- cropIDs[cropcnt]

  cropfile <- file.path(paths$init_weight_dir,  paste0("fct", CropNm, "1.txt"))
  cropfile_final <- file.path(paths$final_weight_dir, paste0(CropNm, "Acres.txt"))

  stop_if_missing(cropfile,       paste0("Initial weight file for ", CropNm))
  stop_if_missing(cropfile_final, paste0("Final weight file for ", CropNm))

  AcresData <- read.csv(cropfile, sep = "", header = FALSE)
  WeightsData <- read.csv(cropfile_final, sep = "", header = FALSE)

  AcresData <- AcresData %>% select(where(function(x) all(!is.na(x))))

  colnames(AcresData)[1:5] <- c("FarmID", "State", "CRD", "FarmSize", "Acres")
  colnames(WeightsData) <- c("FarmID", "W1", "W2")

  AcresData <- left_join(AcresData, WeightsData, by = "FarmID")
  AcresData <- AcresData %>% filter(State == state)

  # --- Universe acres by CRD for this crop/state ---
  AcresByCRD <- AcresByCRD_all %>%
    filter(Crop == CropNm, substr(CurrYr, start = 1, stop = 2) == stateNm)

  Curr <- AcresByCRD$CurrYr

  # If only state-level universe row exists, synthesize CRD list from sample
  if (length(Curr) == 1) {
    Curr <- paste(sort(unique(AcresData$CRD)), collapse = ", ")
    AcresByCRD$CurrYr <- paste(stateNm, "- CRD ", Curr)
    AcresByCRD$CRD <- Curr[1]
  } else {
    Curr_tmp <- substr_right(AcresByCRD$CurrYr, 5)
    Curr <- substring(AcresByCRD$CurrYr, 10)
    Curr <- Curr[Curr_tmp != "Delta"]
    AcresByCRD <- AcresByCRD[Curr_tmp != "Delta", ]
  }

  # Harmonize CRD coding to match universe grouping
  for (cnt in seq_along(Curr)) {
    tmpvec <- tryCatch(as.integer(strsplit(Curr[cnt], ",")[[1]]), warning = function(w) w)
    if (inherits(tmpvec, "warning")) {
      message("CRD parse failed for CropNm=", CropNm)
      stop(tmpvec)
    }

    crd0 <- tmpvec[1]
    AcresByCRD$CRD[cnt] <- crd0

    for (crd in tmpvec) {
      AcresData$CRD[AcresData$CRD == crd] <- crd0
    }
  }

  # Initial weights within CRD
  AcresByCRD <- merge(AcresByCRD, AcresData %>% count(CRD))
  AcresByCRD$InitW <- AcresByCRD$Farms / AcresByCRD$n
  AcresData <- merge(AcresData, AcresByCRD %>% select(CRD, InitW))

  # --- Build calibration columns: Acres by CRD (crop-specific) ---
  if (nrow(AcresByCRD) > 1) {
    crd <- as.factor(AcresData$CRD)
    modelmat <- model.matrix(~ -1 + crd)
    colnames(modelmat) <- sapply(sort(unique(AcresByCRD$CRD)), function(x) sprintf(paste0(CropNm, "%d"), x))
    AcresData <- cbind(AcresData, modelmat * AcresData$Acres)
  } else {
    AcresData <- cbind(AcresData, AcresData$Acres)
    colnames(AcresData)[ncol(AcresData)] <- paste0(CropNm, "10")
  }

  # --- Build calibration columns: Farm counts by size (crop-specific) ---
  if (length(unique(AcresData$FarmSize)) > 1) {
    farmsize <- as.factor(AcresData$FarmSize)
    modelmat <- model.matrix(~ -1 + farmsize)
    colnames(modelmat) <- sapply(sort(unique(AcresData$FarmSize)), function(x) sprintf(paste0(CropNm, "_%d"), x))
    AcresData <- cbind(AcresData, modelmat)
  } else {
    AcresData <- cbind(AcresData, 1)
    colnames(AcresData)[ncol(AcresData)] <- paste0(CropNm, "_", unique(AcresData$FarmSize))
  }

  # --- Special-case auxiliary columns (kept from original) ---
  if (!is.null(AcresData$V6) &&
      ((stateNm == "MN" && CropNm == "SweetCorn") ||
       (stateNm == "WI" && CropNm == "SweetCorn") ||
       (stateNm == "WI" && CropNm == "SnapBeans") ||
       (stateNm == "IL" && CropNm == "Pumpkin") ||
       (stateNm == "WI" && CropNm == "Peas") ||
       (stateNm == "MN" && CropNm == "Peas"))) {

    V6 <- as.factor(AcresData$V6)
    modelmat <- model.matrix(~ -1 + V6)
    colnames(modelmat) <- sapply(sort(unique(AcresData$V6)), function(x) sprintf(paste0("Proc", CropNm, "%d"), x))
    AcresData <- cbind(AcresData, modelmat * AcresData$Acres)
  }

  if (!is.null(AcresData$V6) && stateNm == "CA" && CropNm == "Grapes") {
    names(AcresData)[names(AcresData) == "V6"] <- "RaisinGrapes"
    names(AcresData)[names(AcresData) == "V7"] <- "WineGrapes"
    names(AcresData)[names(AcresData) == "V8"] <- "TableGrapes"
  }

  if (!is.null(AcresData$V6) && CropNm == "Soybean" && stateNm %in% c("MO", "TN", "AR", "LA", "MS")) {
    V6 <- as.factor(AcresData$V6)
    modelmat <- model.matrix(~ -1 + V6)
    DeltaFarms <- modelmat
    colnames(DeltaFarms) <- c("DeltaFarms", "NonDeltaFarms")
    DeltaAcres <- modelmat * AcresData$Acres
    colnames(DeltaAcres) <- c("DeltaAcres", "NonDeltaAcres")
    AcresData <- cbind(AcresData, DeltaFarms, DeltaAcres)
  }

  AcresData$Crop <- CropID

  # --- Impute FarmID for any mismatched points (kept from original) ---
  Acres1 <- AcresData$Acres[!(AcresData$FarmID %in% data$FarmID)]
  Acres2 <- data$Acres[(data$Crop == CropID) & !(data$FarmID %in% AcresData$FarmID)]

  if (length(Acres1) == length(Acres2) && length(Acres1) > 0) {
    AcresData$FarmID[!(AcresData$FarmID %in% data$FarmID)][order(Acres1)] <-
      data$FarmID[data$Crop == CropID][!(data$FarmID[data$Crop == CropID] %in% AcresData$FarmID)][order(Acres2)]
  }

  AcresByCRD0 <- rbind(AcresByCRD0, AcresByCRD)

  # Merge crop-specific aux/design cols into the main record-level dataset
  data <- left_join(data, AcresData, by = c("State", "Crop", "Acres", "FarmID"))
}

## ---- CLEANUP: COLLAPSE DUPLICATE COLUMNS ----------------------------------
# After the repeated left_join(), some columns appear as InitW.x / InitW.y, etc.
# Keep your original approach (row-sum across duplicates).

data <- collapse_dup_cols(data, pattern = "^InitW",   new_name = "InitW")
data <- collapse_dup_cols(data, pattern = "^FarmSize", new_name = "FarmSize")
data <- collapse_dup_cols(data, pattern = "^CRD",     new_name = "CRD")

# Drop any records with zero initial weight
if (any(data$InitW == 0, na.rm = TRUE)) {
  message("Dropping ", sum(data$InitW == 0, na.rm = TRUE), " rows with InitW==0")
}
data <- data[data$InitW != 0, ]


## ---- DESIGN / STRATA INFO (FOR SEs) ---------------------------------------
# We use (Crop, FarmSize) strata + FarmsBySize.xlsx to define N_h and n_h,
# then compute design weights d = N_h / n_h and main.R-style variance estimates.

data <- preprocess_farmsize_for_variance(data, stateNm)

design <- build_design_from_farms_by_size(
  data = data,
  FarmsBySize = FarmsBySize,
  CropCodes = CropCodes,
  stateNm = stateNm
)

d_design <- design$d
stratum_id <- design$stratum_id
strata_params <- design$strata_params
N_total <- design$N_total

## ---- BUILD UNIVERSE TOTALS -------------------------------------------------
# Acres totals (by crop x CRD)
AcresByCRD1 <- AcresByCRD0 %>% select(CRD, Acres, Crop)
AcresTotal <- AcresByCRD1$Acres
names(AcresTotal) <- paste0(AcresByCRD1$Crop, AcresByCRD1$CRD)

# Additional universes (kept from original)
if (stateNm %in% c("MN", "WI", "IL")) {
  AdditionalUniv <- readxl::read_xlsx(paths$addl_univ_xlsx, sheet = "Processor crops") %>%
    filter(State == stateNm)

  AddlTotal <- AdditionalUniv$Acres
  names(AddlTotal) <- paste0("Proc", AdditionalUniv$Crop, AdditionalUniv$Processor)
  AcresTotal <- c(AcresTotal, AddlTotal)
}

if (stateNm %in% c("CA")) {
  AdditionalUniv <- readxl::read_xlsx(paths$addl_univ_xlsx, sheet = "GrapeCABreakout")
  AddlTotal <- as.numeric(unlist(AdditionalUniv)[3:5])
  names(AddlTotal) <- names(unlist(AdditionalUniv)[3:5])
  AcresTotal <- c(AcresTotal, AddlTotal)
}

if (stateNm %in% c("MO", "TN", "AR", "LA", "MS")) {
  AdditionalUniv <- readxl::read_xlsx(paths$addl_univ_xlsx, sheet = "SoybeanDeltaRegions") %>%
    filter(State == stateNm)

  AddlTotal <- as.numeric(unlist(AdditionalUniv)[5:6])
  names(AddlTotal) <- names(unlist(AdditionalUniv)[5:6])
  AcresTotal <- c(AcresTotal, AddlTotal)
}

# Farms totals (by crop x size)
FarmsBySize1 <- melt((FarmsBySize %>% filter(State == stateNm))[, c(-2, -3, -ncol(FarmsBySize))])
FarmsBySize1 <- FarmsBySize1[!is.na(FarmsBySize1$value), ]

FarmsTotal <- FarmsBySize1$value
names(FarmsTotal) <- paste0(FarmsBySize1$Crop, "_", FarmsBySize1$variable)

if (stateNm %in% c("MO", "TN", "AR", "LA", "MS")) {
  AdditionalUniv <- readxl::read_xlsx(paths$addl_univ_xlsx, sheet = "SoybeanDeltaRegions") %>%
    filter(State == stateNm)

  AddlTotal <- as.numeric(unlist(AdditionalUniv)[3:4])
  names(AddlTotal) <- names(unlist(AdditionalUniv)[3:4])
  FarmsTotal <- c(FarmsTotal, AddlTotal)
}

# Sales universe totals (truth for comparison)
SalesUniverse1 <- SalesUniverse %>%
  filter(State == state) %>%
  select(PesticideType, CompanyDesc, UniverseExpenditures_Actual)

SalesUniverse11 <- SalesUniverse1 %>%
  group_by(PesticideType) %>%
  summarize(sum = sum(UniverseExpenditures_Actual), .groups = "drop")

SalesUniverse12 <- SalesUniverse1 %>%
  group_by(CompanyDesc) %>%
  summarize(sum = sum(UniverseExpenditures_Actual), .groups = "drop")

salestotal1 <- SalesUniverse11$sum
names(salestotal1) <- SalesUniverse11$PesticideType

salestotal2 <- SalesUniverse12$sum
names(salestotal2) <- SalesUniverse12$CompanyDesc

# Keep only totals that exist as columns in the data
salestotal1 <- salestotal1[names(salestotal1) %in% names(data)]
salestotal2 <- salestotal2[names(salestotal2) %in% names(data)]
AcresTotal  <- AcresTotal[names(AcresTotal) %in% names(data)]

# Combine all totals into one named vector
#   - We'll calibrate only on a subset (Acres + Farms), but use sales totals for evaluation.
total <- c(salestotal1, salestotal2, AcresTotal, FarmsTotal)

## ---- BUILD X MATRIX FOR CALIBRATION ---------------------------------------
# Base/design weights (from Crop x FarmSize strata)
d <- d_design

# Xs includes columns that have known totals OR are needed (missing FarmTotal cols set to 0)
Xs <- data[names(total)[names(total) %in% names(data)]]
Xs[names(FarmsTotal)[!(names(FarmsTotal) %in% names(Xs))]] <- 0
Xs <- Xs %>% replace(is.na(.), 0)

# Merge small farm-size cells for stability
collapsed <- collapse_sparse_farmsize_cells(
  Xs = Xs,
  FarmsTotal = FarmsTotal,
  total = total,
  cropNames = cropNames,
  min_n = min_cell_n
)

Xs <- collapsed$Xs
FarmsTotal <- collapsed$FarmsTotal
total <- collapsed$total

## ---- CHOOSE CALIBRATION VARIABLES -----------------------------------------
# Default: calibrate on (Acres totals + Farms-by-size totals)
calib_vars <- c(names(AcresTotal), names(FarmsTotal))
# calib_vars <- c(names(AcresTotal))
calib_vars <- calib_vars[calib_vars %in% names(Xs)]

Xs_calib <- as.matrix(Xs[, calib_vars, drop = FALSE])
X_total  <- total[calib_vars]

# Add intercept so regression/calibration includes a constant term (main.R-style)
Xs_calib <- cbind(Intercept = 1, Xs_calib)
X_total  <- c(Intercept = N_total, X_total)

# N_total already computed from (Crop, FarmSize) strata above


## ---- ESTIMATOR COMPARISON (SALES TOTALS) ----------------------------------
# Estimators included:
#   - HT
#   - REG        (OLS GREG)
#   - SREG       (sample-split OLS GREG; K-fold cross-fit)
#   - REG.Spline (GAM spline over CRD, by crop; model-assisted)
#   - SREG.Spline (sample-split GAM)
#   - REG.RF     (Random forest via ranger; model-assisted)
#   - SREG.RF    (sample-split RF)
#
# Targets: pesticide type totals + company totals (truth from SalesUniverse)

# Choose targets (truth from SalesUniverse)
targets_pesticide <- names(salestotal1)
if (!is.null(targets_of_interest)) {
  targets_pesticide <- intersect(targets_pesticide, targets_of_interest)
  if (length(targets_pesticide) == 0) {
    warning("No requested targets_of_interest found in pesticide totals; using all pesticide types.")
    targets_pesticide <- names(salestotal1)
  }
}

y_targets <- targets_pesticide
if (isTRUE(include_company_targets)) {
  y_targets <- c(y_targets, names(salestotal2))
}

# One shared K-fold assignment across all targets (fair comparison)
fold_id <- make_folds(n = nrow(data), K = ss_K, seed = ss_seed)

# --- Prepare predictors for nonparametric estimators (Crop x CRD) -----------
# We restrict nonparametric predictors to variables for which we have a
# population distribution: Crop x CRD with known Farms counts from AcresByCRD.

# Sample-side predictors
CropCodesvec_all <- CropCodes$CropDescription
names(CropCodesvec_all) <- as.character(CropCodes$CropCode)

data_np <- data %>%
  mutate(
    CropNm = factor(CropCodesvec_all[as.character(Crop)]),
    CRD_fac = factor(as.character(CRD))
  )

# Universe-side cells (Crop name in AcresByCRD0 is the Crop description)
pop_cells <- AcresByCRD0 %>%
  transmute(Crop = as.character(Crop),
            CRD = as.numeric(CRD),
            Farms = as.numeric(Farms)) %>%
  group_by(Crop, CRD) %>%
  summarize(Farms = max(Farms, na.rm = TRUE), .groups = "drop")

# Align factor levels (union of sample and universe)
crop_lvls <- sort(unique(c(levels(data_np$CropNm), pop_cells$Crop)))
crd_lvls <- sort(unique(c(levels(data_np$CRD_fac), as.character(pop_cells$CRD))))

data_np <- data_np %>%
  mutate(
    CropNm = factor(as.character(CropNm), levels = crop_lvls),
    CRD_fac = factor(as.character(CRD_fac), levels = crd_lvls),
    CRD_num = suppressWarnings(as.numeric(as.character(CRD_fac)))
  )

pop_np <- data.frame(
  CropNm = factor(pop_cells$Crop, levels = crop_lvls),
  CRD_fac = factor(as.character(pop_cells$CRD), levels = crd_lvls),
  Farms = as.numeric(pop_cells$Farms)
)
pop_np$CRD_num <- suppressWarnings(as.numeric(as.character(pop_np$CRD_fac)))

# Pre-allocate results
res_list <- vector("list", length(y_targets))

for (j in seq_along(y_targets)) {
  yname <- y_targets[j]
  y <- data[[yname]]
  y[is.na(y)] <- 0

  true_total <- total[[yname]]

  ## ---- Point estimates ----
  est_ht <- estimate_total_ht(y, d)

  # OLS REG / SREG (main regressors from calibration matrix)
  reg_fit <- estimate_total_greg(
    y, X_s = Xs_calib, d_s = d, X_total = X_total,
    wls = ss_wls, return_resid = TRUE
  )
  est_reg <- reg_fit$total

  sreg_fit <- estimate_total_sreg(
    y, X_s = Xs_calib, d_s = d, X_total = X_total,
    fold_id = fold_id, wls = ss_wls, return_resid = TRUE
  )
  est_sreg <- sreg_fit$total

  # Nonparametric spline (GAM) on (CropNm, CRD)
  gam_fit <- estimate_total_gam_greg(
    y = y, df_s = data_np[, c("CropNm", "CRD_fac", "CRD_num")], d_s = d,
    pop_df = pop_np[, c("CropNm", "CRD_fac", "CRD_num", "Farms")],
    wls = ss_wls, gam_k = gam_k,
    seed = ss_seed + 100L + 10L * j,
    return_resid = TRUE
  )
  est_spline <- gam_fit$total

  sreg_gam_fit <- estimate_total_gam_sreg(
    y = y, df_s = data_np[, c("CropNm", "CRD_fac", "CRD_num")], d_s = d,
    pop_df = pop_np[, c("CropNm", "CRD_fac", "CRD_num", "Farms")],
    fold_id = fold_id,
    wls = ss_wls, gam_k = gam_k,
    seed = ss_seed + 1000L + 10L * j,
    return_resid = TRUE
  )
  est_sreg_spline <- sreg_gam_fit$total

  # Nonparametric random forest (ranger) on (CropNm, CRD)
  rf_fit <- estimate_total_rf_greg(
    y = y, df_s = data_np[, c("CropNm", "CRD_fac")], d_s = d,
    pop_df = pop_np[, c("CropNm", "CRD_fac", "Farms")],
    wls = ss_wls,
    num_trees = rf_num_trees, mtry = rf_mtry, min_node_size = rf_min_node_size,
    num_threads = rf_num_threads, verbose = rf_verbose,
    seed = ss_seed + 20000L + 10L * j,
    return_resid = TRUE
  )
  est_rf <- rf_fit$total

  sreg_rf_fit <- estimate_total_rf_sreg(
    y = y, df_s = data_np[, c("CropNm", "CRD_fac")], d_s = d,
    pop_df = pop_np[, c("CropNm", "CRD_fac", "Farms")],
    fold_id = fold_id,
    wls = ss_wls,
    num_trees = rf_num_trees, mtry = rf_mtry, min_node_size = rf_min_node_size,
    num_threads = rf_num_threads, verbose = rf_verbose,
    seed = ss_seed + 30000L + 10L * j,
    return_resid = TRUE
  )
  est_sreg_rf <- sreg_rf_fit$total

  ## ---- Standard errors (design-based) ----
  # HT: use y itself
  se_ht <- sqrt(varhat_stratified(e = y, stratum_id = stratum_id, strata_params = strata_params))

  # OLS REG: residual linearization
  se_reg <- sqrt(varhat_stratified(e = reg_fit$resid, stratum_id = stratum_id, strata_params = strata_params))

  # OLS SREG: cross-fitted residuals (main.R style)
  se_sreg <- sqrt(varhat_stratified(e = sreg_fit$resid, stratum_id = stratum_id, strata_params = strata_params))

  # GAM / RF: residual linearization; SREG versions: cross-fitted residuals
  se_spline <- if (all(is.na(gam_fit$resid))) NA_real_ else
    sqrt(varhat_stratified(e = gam_fit$resid, stratum_id = stratum_id, strata_params = strata_params))

  se_sreg_spline <- if (all(is.na(sreg_gam_fit$resid))) NA_real_ else
    sqrt(varhat_stratified(e = sreg_gam_fit$resid, stratum_id = stratum_id, strata_params = strata_params))

  se_rf <- if (all(is.na(rf_fit$resid))) NA_real_ else
    sqrt(varhat_stratified(e = rf_fit$resid, stratum_id = stratum_id, strata_params = strata_params))

  se_sreg_rf <- if (all(is.na(sreg_rf_fit$resid))) NA_real_ else
    sqrt(varhat_stratified(e = sreg_rf_fit$resid, stratum_id = stratum_id, strata_params = strata_params))

  res_list[[j]] <- data.frame(
    target = yname,
    true_total = as.numeric(true_total),

    # HT = as.numeric(est_ht),
    # HT_se = as.numeric(se_ht),

    GREG = as.numeric(est_reg),
    GREG_se = as.numeric(se_reg),

    SREG = as.numeric(est_sreg),
    SREG_se = as.numeric(se_sreg),

    GREG.Spline = as.numeric(est_spline),
    GREG.Spline_se = as.numeric(se_spline),

    SREG.Spline = as.numeric(est_sreg_spline),
    SREG.Spline_se = as.numeric(se_sreg_spline),

    GREG.RF = as.numeric(est_rf),
    GREG.RF_se = as.numeric(se_rf),

    SREG.RF = as.numeric(est_sreg_rf),
    SREG.RF_se = as.numeric(se_sreg_rf)
  )
}

res_est <- do.call(rbind, res_list)

# Add relative errors vs known universe totals (for quick diagnostics)
estimator_cols <- c("GREG", "SREG", "GREG.Spline", "SREG.Spline", "GREG.RF", "SREG.RF")
for (nm in estimator_cols) {
  if (!nm %in% names(res_est)) next
  res_est[[paste0(nm, "_relerr")]] <- (res_est[[nm]] / res_est$true_total) - 1
}

print(res_est)

## ---- VISUALIZE ESTIMATOR PERFORMANCE (FOCUS TARGETS) ----------------------
# By default, we visualize only Herbicide / Fungicide / Insecticide.
# Uses point estimates and (design-based) SEs computed above.

make_plots <- TRUE
make_extra_plots <- TRUE  # set FALSE to show only the totals-by-estimator plot

# Targets to plot
plot_targets <- targets_of_interest

if (make_plots) {
  estimators <- estimator_cols
  se_cols <- paste0(estimators, "_se")

  # Restrict plotting set
  res_plot <- res_est
  if (!is.null(plot_targets)) {
    keep <- res_plot$target %in% plot_targets
    if (!any(keep)) {
      warning("No plot_targets found in res_est; plotting all available targets.")
    } else {
      res_plot <- res_plot[keep, , drop = FALSE]
    }
  }

  # Order targets: keep user-specified order when available
  if (!is.null(plot_targets)) {
    target_order <- plot_targets[plot_targets %in% res_plot$target]
  } else {
    target_order <- res_plot$target[order(res_plot$true_total, decreasing = TRUE)]
  }

  # long format (estimate)
  df_est <- res_plot %>%
    select(target, true_total, all_of(estimators)) %>%
    pivot_longer(cols = all_of(estimators), names_to = "estimator", values_to = "estimate")

  # long format (se)
  df_se <- res_plot %>%
    select(target, all_of(se_cols)) %>%
    pivot_longer(cols = all_of(se_cols), names_to = "estimator", values_to = "se") %>%
    mutate(estimator = sub("_se$", "", estimator))

  df_long <- left_join(df_est, df_se, by = c("target", "estimator")) %>%
    filter(!is.na(estimate))

  df_long <- df_long %>%
    mutate(
      estimator = factor(estimator, levels = estimators),
      target = factor(target, levels = target_order)
    )

  z <- qnorm(0.975)
  df_long <- df_long %>%
    mutate(
      ci_low = estimate - z * se,
      ci_high = estimate + z * se,
      ratio = estimate / true_total,
      relerr = ratio - 1,
      ci_width_rel = (ci_high - ci_low) / true_total,
      covered = (ci_low <= true_total) & (ci_high >= true_total)
    )

  truth_df <- res_plot %>%
    select(target, true_total) %>%
    distinct() %>%
    mutate(target = factor(target, levels = target_order))

  # --- Figure 1: Totals by estimator (with 95% CI) -------------------------
  p_totals <- ggplot(df_long, aes(x = estimator, y = estimate)) +
    geom_hline(data = truth_df, aes(yintercept = true_total), linetype = "dashed", color = "red") +
    geom_point() +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15) +
    # facet_wrap(~ target, scales = "free_y") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      # title = paste0(stateNm, ": Estimated SALES totals by estimator (K = ", ss_K, ")"),
      title = stateNm,
      # subtitle = "Points are estimates; bars are 95% CI; red dashed line is the benchmark total",
      x = NULL,
      y = "Estimated Total Herbicide"
    )

  print(p_totals)

  # --- Optional extra figures (still limited to plot_targets) --------------
  if (isTRUE(make_extra_plots)) {

    # Figure 2: Relative error (%)
    p_relerr <- ggplot(df_long, aes(x = estimator, y = 100 * relerr)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point() +
      facet_wrap(~ target, scales = "free_y") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = paste0(stateNm, ": Relative error by estimator (K = ", ss_K, ")"),
        subtitle = "(Estimate / Truth - 1) × 100%",
        x = NULL,
        y = "Relative error (%)"
      )

    # print(p_relerr)

    # Figure 3: CI coverage heatmap (does 95% CI include truth?)
    p_cov <- ggplot(df_long, aes(x = estimator, y = target, fill = covered)) +
      geom_tile(color = "white") +
      theme_bw() +
      labs(
        title = paste0(stateNm, ": 95% CI coverage (K = ", ss_K, ")"),
        subtitle = "TRUE means the 95% CI contains the SalesUniverse total",
        x = NULL,
        y = NULL
      )

    # print(p_cov)

    # Figure 4: Relative CI width (CI length / truth)
    p_ciwidth <- ggplot(df_long, aes(x = estimator, y = ci_width_rel)) +
      geom_point() +
      facet_wrap(~ target, scales = "free_y") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = paste0(stateNm, ": Relative CI width by estimator (K = ", ss_K, ")"),
        subtitle = "(CI_high - CI_low) / Truth",
        x = NULL,
        y = "Relative CI width"
      )

    # print(p_ciwidth)
  }
}
