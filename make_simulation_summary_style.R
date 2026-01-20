################################################################################
# make_simulation_summary_style.R
#
# What this script does
#  1) Loads ONLY the function definitions from your main.R (no driver execution)
#  2) Runs 4 scenarios:
#       - Stratified / noninformative  (r = 0)
#       - Stratified / informative     (r = -0.75)
#       - Rejective  / noninformative  (r = 0)
#       - Rejective  / informative     (r = -0.75)
#  3) Builds a “screenshot-style” 2x2 figure:
#       (A) Bias distribution boxplots
#       (B) Relative RMSE heatmap
#       (C) Mean 95% CI width heatmap
#       (D) Coverage probability (%) heatmap
#  4) Saves PNG/PDF WITHOUT patchwork (avoids the guide alignment error)
#
# How to run:
#   Put this script in the same folder as main.R, then:
#     source("make_simulation_summary_style.R")
################################################################################

# ----------------------------
# Packages (master)
# ----------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(foreach)
  library(doParallel)
  library(doRNG)
  library(mvtnorm)
  library(grid)      # safe 2x2 saving (no patchwork)
})

# Needed by your estimator/sampling functions
suppressPackageStartupMessages({
  library(glmnet)
  library(sampling)
})

# ----------------------------
# User controls
# ----------------------------
SIMNUM   <- 50   # start with 10/50. Increase to 100/500 once stable.

N        <- 1000
n        <- 300
p        <- 90
s        <- 5
seq_K    <- 10

r_noninf <- 0
r_inf    <- -0.75

workers  <- max(1L, min(parallel::detectCores() - 1L, 8L, SIMNUM))

out_png  <- "simulation_summary_style.png"
out_pdf  <- "simulation_summary_style.pdf"

# ----------------------------
# Helper: format numbers like the screenshot
# ----------------------------
fmt_trim <- function(x, digits = 2) {
  out <- formatC(x, format = "f", digits = digits)
  out <- sub("0+$", "", out)
  out <- sub("\\.$", "", out)
  out
}

# ----------------------------
# 1) Load ONLY the function block from main.R (avoid running driver)
#
# Your main.R has:
#   ## -------- Utilities --------
#   ... (functions)
#   ## Parallel setup (shared)
# We parse ONLY the Utilities block.
# ----------------------------
source_main_functions <- function(path = "main.R") {
  lines <- readLines(path, warn = FALSE)
  
  start <- grep("^##\\s*--------\\s*Utilities\\s*--------\\s*$", lines)[1]
  if (is.na(start)) stop("Could not find: '## -------- Utilities --------' in main.R")
  
  end <- grep("^##\\s*Parallel setup\\s*\\(shared\\)\\s*$", lines)[1]
  if (is.na(end)) stop("Could not find: '## Parallel setup (shared)' in main.R")
  
  code <- paste(lines[start:(end - 1)], collapse = "\n")
  eval(parse(text = code), envir = .GlobalEnv)
  
  # Required objects we will export to workers
  EXPORT_VARS <- c(
    "ar1_cor",
    "sample_design_stratified", "quad_Omega_strata",
    "make_pi_from_z", "rejective_sample", "sample_design_rejective",
    "var_poisson_approx",
    "run_estimators"
  )
  
  missing <- EXPORT_VARS[!vapply(EXPORT_VARS, exists, logical(1), envir = .GlobalEnv)]
  if (length(missing)) {
    stop("Loaded Utilities block but missing: ", paste(missing, collapse = ", "))
  }
  
  message("Loaded Utilities block from main.R. Will export to workers: ",
          paste(EXPORT_VARS, collapse = ", "))
  
  EXPORT_VARS
}

EXPORT_VARS <- source_main_functions("main.R")

# ----------------------------
# 2) Convert replicate matrices -> tidy long
# ----------------------------
make_long_from_mats <- function(res_mat, se_mat, truth, scheme_label) {
  if (is.null(res_mat) || nrow(res_mat) == 0 || ncol(res_mat) == 0)
    stop("res_mat empty: all reps failed or returned no columns.")
  if (is.null(colnames(res_mat)) || is.null(colnames(se_mat)))
    stop("Missing colnames on result matrices; cannot label estimators.")
  
  est_df <- as.data.frame(res_mat, check.names = FALSE) %>% mutate(rep = row_number())
  se_df  <- as.data.frame(se_mat,  check.names = FALSE) %>% mutate(rep = row_number())
  
  est_long <- est_df %>% pivot_longer(-rep, names_to = "estimator", values_to = "estimate")
  se_long  <- se_df  %>% pivot_longer(-rep, names_to = "estimator", values_to = "se")
  
  est_long %>%
    left_join(se_long, by = c("rep","estimator")) %>%
    mutate(truth = truth, scheme = scheme_label)
}

# ----------------------------
# 3) Run one scenario with PSOCK workers
#    Key fix: explicitly export functions to workers.
# ----------------------------
run_one_scenario <- function(SIMNUM,
                             sampling = c("stratified","rejective"),
                             r_value,
                             N, n, p, s, seq_K,
                             seed_pop = 2,
                             seed_reps = 11,
                             workers = 4,
                             scheme_label = NULL,
                             EXPORT_VARS) {
  
  sampling <- match.arg(sampling)
  if (is.null(scheme_label)) scheme_label <- paste0(sampling, "_r=", r_value)
  
  # Reduce oversubscription
  Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
  
  # ----- Population (fixed within scenario; same across scenarios if seed_pop identical)
  set.seed(seed_pop)
  len <- round(N/4)
  
  X0 <- mvtnorm::rmvnorm(N, mean = rep(0, p), sigma = ar1_cor(p, 0.2)) + 2
  e  <- rnorm(N)
  
  beta <- c(rep(1, s), rep(0, p - s))
  mu   <- X0 %*% beta
  y    <- drop(mu + e)
  t_y  <- sum(y)
  
  X <- cbind(1, X0)
  
  z <- r_value * e + rnorm(N, 0, sqrt(1 - r_value^2))
  z_sorted_idx <- order(z)
  
  # ----- Cluster
  cl <- parallel::makeCluster(workers, type = "PSOCK")
  doParallel::registerDoParallel(cl)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  
  # Packages on workers
  parallel::clusterEvalQ(cl, {
    library(glmnet)
    library(sampling)
    NULL
  })
  
  # EXPORT the required functions to workers (THIS FIXES your error)
  parallel::clusterExport(cl, varlist = EXPORT_VARS, envir = .GlobalEnv)
  
  doRNG::registerDoRNG(seed_reps)
  
  final_res <- foreach(simnum = 1:SIMNUM, .errorhandling = "pass") %dopar% {
    if (sampling == "stratified") {
      samp <- sample_design_stratified(N=N, n=n, len=len, z_sorted_idx=z_sorted_idx, y=y, X=X)
      out  <- run_estimators(
        sampling = sampling,
        N=N, n=n, s=s, y=y, X=X, mu=mu,
        Index=samp$Index, d1=samp$d1, d1_s=samp$d1_s, y_s=samp$y_s, X_s=samp$X_s,
        n_h=samp$n_h, st=samp$st, cumn_h=samp$cumn_h, len=samp$len,
        seq_K=round(seq_K), wls=FALSE
      )
    } else {
      samp <- sample_design_rejective(N=N, n=n, z=z, y=y, X=X, gamma=1)
      out  <- run_estimators(
        sampling = sampling,
        N=N, n=n, s=s, y=y, X=X, mu=mu,
        Index=samp$Index, pi=samp$pi, d1=samp$d1, d1_s=samp$d1_s, y_s=samp$y_s, X_s=samp$X_s,
        seq_K=round(seq_K), wls=FALSE
      )
    }
    
    list(y_res = out$y_res[-1], se_res = out$sigma_res[-1])
  }
  
  ok <- vapply(final_res, function(x) is.list(x) && all(c("y_res","se_res") %in% names(x)), logical(1))
  if (!any(ok)) {
    errs <- final_res[!ok]
    msg  <- paste0(capture.output(str(errs[[1]])), collapse = "\n")
    stop("All reps failed in ", scheme_label, ". Example error object:\n", msg)
  }
  
  res_mat <- do.call(rbind, lapply(final_res[ok], `[[`, "y_res"))
  se_mat  <- do.call(rbind, lapply(final_res[ok], `[[`, "se_res"))
  
  make_long_from_mats(res_mat, se_mat, truth = t_y, scheme_label = scheme_label)
}

# ----------------------------
# 4) Build 4 panels (A,B,C,D) as ggplot objects
# ----------------------------
make_summary_panels <- function(sim_long,
                                target_estimator = "Diff",
                                estimator_levels = NULL,
                                scheme_levels = NULL) {
  
  stopifnot(all(c("rep","scheme","estimator","estimate","se","truth") %in% names(sim_long)))
  
  dat <- sim_long %>%
    mutate(
      estimator = as.character(estimator),
      scheme    = as.character(scheme)
    )
  
  if (!is.null(scheme_levels)) {
    dat$scheme <- factor(dat$scheme, levels = scheme_levels)
  } else {
    dat$scheme <- factor(dat$scheme)
  }
  
  # Handle estimator ordering safely: keep only levels that exist
  if (!is.null(estimator_levels)) {
    keep <- estimator_levels[estimator_levels %in% unique(dat$estimator)]
    if (length(keep) == 0) keep <- sort(unique(dat$estimator))
    dat$estimator <- factor(dat$estimator, levels = keep)
    estimator_order <- keep
  } else {
    dat$estimator <- factor(dat$estimator)
    estimator_order <- levels(dat$estimator)
  }
  
  dat <- dat %>%
    mutate(
      error    = estimate - truth,
      ci_width = 2 * 1.96 * se,
      covered  = abs(error) <= 1.96 * se
    )
  
  # (A) Bias distribution
  pA <- ggplot(dat, aes(x = scheme, y = error, fill = estimator)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.7, color = "deeppink2") +
    geom_boxplot(width = 0.65, outlier.shape = NA,
                 position = position_dodge(width = 0.75)) +
    # geom_point(aes(group = estimator),
    #            position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
    #            size = 1.1, alpha = 0.55, color = "grey35") +
    labs(title = "(A) Estimation Bias", fill = "Method") +
    theme_bw(base_size = 12) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold")
    )
  
  # Summaries for heatmaps
  sum_df <- dat %>%
    group_by(scheme, estimator) %>%
    summarise(
      rmse      = sqrt(mean(error^2)),
      mean_ci_w = mean(ci_width),
      cpct      = 100 * mean(covered),
      .groups   = "drop"
    )
  
  # Relative RMSE
  rrmse_df <- sum_df %>%
    group_by(scheme) %>%
    mutate(rrmse = rmse / rmse[estimator == target_estimator][1]) %>%
    ungroup()
  
  heatmap_plot <- function(df, value_col, title, legend_title,
                           digits = 2, limits = NULL,
                           low = "white", high = "orange") {
    ggplot(df, aes(x = scheme, y = estimator, fill = .data[[value_col]])) +
      geom_tile(color = "white") +
      geom_text(aes(label = fmt_trim(.data[[value_col]], digits = digits)), size = 3) +
      scale_fill_gradient(low = low, high = high, name = legend_title, limits = limits) +
      scale_y_discrete(limits = rev(estimator_order)) +
      labs(title = title) +
      theme_minimal(base_size = 12) +
      theme(
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.title = element_text(face = "bold")
      )
  }
  
  pB <- heatmap_plot(
    rrmse_df, "rrmse",
    "(B) Relative Root Mean Square Error\n(compared to the target estimator)",
    "RRMSE",
    digits = 2,
    low = "orange", high = "white"
  )
  
  pC <- heatmap_plot(
    sum_df, "mean_ci_w",
    "(C) 95% Confidence Interval Width",
    "CI width",
    digits = 2,
    low = "slateblue", high = "white"
  )
  
  pD <- heatmap_plot(
    sum_df, "cpct",
    "(D) Coverage Probability (%)\n(nominal level: 95%)",
    "CP%",
    digits = 1,
    limits = c(50, 100),
    low = "white", high = "purple3"
  )
  
  list(pA = pA, pB = pB, pC = pC, pD = pD)
}

# ----------------------------
# 5) Save 2x2 safely (NO patchwork)
# ----------------------------
save_2x2_grid_png <- function(filename, panels, width = 16, height = 5.8, dpi = 300) {
  grDevices::png(filename, width = width, height = height, units = "in", res = dpi)
  grid::grid.newpage()
  lay <- grid::grid.layout(nrow = 2, ncol = 2)
  grid::pushViewport(grid::viewport(layout = lay))
  
  print(panels$pA, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(panels$pC, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(panels$pB, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(panels$pD, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
  
  grDevices::dev.off()
}

save_2x2_grid_pdf <- function(filename, panels, width = 16, height = 5.8) {
  grDevices::pdf(filename, width = width, height = height)
  grid::grid.newpage()
  lay <- grid::grid.layout(nrow = 2, ncol = 2)
  grid::pushViewport(grid::viewport(layout = lay))
  
  print(panels$pA, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(panels$pC, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(panels$pB, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(panels$pD, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
  
  grDevices::dev.off()
}

# ----------------------------
# 6) Run the 4 schemes
# ----------------------------
scheme_levels <- c(
  "Strat / noninfo",
  "Strat / info",
  "Reject / noninfo",
  "Reject / info"
)

# Use your paper order if these names match your run_estimators() output
estimator_levels <- c("Diff","GREG.oracle","GREG","SREG","GREG.Lasso","SREG.Lasso")

message("Running simulations: SIMNUM=", SIMNUM, "  workers=", workers)

schemes <- bind_rows(
  run_one_scenario(SIMNUM, "stratified", r_noninf,
                   N=N, n=n, p=p, s=s, seq_K=seq_K, workers=workers, EXPORT_VARS=EXPORT_VARS,
                   scheme_label = "Strat / noninfo"),
  run_one_scenario(SIMNUM, "stratified", r_inf,
                   N=N, n=n, p=p, s=s, seq_K=seq_K, workers=workers, EXPORT_VARS=EXPORT_VARS,
                   scheme_label = "Strat / info"),
  run_one_scenario(SIMNUM, "rejective", r_noninf,
                   N=N, n=n, p=p, s=s, seq_K=seq_K, workers=workers, EXPORT_VARS=EXPORT_VARS,
                   scheme_label = "Reject / noninfo"),
  run_one_scenario(SIMNUM, "rejective", r_inf,
                   N=N, n=n, p=p, s=s, seq_K=seq_K, workers=workers, EXPORT_VARS=EXPORT_VARS,
                   scheme_label = "Reject / info")
)

# ----------------------------
# 7) Build panels and save
# ----------------------------
panels <- make_summary_panels(
  schemes,
  target_estimator = "Diff",
  estimator_levels = estimator_levels,
  scheme_levels    = scheme_levels
)

save_2x2_grid_png(out_png, panels, width = 16, height = 5.8, dpi = 300)
save_2x2_grid_pdf(out_pdf, panels, width = 16, height = 5.8)

message("Saved: ", out_png)
message("Saved: ", out_pdf)

# Optional: show a panel in the R plotting window
print(panels$pA)
print(panels$pB)
print(panels$pC)
print(panels$pD)
