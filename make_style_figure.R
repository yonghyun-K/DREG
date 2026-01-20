suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(patchwork)
  library(foreach)
  library(doParallel)
  library(doRNG)
  library(mvtnorm)
  library(glmnet)
  library(sampling)
})

# ---- 1) Load your functions without running the driver ----
options(RUN_DRIVER = FALSE)
source("main.R")   # expects you've added the guard described above

# ---- 2) Small helper: pretty numeric labels like the screenshot ----
fmt_trim <- function(x, digits = 2) {
  out <- formatC(x, format = "f", digits = digits)
  out <- sub("0+$", "", out)
  out <- sub("\\.$", "", out)
  out
}

# ---- 3) Plot builder: 2x2 panel (A)-(D) ----
make_summary_figure <- function(sim_long,
                                target_estimator = "SREG",
                                estimator_levels = NULL,
                                scheme_levels = NULL,
                                method_colors = NULL) {
  
  stopifnot(all(c("rep","scheme","estimator","estimate","se","truth") %in% names(sim_long)))
  
  dat <- sim_long %>%
    mutate(
      estimator = as.character(estimator),
      scheme    = as.character(scheme)
    )
  
  if (!is.null(estimator_levels)) dat$estimator <- factor(dat$estimator, levels = estimator_levels)
  if (!is.null(scheme_levels))    dat$scheme    <- factor(dat$scheme, levels = scheme_levels)
  
  dat <- dat %>%
    mutate(
      error    = estimate - truth,
      ci_width = 2 * 1.96 * se,
      covered  = abs(error) <= 1.96 * se
    )
  
  # (A) boxplots of estimation error ("bias distribution")
  pA <- ggplot(dat, aes(x = scheme, y = error, fill = estimator)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.7, color = "deeppink2") +
    geom_boxplot(
      width = 0.65, outlier.shape = NA,
      position = position_dodge(width = 0.75)
    ) +
    geom_point(
      aes(group = estimator),
      position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
      size = 1.1, alpha = 0.55, color = "grey35"
    ) +
    labs(title = "(A) Estimation Bias", fill = "Method") +
    theme_bw(base_size = 12) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold")
    )
  
  if (!is.null(method_colors)) pA <- pA + scale_fill_manual(values = method_colors)
  
  # Summaries for heatmaps
  sum_df <- dat %>%
    group_by(scheme, estimator) %>%
    summarise(
      rmse      = sqrt(mean(error^2)),
      mean_ci_w = mean(ci_width),
      cpct      = 100 * mean(covered),
      .groups   = "drop"
    )
  
  # Relative RMSE vs target within each scheme
  rrmse_df <- sum_df %>%
    group_by(scheme) %>%
    mutate(rrmse = rmse / rmse[estimator == target_estimator][1]) %>%
    ungroup()
  
  estimator_order <- if (!is.null(estimator_levels)) estimator_levels else levels(factor(dat$estimator))
  
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
  
  # (B) Relative RMSE
  pB <- heatmap_plot(
    rrmse_df, "rrmse",
    "(B) Relative Root Mean Square Error\n(compared to the target estimator)",
    "RRMSE",
    digits = 2,
    low = "white", high = "orange"
  )
  
  # (C) Mean CI width
  pC <- heatmap_plot(
    sum_df, "mean_ci_w",
    "(C) 95% Confidence Interval Width",
    "CI width",
    digits = 2,
    low = "white", high = "slateblue"
  )
  
  # (D) Coverage %
  pD <- heatmap_plot(
    sum_df, "cpct",
    "(D) Coverage Probability (%)\n(nominal level: 95%)",
    "CP%",
    digits = 1,
    limits = c(0, 100),
    low = "white", high = "purple3"
  )
  
  # Layout matches your screenshot: A|C on top, B|D bottom
  (pA | pC) / (pB | pD)
}

# ---- 4) Scenario runner: one (sampling, r) combination -> replicate-level tidy data ----
run_one_scenario <- function(SIMNUM,
                             sampling = c("stratified","rejective"),
                             r_value,
                             p = 90,
                             N = 1000, n = 300, s = 5,
                             seq_K = 10,
                             seed_pop = 2,
                             seed_sampling = 123,
                             seed_reps = 11,
                             workers = max(1L, parallel::detectCores() - 1L),
                             scheme_label = NULL) {
  
  sampling <- match.arg(sampling)
  
  # --- Population (kept fixed across scenarios if seed_pop fixed) ---
  set.seed(seed_pop)
  len <- round(N / 4)
  
  X0 <- mvtnorm::rmvnorm(
    n = N,
    mean = rep(0, p),
    sigma = ar1_cor(p, 0.2)   # uses your function from main.R
  ) + 2
  
  e <- rnorm(N, 0, 1)
  
  X_base <- X0[, 1:p, drop = FALSE]
  beta <- c(rep(1, s), rep(0, p - s))
  mu <- X_base %*% beta
  y <- drop(mu + e)
  t_y <- sum(y)
  
  X <- cbind(1, X_base)  # add intercept (same as your code)
  
  # --- Informativeness via z (depends on r) ---
  set.seed(seed_sampling)
  z <- r_value * e + rnorm(n = N, mean = 0, sd = sqrt(1 - r_value^2))
  z_sorted_idx <- order(z)
  
  if (is.null(scheme_label)) {
    scheme_label <- paste0(
      ifelse(sampling == "stratified", "Stratified", "Rejective"),
      " (r=", r_value, ")"
    )
  }
  
  # --- Parallel setup ---
  cl <- parallel::makeCluster(workers, type = "PSOCK")
  doParallel::registerDoParallel(cl)
  on.exit({ try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)
  
  doRNG::registerDoRNG(seed_reps)
  
  final_res <- foreach(
    simnum = 1:SIMNUM,
    .errorhandling = "pass",
    .packages = c("glmnet","sampling")
  ) %dopar% {
    
    if (sampling == "stratified") {
      samp <- sample_design_stratified(
        N = N, n = n, len = len, z_sorted_idx = z_sorted_idx, y = y, X = X
      )
      
      out <- run_estimators(
        sampling = sampling,
        N=N, n=n, s=s, y=y, X=X, mu=mu,
        Index=samp$Index, d1=samp$d1, d1_s=samp$d1_s, y_s=samp$y_s, X_s=samp$X_s,
        n_h=samp$n_h, st=samp$st, cumn_h=samp$cumn_h, len=samp$len,
        seq_K=round(seq_K), wls = FALSE
      )
      
    } else {
      samp <- sample_design_rejective(N = N, n = n, z = z, y = y, X = X, gamma = 1)
      
      out <- run_estimators(
        sampling = sampling,
        N=N, n=n, s=s, y=y, X=X, mu=mu,
        Index=samp$Index, pi=samp$pi, d1=samp$d1, d1_s=samp$d1_s, y_s=samp$y_s, X_s=samp$X_s,
        seq_K=round(seq_K), wls = FALSE
      )
    }
    
    list(y_res = out$y_res, se_res = out$sigma_res)
  }
  
  ok <- vapply(final_res, function(x) is.list(x) && all(c("y_res","se_res") %in% names(x)), logical(1))
  final_res <- final_res[ok]
  
  res_mat <- do.call(rbind, lapply(final_res, function(x) x$y_res))
  se_mat  <- do.call(rbind, lapply(final_res, function(x) x$se_res))
  
  # Convert to tidy long
  est_df <- as.data.frame(res_mat) %>% mutate(rep = row_number())
  se_df  <- as.data.frame(se_mat)  %>% mutate(rep = row_number())
  
  long_est <- est_df %>%
    pivot_longer(-rep, names_to = "estimator", values_to = "estimate")
  long_se  <- se_df %>%
    pivot_longer(-rep, names_to = "estimator", values_to = "se")
  
  long_est %>%
    left_join(long_se, by = c("rep","estimator")) %>%
    mutate(truth = t_y, scheme = scheme_label)
}

# -------------------------
# RUN: 4 schemes = 2 sampling × 2 informativeness
# -------------------------
SIMNUM  <- 500
p_fixed <- 90

r_noninf <- 0
r_inf    <- -0.75   # your requested informative setting

schemes <- bind_rows(
  run_one_scenario(SIMNUM, "stratified", r_noninf, p = p_fixed, scheme_label = "Stratified / noninformative (r=0)"),
  run_one_scenario(SIMNUM, "stratified", r_inf,    p = p_fixed, scheme_label = "Stratified / informative (r=-0.75)"),
  run_one_scenario(SIMNUM, "rejective",  r_noninf, p = p_fixed, scheme_label = "Rejective / noninformative (r=0)"),
  run_one_scenario(SIMNUM, "rejective",  r_inf,    p = p_fixed, scheme_label = "Rejective / informative (r=-0.75)")
)

# Order axes (recommended)
scheme_levels <- c(
  "Stratified / noninformative (r=0)",
  "Stratified / informative (r=-0.75)",
  "Rejective / noninformative (r=0)",
  "Rejective / informative (r=-0.75)"
)

# Put your SREG estimators where you want them in the y-axis order
estimator_levels <- c("HT","Diff","GREG.oracle","GREG","SREG","GREG.Lasso","SREG.Lasso")

fig <- make_summary_figure(
  schemes,
  target_estimator  = "SREG",        # <-- change if you want a different reference
  estimator_levels  = estimator_levels,
  scheme_levels     = scheme_levels
)

ggsave("simulation_summary_style.png", fig, width = 16, height = 5.8, dpi = 300)
print(fig)
