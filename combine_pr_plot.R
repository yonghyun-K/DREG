# combine_pr_plot.R
#
# Build the side-by-side (vary-p | vary-r) figure used in Section 6 of the paper
# (e.g. fig/plot_pr_stratified_ols.png) from two main.R runs, with a polished
# publication style.
#
# Usage:
#   Rscript combine_pr_plot.R --p_dir <dir for vary_over=p run> \
#                             --r_dir <dir for vary_over=r run> \
#                             [--out  <output path, .png or .pdf>] \
#                             [--width 11] [--height 4]     # inches
#                             [--font serif|sans] \
#                             [--style minimal|classic] \
#                             [--palette tableau|nejm|lancet|bright|paired|okabe] \
#                             [--base_size 15]              # font size dial
#
# Each <dir> is the --out directory of a main.R run; the script reads
# <dir>/plot_df.rds, which main.R saves automatically.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

## ---- CLI ---------------------------------------------------------------
.parse_flags <- function(args) {
  out <- list(); i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (!startsWith(a, "--")) stop("Expected --flag, got: ", a)
    key <- sub("^--", "", a)
    if (i + 1L > length(args)) stop("Missing value for --", key)
    out[[key]] <- args[i + 1L]; i <- i + 2L
  }
  out
}

opts <- .parse_flags(commandArgs(trailingOnly = TRUE))
if (is.null(opts[["p_dir"]]) || is.null(opts[["r_dir"]])) {
  stop("Must supply both --p_dir and --r_dir.")
}

out_path  <- if (!is.null(opts[["out"]]))       opts[["out"]]                  else "plot_pr_combined.png"
out_w     <- if (!is.null(opts[["width"]]))     as.numeric(opts[["width"]])    else 11   # inches
out_h     <- if (!is.null(opts[["height"]]))    as.numeric(opts[["height"]])   else 4    # inches
font_fam  <- if (!is.null(opts[["font"]]))      opts[["font"]]                 else "serif"
style     <- if (!is.null(opts[["style"]]))     opts[["style"]]                else "minimal"
palette   <- if (!is.null(opts[["palette"]]))   opts[["palette"]]              else "tableau"
base_size <- if (!is.null(opts[["base_size"]])) as.numeric(opts[["base_size"]]) else 15

stopifnot(style %in% c("minimal", "classic"))
stopifnot(font_fam %in% c("serif", "sans", "mono"))

## ---- Data --------------------------------------------------------------
read_plot_df <- function(dir) {
  f <- file.path(dir, "plot_df.rds")
  if (!file.exists(f)) {
    stop("Could not find ", f,
         ". Re-run main.R (recent version) with --out ", dir, ".")
  }
  readRDS(f)
}
df_p <- read_plot_df(opts[["p_dir"]])
df_r <- read_plot_df(opts[["r_dir"]])

# Numeric ordering on x (factor preserves levels from the simulation runs).
df_p$x <- factor(df_p$x, levels = sort(unique(as.numeric(as.character(df_p$x)))))
df_r$x <- factor(df_r$x, levels = sort(unique(as.numeric(as.character(df_r$x)))))

## ---- Style -------------------------------------------------------------
# Palette registry. Pick one with --palette. All four entries are placed in
# the same order: GREG, DREG, GREG.Lasso, DREG.Lasso. The pairing strategy
# differs per palette:
#  - tableau / nejm / lancet / bright: four maximally distinct hues.
#  - paired: warm family for GREG/DREG (baseline + proposed OLS), cool family
#    for the Lasso variants. Within a family, light = OLS, dark = Lasso.
#  - okabe: the previous Okabe-Ito colorblind-safe choice.
palettes <- list(
  # Tableau 10 — the data-visualization gold standard. Polished and modern.
  tableau = c(GREG = "#4E79A7", DREG = "#F28E2B",
              GREG.Lasso = "#59A14F", DREG.Lasso = "#E15759"),
  # NEJM journal palette — vivid but mature, common in statistical papers.
  nejm    = c(GREG = "#BC3C29", DREG = "#0072B5",
              GREG.Lasso = "#E18727", DREG.Lasso = "#20854E"),
  # Lancet journal palette — slightly bolder, very legible in print.
  lancet  = c(GREG = "#00468B", DREG = "#ED0000",
              GREG.Lasso = "#42B540", DREG.Lasso = "#0099B4"),
  # Paul Tol "bright" — colorblind-safe but more saturated than Okabe-Ito.
  bright  = c(GREG = "#4477AA", DREG = "#EE6677",
              GREG.Lasso = "#228833", DREG.Lasso = "#AA3377"),
  # Semantic pairing: GREG-family warm, Lasso-family cool; OLS lighter, Lasso deeper.
  paired  = c(GREG = "#E8A33D", DREG = "#C24634",
              GREG.Lasso = "#7AB8C7", DREG.Lasso = "#2C5784"),
  # Previous default (Okabe-Ito), colorblind-safe.
  okabe   = c(GREG = "#D55E00", DREG = "#009E73",
              GREG.Lasso = "#0072B2", DREG.Lasso = "#CC79A7")
)

if (!palette %in% names(palettes)) {
  stop("Unknown --palette '", palette,
       "'. Choose one of: ", paste(names(palettes), collapse = ", "))
}
est_cols <- palettes[[palette]]

base_theme <- function(family = "serif", base_size = 15) {
  # All text sizes derive from base_size so a single dial controls the figure.
  sz_axis_text  <- base_size - 2     # tick labels
  sz_axis_title <- base_size         # axis titles (p, r)
  sz_y_title    <- base_size - 1     # y-axis description (slightly smaller)
  sz_leg_title  <- base_size - 1     # "Estimator" header
  sz_leg_text   <- base_size - 2     # legend entries
  sz_tag        <- base_size         # (a) / (b) panel tags

  thm <- theme_minimal(base_size = base_size, base_family = family) +
    theme(
      # Lose the box; keep light horizontal-only guides for value reading.
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25),
      panel.border       = element_blank(),
      # Proper axis lines: dominant horizontal/vertical structure.
      axis.line.x        = element_line(color = "grey25", linewidth = 0.45),
      axis.line.y        = element_line(color = "grey25", linewidth = 0.45),
      axis.ticks         = element_line(color = "grey25", linewidth = 0.4),
      axis.ticks.length  = unit(3.5, "pt"),
      axis.text          = element_text(color = "grey20", size = sz_axis_text),
      axis.title.x       = element_text(face = "italic", size = sz_axis_title,
                                        margin = margin(t = 8)),
      axis.title.y       = element_text(face = "plain", size = sz_y_title,
                                        margin = margin(r = 8)),
      # Compact, readable legend.
      legend.title       = element_text(face = "bold", size = sz_leg_title),
      legend.text        = element_text(size = sz_leg_text),
      legend.key.width   = unit(1.7, "lines"),
      legend.key.height  = unit(1.2, "lines"),
      legend.spacing.y   = unit(3, "pt"),
      legend.box.spacing = unit(10, "pt"),
      # Panel tag inside the plot, top-left.
      plot.tag           = element_text(face = "bold", size = sz_tag),
      plot.tag.position  = c(0.015, 0.985),
      plot.margin        = margin(10, 10, 6, 10)
    )
  if (style == "classic") {
    # Strip horizontal guide lines too — pure axes only.
    thm <- thm + theme(panel.grid.major.y = element_blank())
  }
  thm
}

make_panel <- function(df, xlab, tag = NULL, show_y_title = TRUE) {
  p <- ggplot(df,
              aes(x = x, y = value,
                  group    = interaction(estimator, metric),
                  color    = estimator,
                  linetype = metric)) +
    geom_line(linewidth = 1) +
    geom_point(data = subset(df, metric == "RMSE"), size = 2) +
    scale_color_manual(values = est_cols, name = "Estimator") +
    scale_linetype_manual(values = c(RMSE = "solid", Bias = "dashed"),
                          name = NULL,
                          breaks = c("RMSE", "Bias")) +
    labs(x = xlab,
         y = if (show_y_title) "RMSE (solid) / Bias (dashed)" else NULL) +
    base_theme(family = font_fam, base_size = base_size) +
    guides(
      linetype = guide_legend(
        order = 1,
        override.aes = list(color = "grey20", linewidth = 1)
      ),
      color = guide_legend(order = 2)
    )
  if (!is.null(tag)) p <- p + labs(tag = tag)
  p
}

p_left  <- make_panel(df_p, "p", tag = "(a)", show_y_title = TRUE)
p_right <- make_panel(df_r, "r", tag = "(b)", show_y_title = FALSE)

## ---- Combine -----------------------------------------------------------
combined <- (p_left + p_right) +
  plot_layout(guides = "collect") &
  theme(legend.position    = "right",
        legend.box.spacing = unit(8, "pt"))

## ---- Save --------------------------------------------------------------
ext <- tolower(tools::file_ext(out_path))
dev <- switch(ext,
              "png"  = "png",
              "pdf"  = "pdf",
              "tiff" = "tiff",
              "svg"  = "svg",
              "png")

# Use cairo on Linux for nicer PNG anti-aliasing of text when available.
extra <- list()
if (dev == "png" && capabilities("cairo")) extra$type <- "cairo"

do.call(ggsave,
        c(list(filename = out_path,
               plot     = combined,
               width    = out_w,
               height   = out_h,
               dpi      = 300,
               units    = "in",
               device   = dev),
          extra))

message("Wrote ", out_path,
        " (", out_w, "x", out_h, " in, ", dev,
        ", font='", font_fam, "', style='", style,
        "', palette='", palette, "', base_size=", base_size, ")")
