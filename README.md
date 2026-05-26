# DREG: Debiased Regression Estimation for High-Dimensional Survey Data

R code accompanying the paper on debiased / cross-fit regression estimators for
finite-population totals under complex survey designs.

The simulation study compares the following estimators of a population total under
stratified and rejective (Poisson-type) sampling:

- `HT` — Horvitz–Thompson
- `Diff` — difference estimator with the true mean function
- `GREG` (and `GREG.oracle`) — generalized regression with full / oracle covariates
- `GREG.Lasso` — GREG with a CV-tuned lasso working model
- `SREG`, `SREG.Lasso` — sample-split (K-fold cross-fit) versions of the above

The K-sensitivity and nonlinear tasks additionally include nonparametric working
models (regression splines and random forests).

## Repository layout

```
DREG-main/
├── main.R                          # Main simulation script (grid / K / nonlinear)
├── kynetec_calib_organized_v10.R   # Applied calibration workflow (Kynetec data)
├── Results/                        # Pre-computed simulation outputs
├── Archive/                        # Earlier drafts, slides, notes
├── LICENSE                         # MIT
└── README.md
```

## Requirements

- R (≥ 4.1 recommended)
- CRAN packages: `mvtnorm`, `glmnet`, `ggplot2`, `foreach`, `doParallel`, `doRNG`,
  `sampling` (plus `mgcv` and `ranger` for the nonlinear task)

Install with:

```r
install.packages(c("mvtnorm","glmnet","ggplot2","foreach","doParallel",
                   "doRNG","sampling","mgcv","ranger"))
```

## Reproducing the simulations

`main.R` is driven by command-line key/value flags and is written to be
PSOCK/SLURM-friendly. From a shell:

```bash
# Main grid figures (vary p, vary r); OLS working model; 500 reps
Rscript main.R --task=grid --simnum=500

# WLS variants
Rscript main.R --task=grid --fit=wls --simnum=500

# Sensitivity to the number of folds K (fixed p=90, r=-0.75)
Rscript main.R --task=K --K_grid=2,5,10,20 --simnum=500

# Nonlinear DGP + spline / RF working models
Rscript main.R --task=nonlinear --simnum=200 --K=5 \
               --np_dim=5 --spline_df=4 --rf_trees=200

# Everything at once
Rscript main.R --task=all --fit=both --simnum=500
```

Useful flags: `--design={all,stratified,rejective}`, `--vary={all,p,r}`,
`--outdir=fig`, `--workers=<int>`, `--cluster={psock,fork,auto}`,
`--seed_pop=2`, `--seed_sim=11`. Outputs (`.png` figures and `.rds` results)
are written to `fig/` by default.

## Applied analysis

`kynetec_calib_organized_v10.R` is a self-contained workflow that builds
calibration weights and compares point estimators (HT, GREG, sample-split
GREG, and spline / random-forest variants) on Kynetec crop-survey data.
Edit the `project_dir` and `stateNm` variables in the USER CONFIG block at
the top before running.

## License

MIT — see `LICENSE`.

## Citation

If you use this code, please cite the accompanying paper (citation to be added
upon publication).
