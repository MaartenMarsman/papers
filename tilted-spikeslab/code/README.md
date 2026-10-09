# Determinant-tilted spike-and-slab prior for Gaussian graphical models: code

Simulation code, verification scripts, archived results, and figure and table generators for

> *The determinant-tilted spike-and-slab prior for the Gaussian graphical model: construction, properties, and Bayesian inference* (Marsman et al., in preparation).

The manuscript is not part of this repository. Section and label names below (for example `sec:hier-ratio`, `fig:cone-interior`, checks V1 to V8) refer to the manuscript.

## Contents

1. [Overview](#overview)
2. [Layout](#layout)
3. [Requirements](#requirements)
4. [Running a script](#running-a-script)
5. [The sampler](#the-sampler) (with a minimal example)
6. [Run order and dependencies](#run-order-and-dependencies)
7. [Script index](#script-index)
8. [Where each manuscript figure, table, and quoted number comes from](#where-each-manuscript-figure-table-and-quoted-number-comes-from)
9. [Reading the verification output](#reading-the-verification-output)
10. [Known issues and limitations](#known-issues-and-limitations)
11. [Reproducibility notes](#reproducibility-notes)

## Overview

The prior distribution places an entry-wise discrete spike-and-slab prior on the precision matrix $K$ of a Gaussian graphical model (point mass at zero for absent edges, a $N(0, \sigma^2)$ slab for present edges, an exponential prior with rate $\beta$ on the diagonal), restricts it to the positive definite cone, and multiplies it by the determinant tilt $|K|^{\delta}$. In the standardized parametrization of the manuscript the two scales enter only through $\eta = \beta\sigma$; every script sets `sigma = 1` and passes `beta = eta`.

Inference uses a within-model row-block Gibbs sweep (conjugate, positive definiteness automatic), between-model edge toggles along the permuted-Cholesky edge curve with coefficients read from $\Sigma = K^{-1}$ in constant time, and, for the hierarchical specification, a coupled Monte Carlo estimator of the ratio of normalizing constants $Z(\Gamma^-)/Z(\Gamma^+)$. Two specifications of the graph layer are supported: the *joint* specification (one global constant, exactly computable, bounded sparsening of the graph prior) and the *hierarchical* specification (per-graph constants, keeps the stated graph prior, needs the ratio estimator).

## Layout

```
R/scripts/     analysis, verification, and figure scripts (run from this directory)
R/src/         C++ sources, compiled on the fly with Rcpp::sourceCpp()
R/results/     archived outputs (.rds) and console logs (.log) of the archived runs
figures/       figures (PDF) produced by the fig_*.R scripts
tables/        LaTeX tables produced by the scripts
```

## Requirements

- R 4.x with a working C++ toolchain (Xcode command line tools on macOS, Rtools on Windows, gcc on Linux). The archived runs were made on an Apple silicon Mac; package versions were not recorded.
- Packages used by the current scripts: `Rcpp` and `RcppArmadillo` (every C++ source), `ggplot2`, `patchwork`, `viridisLite`, `scales` (attached by `fig_style.R`, which every figure script sources), `gridExtra` (some figure scripts). `parallel` and `grid` are part of base R. The exploratory script `delta_default_investigation.R` additionally needs the `bgms` package.

```r
install.packages(c("Rcpp", "RcppArmadillo", "ggplot2", "patchwork", "viridisLite", "scales", "gridExtra"))
```

The long-running sweeps call `parallel::mclapply(..., mc.cores = 12L)` with the core count written inside the call (only `monotonicity_check.R` has a top-level variable). On Windows, which lacks forking, change those calls to `mc.cores = 1L`.

## Running a script

Scripts resolve paths relative to a project root that defaults to the working directory, so run them from this directory:

```sh
cd papers/tilted-spikeslab/code
Rscript R/scripts/verify_sampler_port.R
```

To run from elsewhere, set the environment variable `SPIKESLAB_ROOT` to this directory. Each script sets its random seed, writes to `R/results/`, `figures/`, or `tables/`, and overwrites the archived file of the same name. The archived `.log` files were produced by redirecting the console, `Rscript R/scripts/<name>.R > R/results/<name>.log 2>&1`; not every script was run that way, so some results have no log.

## The sampler

### C++ implementation: `R/src/spikeslab_sampler.cpp`

The sampler used in every study. It maintains $\Sigma = K^{-1}$ along the chain, refreshes it by rank-two Sherman–Morrison–Woodbury updates after each row update and each accepted toggle, and recomputes it exactly every `refresh_every` sweeps (and after a failed Cholesky factorization or a diagonal slack below $10^{-6}$). The header comment of the source file documents every argument and every returned field. In brief:

- `sampler_chain_cpp(S, n_obs, q, delta, beta, sigma, p_inc, spec, log_z, n_iter, n_burn, refresh_every, track_masks, n_aux = 0, stat_thin = 0, time_blocks = FALSE)` runs the complete chain. `S` is the scatter matrix $Y^\top Y$ (the zero matrix with `n_obs = 0` for a prior-level run). `spec` is `0` for the joint specification, `1` for the hierarchical specification with an exact table `log_z` of per-graph constants (length $2^{q(q-1)/2}$, feasible for $q \le 5$), and `2` for the hierarchical specification with the coupled estimator, which then needs `n_aux > 0` auxiliary draws per candidate. Pass `log_z = numeric(0)` for `spec` 0 and 2. The chain starts at the empty graph with $K = ((\delta+1)/\beta) I$ and does not return draws of $K$; the returned `edge_freq` holds the post-burn-in inclusion frequency of each edge in its upper triangle, which is the posterior edge-inclusion probability.
- `akm_ratio_cpp(G_minus, i, j, delta, beta, sigma, N, coupled = TRUE, return_draws = FALSE)` is the coupled estimate of $Z(\Gamma^-)/Z(\Gamma^+)$ for adding edge $(i, j)$ to `G_minus`, with `N` draws per arm.

Only the exponential diagonal (Gamma shape one) is implemented; the independence proposal for other shapes in the manuscript's appendix has no code here. `track_masks = TRUE` and `spec = 1` index graphs by a 64-bit mask and are refused for $q(q-1)/2 > 62$.

Minimal example on simulated data:

```r
library(Rcpp)
sourceCpp("R/src/spikeslab_sampler.cpp")
set.seed(1)
q <- 8; n <- 50
Ktrue <- diag(q); Ktrue[cbind(1:(q-1), 2:q)] <- Ktrue[cbind(2:q, 1:(q-1))] <- 0.4
Y <- matrix(rnorm(n * q), n, q) %*% chol(solve(Ktrue))
S <- crossprod(Y)
# joint specification, eta = 1 (beta = 1, sigma = 1), tilt delta = 1, p_inc = 0.2
fit <- sampler_chain_cpp(S, n, q, delta = 1, beta = 1, sigma = 1, p_inc = 0.2,
                         spec = 0L, log_z = numeric(0), n_iter = 20000L, n_burn = 2000L,
                         refresh_every = 50L, track_masks = FALSE)
round(fit$edge_freq, 2)          # posterior edge-inclusion probabilities (upper triangle)
fit$E_m; fit$acc_rate            # posterior mean edge count; toggle acceptance rate
# hierarchical specification with 200 auxiliary draws per candidate
fitH <- sampler_chain_cpp(S, n, q, 1, 1, 1, 0.2, spec = 2L, log_z = numeric(0),
                          n_iter = 5000L, n_burn = 500L, refresh_every = 50L,
                          track_masks = FALSE, n_aux = 200L)
```

### R reference implementation: `R/scripts/sampler_reference.R`

A pure R implementation of the same algorithm, written for correctness before speed. It defines `sample_row_block()`, `toggle_edge()`, and `run_chain()`, recomputes $\Sigma$ by `solve()` at every step, and accepts an initial graph. It is sourced, not run: the verification scripts use `run_chain()` as the R-level chain (checks V4 and V5) and `sample_row_block()` as the fixed-graph prior sampler (V6, V7, and the hierarchical calibration studies).

## Run order and dependencies

Most scripts are self-contained. The dependencies are:

```
monotonicity_check.cpp ──> joint_spec_sparsening_J1.R ──> R/results/joint_spec_sparsening_J1_z.rds
                                (brute-force constants at q = 4, 5; about 10 min on 12 cores)

joint_spec_sparsening_J1_z.rds ──> verify_sampler_prior_level.R (V4 R, V5)
                               ──> verify_sampler_port.R        (V4 C++)
                               ──> verify_posterior_enumerated_V6.R (V6) ──> hier_plugin_decision_B.R
                               ──> verify_sbc_V7.R (V7)
                               ──> verify_akm_ratio.R
                               ──> joint_spec_sparsening_J2.R (reads the J1 summaries)
cone_interior_diagnostics_expanded.R ──> cone_interior_expanded.rds ──> fig_cone_interior_expanded.R
delta_scaling_regime_sweep.R ──> delta_scaling_regime_sweep.rds ──> fig_delta_scaling_regime.R
marginal_offdiag_sweep.R ──> marginal_offdiag.rds ──┬─> fig_realized_conditional_slab.R
partial_correlation_sweep.R ──> partial_correlation.rds ─┴─> fig_specified_vs_realized.R
joint_spec_sparsening_J1_chain_table.R ──> tables/joint_sparsening.tex
timing_scaling.R ──> tables/timing_scaling.tex
fig_style.R ──> every fig_*.R script (sourced)
```

All inputs are archived, so any script can be run on its own against the archived inputs. Recorded run times (12 cores): J1 about 10 minutes; the regime sweep about 2 minutes; the compiled chain runs at about 14,000 iterations per second at $q = 20$ and 1,000 at $q = 50$ at the prior level; V8 chains take seconds to minutes per chain depending on `N` (see its log). The perfect-sampler sweeps (`marginal_offdiag_sweep.R`, `partial_correlation_sweep.R`) and V6 are long-running (hours); their run times were not recorded.

## Script index

Listed by the manuscript section they support. Inputs are archived results the script reads; outputs are what it writes.

### Verification at small dimension (checks V1 to V8)

| Check | Script | Inputs | Outputs |
|---|---|---|---|
| V1 exhaustive monotonicity check | `monotonicity_check.R` (uses `R/src/monotonicity_check.cpp`) | | `monotonicity_check.rds`, `.log` |
| V2 three-estimator agreement | `verify_conditional_gaussian_route.R` | | `verify_conditional_gaussian_route.rds`, `.log` |
| V3 prior row-block Gibbs sampler against importance sampling | the validation gate inside `delta_scaling_regime_sweep.R` | | `delta_scaling_regime_sweep.rds`, `.log` |
| V4 complete chain at the prior level (R reference) and V5 within-model update at the posterior level | `verify_sampler_prior_level.R` | `joint_spec_sparsening_J1_z.rds` | `verify_sampler_prior_level.rds`, `.log` |
| V4 complete chain at the prior level (C++ port) | `verify_sampler_port.R` | `joint_spec_sparsening_J1_z.rds` | `verify_sampler_port.rds`, `.log` |
| V6 complete chain against the enumerated posterior | `verify_posterior_enumerated_V6.R` | `joint_spec_sparsening_J1_z.rds` | `verify_posterior_enumerated_V6.rds`, `.log` |
| V7 simulation-based calibration | `verify_sbc_V7.R` | `joint_spec_sparsening_J1_z.rds` | `verify_sbc_V7.rds`, `.log` |
| V8 hierarchical plug-in chain | `verify_estimated_ratio_chain_V8.R` | (the V4 reference values are written into the script) | `verify_estimated_ratio_chain_V8.rds`, `.log` |
| Ratio estimator against analytic and brute-force values | `verify_akm_ratio.R` | `joint_spec_sparsening_J1_z.rds` | `verify_akm_ratio.rds`, `.log` |

### Prior geometry and tilt calibration

| Purpose | Script | Inputs | Outputs |
|---|---|---|---|
| Loss of mass and boundary concentration of the untilted construction (data of figure `fig:cone-interior`; Gamma(3, 3) diagonal, as the caption says) | `cone_interior_diagnostics_expanded.R` (uses `R/src/Z_prior_diag.cpp`) | | `cone_interior_expanded.rds` |
| Figure `fig:cone-interior` | `fig_cone_interior_expanded.R` | `cone_interior_expanded.rds` | `figures/cone_interior.pdf` |
| Tilt required to hold the median smallest eigenvalue, by regime; contains the V3 gate | `delta_scaling_regime_sweep.R` | | `delta_scaling_regime_sweep.rds`, `.log` |
| Figure `fig:delta-scaling`, the inversion to the required tilt, and the console numbers for $c(0)$ and the partial-correlation dispersion quoted in the monotonicity section | `fig_delta_scaling_regime.R` | `delta_scaling_regime_sweep.rds` | `figures/delta_scaling_regime.pdf`, `delta_scaling_regime_inversion.rds` |
| Exact draws of an included entry and its rest-block geometry, conditional on the edge being present (perfect rejection sampler) | `marginal_offdiag_sweep.R` (uses `R/src/Z_prior_gwishart_slab_k12.cpp`) | | `marginal_offdiag.rds` |
| The same with both diagonals, for the induced partial correlation | `partial_correlation_sweep.R` (uses `R/src/Z_prior_gwishart_slab_rho.cpp`) | | `partial_correlation.rds`, `partial_correlation_sweep.log` |
| Figure `fig:realized-slab` | `fig_realized_conditional_slab.R` | `marginal_offdiag.rds` | `figures/realized_conditional_slab.pdf` |
| Figure `fig:specified-realized` | `fig_specified_vs_realized.R` | `marginal_offdiag.rds`, `partial_correlation.rds` | `figures/specified_vs_realized.pdf` |
| Table of total variation distances and variance ratios quoted in the text (the figure is not used) | `fig_marginal_specified_vs_realized.R` | `marginal_offdiag.rds` | `tables/marginal_specified_vs_realized.tex`, `figures/marginal_specified_vs_realized.pdf` |
| Table of partial-correlation dispersion and tail mass quoted in the text (the figure is not used) | `fig_partial_correlation.R` | `partial_correlation.rds` | `tables/partial_correlation_prior.tex`, `figures/partial_correlation_prior.pdf` |

### Sparsening under the joint specification

| Purpose | Script | Inputs | Outputs |
|---|---|---|---|
| J1: exact consequences at $q = 4, 5$ by brute-force constants; the bracket check of all per-pair factors | `joint_spec_sparsening_J1.R` (uses `R/src/monotonicity_check.cpp`) | | `joint_spec_sparsening_J1.rds`, `joint_spec_sparsening_J1_z.rds`, `tables/joint_sparsening_enumeration.tex`, `.log` |
| Between-seed spread of the brute-force constants at $\delta = 2$ with $2 \times 10^6$ draws per graph (about 45 min on 12 cores) | `joint_spec_sparsening_J1_delta2_recheck.R` | | `j1_delta2_recheck.rds`, `.log` |
| The manuscript's sparsening table from long chain runs, with the compensated arm | `joint_spec_sparsening_J1_chain_table.R` | | `tables/joint_sparsening.tex`, `joint_spec_sparsening_J1_chain_table.rds`, `.log` |
| J2: sampled consequences at $q = 5$ to $50$ | `joint_spec_sparsening_J2.R` | `joint_spec_sparsening_J1.rds` | `joint_spec_sparsening_J2.rds`, `.log` |
| J2b: extension to $q = 100$ in the sparse regime | `joint_spec_sparsening_J2b_q100.R` | | `joint_spec_sparsening_J2b_q100.rds`, `.log` |
| Drift map of the maintained $\Sigma$ across cells | `joint_spec_sparsening_J2_driftmap.R` | | `joint_spec_sparsening_J2_driftmap.rds`, `.log` |

### Accuracy and cost of the hierarchical approximation

| Purpose | Script | Inputs | Outputs |
|---|---|---|---|
| H1: bias, RMSE, and coupling correlation of the estimated log ratio | `hier_ratio_accuracy_A.R` | | `hier_ratio_accuracy_A.rds`, `.log` |
| H4: estimator precision and coupling at $q$ up to 100 | `hier_ratio_accuracy_at_scale_A2.R` | | `hier_ratio_accuracy_at_scale_A2.rds`, `.log` |
| H2: plug-in chain against the exact-ratio chain at the posterior level | `hier_plugin_decision_B.R` | `joint_spec_sparsening_J1_z.rds`, `verify_posterior_enumerated_V6.rds` | `hier_plugin_decision_B.rds`, `.log` |
| H3: simulation-based calibration of the plug-in chain versus $N$ | `hier_sbc_vs_N_C.R` | | `hier_sbc_vs_N_C.rds`, `.log` |
| H5: the same at $q = 10$ | `hier_sbc_at_scale_C2.R` | | `hier_sbc_at_scale_C2.rds`, `.log` |

### Computational scaling

| Purpose | Script | Outputs |
|---|---|---|
| Wall-clock cost per sweep, per candidate toggle, and per accepted toggle, from the chain's block timers | `timing_scaling.R` | `timing_scaling.rds`, `.log`, `tables/timing_scaling.tex` |

### Shared helpers

- `fig_style.R`: common figure style (sizes, palette, theme); attaches `ggplot2`, `patchwork`, `viridisLite`, and `scales`.
- `sampler_reference.R`: the R reference sampler (see above).

### Exploratory scripts, not used in the manuscript

`delta_scaling_gwishart.R` and `fig_delta_scaling.R` (tilt calibration by the perfect sampler; replaced by the regime sweep), `c_sensitivity_sweep.R` and `fig_c_sensitivity.R` (sensitivity of the coefficient in $\delta = c \log q$; uses the retired rejection sampler `Z_prior_tilt_resample.cpp`), `delta_default_investigation.R` (uses the `bgms` package as the prior sampler; output not archived), `standardized_rate_geometry.R` (geometry as a function of the diagonal rate for Gamma diagonals), `rho_uniform_calibration.R` and `rho_uniform_q2_check.R` (which specification makes the partial-correlation prior approximately uniform), `joint_spec_compensation.R` (enumeration-based compensation; the manuscript uses the chain-based arm of the chain-table script).

### Superseded scripts from the earlier construction

The project started from a construction that combined the tilt with an eigenvalue floor (parameter `epsilon`) and used a rejection sampler for the prior. Both were dropped. These scripts are kept for the record and produce nothing the manuscript uses: `Z_prior_brute_force.R`, `Z_prior_inspect.R`, `tilt_LO_check.R`, `cone_interior_tilt_resample.R`, `delta_scaling_sweep.R`, `cone_interior_diagnostics.R`, `fig_cone_interior.R` (writes `figures/cone_interior_legacy.pdf`), `fig_tilt_fix.R`, `fig_tilt_fix_expanded.R`, `fig_tilt_fix_resample.R`, `fig_tail_bound.R`, `fig_tail_bound_resample.R`, `fig_marginal_offdiag.R`, `empty_graph_moments.R`, `fig_empty_moments.R`, `verify-scale-injection.R`, and the C++ sources `Z_prior_brute.cpp` and `Z_prior_tilt_resample.cpp`.

## C++ sources

| File | Content | Status |
|---|---|---|
| `spikeslab_sampler.cpp` | The sampler and the coupled ratio estimator | current |
| `monotonicity_check.cpp` | Monte Carlo estimator of $Z(\Gamma; \delta)$ with common random numbers across graphs | current (V1, J1) |
| `Z_prior_diag.cpp` | Draws from the unrestricted entry-wise construction; records the smallest eigenvalue and positive definiteness | current (`fig:cone-interior`) |
| `Z_prior_gwishart_slab.cpp` | Perfect (rejection) sampler for the tilted prior at the exponential diagonal, with the G-Wishart factor as envelope | current (exploratory scripts) |
| `Z_prior_gwishart_slab_k12.cpp`, `Z_prior_gwishart_slab_rho.cpp` | Variants that record one entry, or that entry with its two diagonals, per accepted draw, conditional on the edge being present | current (G3 sweeps) |
| `Z_prior_brute.cpp`, `Z_prior_tilt_resample.cpp` | Brute-force and rejection samplers of the tilt-plus-floor stage | superseded |

## Where each manuscript figure, table, and quoted number comes from

| Manuscript item | Producing script | Archive |
|---|---|---|
| Figure `fig:cone-interior`; Section 2 numbers on the positive-definiteness rate, smallest eigenvalue, and diagonal inflation | `cone_interior_diagnostics_expanded.R`, `fig_cone_interior_expanded.R` | `cone_interior_expanded.rds` |
| Figure `fig:delta-scaling`; required tilt by regime; median smallest eigenvalue at fixed tilt | `delta_scaling_regime_sweep.R`, `fig_delta_scaling_regime.R` | `delta_scaling_regime_sweep.rds`, `delta_scaling_regime_inversion.rds` |
| $c(0)$ values, the tilt required for a sparsening tolerance, and the partial-correlation dispersion at large tilt (monotonicity section) | console output of `fig_delta_scaling_regime.R`; $c(0)$ also in `monotonicity_check.log` | |
| Figures `fig:realized-slab`, `fig:specified-realized` | `fig_realized_conditional_slab.R`, `fig_specified_vs_realized.R` | `marginal_offdiag.rds`, `partial_correlation.rds` |
| Total variation distances and variance ratios of the realized slab | `fig_marginal_specified_vs_realized.R` | `tables/marginal_specified_vs_realized.tex` |
| Partial-correlation standard deviation and tail masses | `fig_partial_correlation.R` | `tables/partial_correlation_prior.tex` |
| Table `tab:joint-sparsening` and the compensated specification | `joint_spec_sparsening_J1_chain_table.R` | `joint_spec_sparsening_J1_chain_table.rds` |
| Bracket check of all per-pair factors; enumeration cross-check | `joint_spec_sparsening_J1.R` | `joint_spec_sparsening_J1.rds` |
| Between-seed spread of the brute-force constants at $\delta = 2$ | `joint_spec_sparsening_J1_delta2_recheck.R` | `j1_delta2_recheck.rds`, `.log` |
| J2 results up to $q = 100$ and the drift diagnostics | `joint_spec_sparsening_J2.R`, `_J2b_q100.R`, `_J2_driftmap.R` | the matching `.rds` and `.log` |
| Hierarchical results H1 to H5 | `hier_*.R` scripts | the matching `.rds` and `.log` |
| Timing table and the scaling numbers | `timing_scaling.R` | `timing_scaling.rds` |
| Verification appendix V1 to V8 | see the verification table above | the matching `.rds` and `.log` |

## Reading the verification output

Each verification script states its pass criterion in its header comment. Two scripts print a verdict: `monotonicity_check.R` prints `OVERALL: PASS` or `FAIL`, and `delta_scaling_regime_sweep.R` prints `validation PASSED` or stops with an error. The others print the distances, noise floors, and tolerances next to each other for inspection; the criteria are those stated in the manuscript's verification appendix (a pooled-chain total variation distance at or below the between-seed noise floor, agreement of moments within Monte Carlo error, uniform ranks). One archived V6 row (joint specification at $\eta = 2$) exceeds its floor by a factor 1.23, which the manuscript reports.

## Known issues and limitations

These were found in a code review on 2026-10-08 and are reflected in the manuscript text.

1. **Perfect-sampler sweeps at $q \ge 20$ (rerun 2026-10-09).** `Z_prior_gwishart_slab_k12.cpp` and `_rho.cpp` discard a proposed graph when no precision matrix is accepted within `max_inner_tries` proposals, which favors graphs with a high acceptance rate when the discarded fraction is large. The original sweeps used a cap of $10^4$; at $q \le 15$ the estimated discarded fraction under that cap is below five percent and those cells are kept. The $q = 20$ and $q = 30$ cells were rerun with `max_inner_tries = 1e6` (set `SPIKESLAB_QS=20,30` to repeat only those cells; the scripts merge them into the archived `.rds`). Kept fractions (`n_graphs_kept / n_graphs_tried`): $q = 20$: 1.00, 0.99, 0.99, 0.94 to 0.96 at $\delta = 0, 1, \tfrac12 \log q, 2$; $q = 30$: 1.00, 0.95, 0.77, 0.68 to 0.69. The $\delta = 2$ cells at $q \ge 20$ exhausted their proposal budgets ($2 \times 10^8$ at $q = 20$, $10^9$ at $q = 30$) with 1500 to 1900 on-edge draws instead of the target 2000 to 3000. No rest-block inversion failed in the rerun (`n_inv_fail = 0`; the recorders fall back to a pivoted inverse). One residual issue: in `marginal_offdiag_sweep.R` a small number of accepted draws had non-finite rest-block completions (68, 24, 8 of about 3000 at $q = 20$ for $\delta = 0, 1, 2$; 339, 229, 157 of about 2000 at $q = 30$), because a NaN acceptance weight compared as accepted. Those draws are excluded from the conditional marginal density estimates (the reported `n` counts finite draws), and a guard `if (!K.is_finite() || !std::isfinite(log_w)) continue;` was added to both recorders on 2026-10-09; the archived rerun predates the guard. The partial-correlation sweep records only the on-edge pair and has no non-finite draws. The manuscript quotes these sweeps for $q \le 20$ and uses $q = 30$ for orientation only.
2. **V6 reference size (resolved).** The archived V6 log printed "(of 15000)" while the script used $6 \times 10^4$ reference draws per graph; a rerun on 2026-10-08 with $6 \times 10^4$ draws reproduced the archived numbers to four decimals, so the archive was made with $6 \times 10^4$ and only the printed string was stale. The archive now holds the rerun.
3. **V8 chain length.** The $N = 10^4$ chains were run for $2 \times 10^4$ iterations instead of $6 \times 10^4$; the manuscript states this.
4. **Table collision (fixed).** `joint_spec_sparsening_J1.R` used to write `tables/joint_sparsening.tex`, the same file as the chain-table script that produces the manuscript's table. It now writes `tables/joint_sparsening_enumeration.tex`.
5. **Archive gaps.** `delta_default_investigation.R` has no archived output. Several scripts have no log.
6. **Sampler scope.** Only the exponential diagonal is implemented; the mask-based options require $q(q-1)/2 \le 62$; there is no option to pass an initial state to the C++ chain.

## Reproducibility notes

- Every script that draws random numbers sets `set.seed()` before the first draw, and the parallel sweeps derive per-task seeds from a base seed written in the script, so a re-run reproduces an archived result up to floating-point differences across platforms and BLAS libraries. The `fig_*.R` scripts are deterministic transformations of archived results; their PDF output depends on the graphics device (`cairo_pdf` with a fallback to `pdf`).
- Monte Carlo standard errors and between-seed noise floors are printed next to the estimates and stored in the `.rds` files.
- The same two seeds are reused across cells in several sweeps (for example 8801 and 8802 in J2); within a cell the replicates are independent, across cells they share random-number streams. This has no bearing on the reported quantities.
