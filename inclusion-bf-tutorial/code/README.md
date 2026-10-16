# Inclusion Bayes factors for psychological networks: code

Analysis scripts, figure and table generators, and the archived simulation-study results for

> *Inclusion Bayes factors for psychological networks: A tutorial on evaluating evidence of edge presence and absence* (Marsman, in preparation).

The manuscript is not part of this repository. Section names below (the worked example, the robustness section, the simulation appendix) refer to the manuscript.

## Contents

1. [Overview](#overview)
2. [Layout](#layout)
3. [Requirements](#requirements)
4. [Data](#data)
5. [Running a script](#running-a-script)
6. [Run order, dependencies, and run times](#run-order-dependencies-and-run-times)
7. [Script index](#script-index)
8. [The simulation study](#the-simulation-study)
9. [Reproducibility notes](#reproducibility-notes)

## Overview

The tutorial analyzes the Wenchuan posttraumatic-stress data (17 symptoms, 362 survivors) twice: once with the estimation workflow (EBICglasso network, bootstrap intervals) and once with the Bayesian workflow of the R package `bgms`, which returns an inclusion Bayes factor for every edge and classifies each edge as present, absent, or undecided. The scripts reproduce every number, figure, and table of that analysis: the model fit and its diagnostics, the edge classification, the comparison with the estimation workflow, the robustness checks (prior sensitivity, an alternative structure prior, posterior predictive and calibration checks), and the centrality and three-symptom examples.

A second group of scripts is the known-truth simulation study behind the appendix on verdict error rates. It compares ordinary four-chain analyses with effectively exact reference analyses of the same simulated data (ordinal, binary, and two Gaussian arms) and grades diagnostic rules by how well they catch the verdicts that a longer run would overturn. The study took about a day of computing; its checkpoints are archived in this repository so that the appendix tables can be regenerated in seconds.

## Layout

```
R/                         analysis, check, and figure scripts (run from this directory)
results/threshold-study/   archived simulation-study checkpoints: manifests, simulated
                           data sets, per-fit results, run logs (shipped, 20 MB)
figures/                   figures (PDF) written by the scripts (regenerated, not shipped)
tables/                    CSV and RDS outputs written by the scripts (regenerated, not shipped)
renv.lock                  the pinned package library (160 packages)
sessionInfo.txt            the R session the manuscript's numbers were produced in
```

## Requirements

- R 4.6 (the results were produced with R 4.6.0 on macOS, Apple silicon; R >= 4.4 should work).
- `bgms` 0.2.0.0 from CRAN, plus `ggplot2`, `qgraph`, `bootnet`, `gridExtra`, `S7`, `mvtnorm`, and `MASS` (ships with R).

To restore the exact package library, run R in this directory and

```r
install.packages("renv")
renv::init(bare = TRUE)   # makes this directory an renv project (creates renv/ and .Rprofile)
renv::restore()           # installs the 160 packages pinned in renv.lock
```

then restart R. Without renv, `install.packages(c("bgms", "ggplot2", "qgraph", "bootnet", "gridExtra", "S7", "mvtnorm"))` installs current versions; the model fits are bit-reproducible only with bgms 0.2.0.0, see [Reproducibility notes](#reproducibility-notes).

`bgms` compiles against the RcppParallel version present at build time. If `library(bgms)` fails to load after a binary install, reinstall it from source: `install.packages("bgms", type = "source")`.

## Data

Both data sets ship with `bgms` and need no download:

- `data("Wenchuan", package = "bgms")`: the 17 posttraumatic-stress symptoms of the worked example.
- `data("ADHD", package = "bgms")`: the binary ADHD screener used to calibrate the binary arm of the simulation study.

## Running a script

Every script resolves paths relative to this directory, so run it from here:

```sh
cd code
Rscript R/01-worked-example.R
```

Scripts write to `figures/` and `tables/` and overwrite earlier output of the same name. Nothing is written outside this directory. Each script sets its own random seeds (see [Reproducibility notes](#reproducibility-notes)).

## Run order, dependencies, and run times

Run `R/01-worked-example.R` first. It fits the model to the Wenchuan data and saves the fit (`tables/fit-wenchuan.rds`) and the edge tables that most other scripts read. After that:

| Script | Needs | Run time |
|---|---|---|
| `R/01-worked-example.R` | nothing | 2 to 3 min (four chains, 15,000 iterations each) |
| `R/02-frequentist-comparison.R` | R/01 | seconds |
| `R/17-fig-intro-network.R` | nothing | seconds |
| `R/18-fig1-evidence-networks.R` | R/01, R/17 (shared node layout) | seconds |
| `R/19-fig-omrf-intuition.R` | nothing | seconds |
| `R/20-fig-centrality.R` | R/01 | seconds |
| `R/21-fig-spike-slab.R` | R/01 | seconds |
| `R/23-run-ppc.R` | R/01 | seconds |
| `R/24-fig-calibration.R` | R/01 | seconds |
| `R/30-ppc-diagnostic-ladder.R` | R/01 | seconds |
| `R/33-bootstrap-vs-posterior.R` | R/01 | about 15 s |
| `R/34-verdict-structure.R` | R/01 | seconds |
| `R/16-run-sensitivity-check.R` | nothing (refits the model) | 3 to 4 min (a fresh fit plus the anchored refits) |
| `R/15-fig-scale-sensitivity.R` | R/16 | seconds |
| `R/22-run-edge-prior-check.R` | R/01 (refits under the Beta-Bernoulli prior) | about 2 min |
| `R/36-robustness-reading.R` | R/16 | seconds |
| `R/38-consolidated-classification.R` | R/01, R/16, R/22 | seconds |
| `R/32-toy-structures.R` | nothing (its own fit, four million draws) | about 8 min |
| `R/35-fig-toy-structures.R` | R/32 | seconds |
| `R/09-rb-pip-prototype.R` | nothing | seconds |
| simulation-study scripts | see [The simulation study](#the-simulation-study) | seconds from the archived checkpoints; about a day from scratch |

Run times are from an Apple M-series laptop. `bgms::bgm()` runs its chains in parallel on the available cores by default.

## Script index

| Script | Reproduces |
|---|---|
| `R/01-worked-example.R` | The worked example: fits the ordinal Markov random field to the Wenchuan data; writes the fitted object (`tables/fit-wenchuan.rds`), the convergence diagnostics (`tables/diagnostics.csv`), the per-edge tables (`tables/edge-classification.csv`, `tables/edge-supplement.csv`, `tables/weight-summary.csv`, `tables/indicator-summary.csv`), and the strength-versus-evidence figure (`figures/fig3-strength-vs-evidence.pdf`). |
| `R/02-frequentist-comparison.R` | The comparison of the EBICglasso network with the Bayesian three-way classification (`tables/glasso-vs-bayes.csv`, `figures/fig-glasso-vs-bayes.pdf`), with the penalty-free checks of the disputed edges quoted in the comparison section. |
| `R/09-rb-pip-prototype.R` | The conjugate check quoted in the Rao-Blackwellization appendix: a two-predictor spike-and-slab regression with an enumerable truth, comparing the raw indicator average with the one-step Rao-Blackwellized estimator across inclusion-probability regimes (`tables/rb-pip-prototype.csv`). Standalone. |
| `R/15-fig-scale-sensitivity.R` | The prior-sensitivity figure (`figures/fig-scale-sensitivity.pdf`), drawn with the package's own plot method from the check object of R/16. |
| `R/16-run-sensitivity-check.R` | The anchored `prior_sensitivity_check()` over 0.4 to 2.5 times the slab scale (`tables/prior-sensitivity-check.rds`, `tables/prior-sensitivity-edges.csv`). |
| `R/17-fig-intro-network.R` | The introduction figure: the field's standard EBICglasso network of the Wenchuan data (`figures/fig-intro-network.pdf`) and the shared node layout (`tables/network-layout.csv`). |
| `R/18-fig1-evidence-networks.R` | Figure 1: the worked-example network split into presence, undecided, and absence panels on one layout (`figures/fig1-edge-evidence.pdf`). |
| `R/19-fig-omrf-intuition.R` | The Box 1 figure: what an edge weight of the ordinal Markov random field does to a focal item's response distribution (`figures/fig-omrf-intuition.pdf`). |
| `R/20-fig-centrality.R` | Posterior strength centrality with structure uncertainty propagated (`figures/fig-centrality.pdf`, `tables/centrality-summary.csv`). |
| `R/21-fig-spike-slab.R` | The concept figure of the Bayesian section: spike-and-slab prior and posterior for three edges spanning the three verdicts (`figures/fig-spike-slab.pdf`). |
| `R/22-run-edge-prior-check.R` | The Beta-Bernoulli structure-prior check of the robustness section: refits with `bgms::beta_bernoulli_prior(1, 1)` and cross-tabulates verdicts against the reported analysis (`tables/edge-prior-check.csv`). |
| `R/23-run-ppc.R` | The posterior predictive (sum-score) check of the robustness section (`figures/fig-ppc.pdf`, `tables/ppc-coverage.csv`). |
| `R/24-fig-calibration.R` | The per-variable PAV-adjusted calibration figure of the robustness section (`figures/fig-calibration.pdf`, `tables/calibration-summary.csv`). |
| `R/30-ppc-diagnostic-ladder.R` | The extended diagnostic ladder behind the predictive check (sufficient-statistic, decomposition, and conditional checks; console output and `figures/fig-ppc-third-order.pdf`). Not shown in the manuscript. |
| `R/32-toy-structures.R` | The three-symptom example of the concepts section: fits the model to three Wenchuan symptoms with four million posterior draws and tabulates the posterior probability of all eight candidate structures (`tables/toy-structures.csv`). The fit itself is not saved. |
| `R/33-bootstrap-vs-posterior.R` | One row of the translation table executed side by side: the bootnet edge-weight bootstrap against the posterior credible intervals for the two focal edges (`tables/bootstrap-edge-cis.csv`). Single-core, so the bootstrap is seed-reproducible. |
| `R/34-verdict-structure.R` | The numbers behind the substantive reading of the verdict panels (criterion structure of the absence panel, cluster concentration of the undecided verdicts, strongest absences and bridge candidates); console output from `tables/edge-supplement.csv`. |
| `R/35-fig-toy-structures.R` | The figure of the three-symptom example (`figures/fig-toy-structures.pdf`). |
| `R/36-robustness-reading.R` | The counts quoted in the robustness section from the prior-sensitivity output (`tables/prior-sensitivity-edges.csv`); console output. |
| `R/38-consolidated-classification.R` | The consolidated per-edge robustness classification across all checks (`tables/consolidated-robustness.csv`). |
| `R/20-threshold-study-design.R` | Simulation study, step 1: calibration constants, generating networks, simulated data sets, and the fit manifest (`results/threshold-study/`). |
| `R/21-threshold-study-pilot.R` | Simulation study, step 2: the timing pilot that projects the cost of the full grid (`--stage build` writes the pilot manifests, `--stage collect` reads their checkpoints back; `tables/threshold-study-pilot-timings.csv`, `-projection.csv`). |
| `R/22-threshold-study-run.R` | Simulation study, step 3: the checkpointed fit worker (`results/threshold-study/fits/<fit_id>.rds`, one file per fit). |
| `R/26-recompute-inclusion-bf.R` | Simulation study, consistency pass: recomputes the stored inclusion Bayes factors of every checkpoint from the per-chain posterior inclusion probabilities, so that all checkpoints are on one scale (`results/threshold-study/inclusion-bf-correction.csv`). |
| `R/23-threshold-study-analysis.R` | Simulation study, step 4: known-truth grading of the diagnostic rules on the calibration arm (`tables/threshold-study-truth.csv`, `-rules-continuous.csv`, `-rules-indicator.csv`, `-rules-weight.csv`, `-rules-nearboundary.csv`, `-catching.csv`, `-stability.csv`, `-markov.csv`, `-bf-sensitivity.csv`). |
| `R/25-threshold-study-transfer.R` | Simulation study, step 5: the stage-1 rules applied unchanged to every arm (`tables/threshold-study-transfer.csv`) and the pooled pair rule that the decision tree quotes (`tables/threshold-study-pair-rule.csv`). |
| `R/27-binary-field-gate.R` | Simulation study follow-up: the per-node field gate that selects the coupling scale of the binary rerun cells (`tables/threshold-study-binary-field-gate.csv`, `results/threshold-study/binary-field-gate-chosen.rds`). |
| `R/28-ggm-rule-search.R` | Simulation study follow-up: Gaussian-graphical-model-specific diagnostic rule candidates (`tables/threshold-study-ggm-rules.csv`). |
| `R/29-followup-design.R` | Simulation study follow-up: designs the gated binary rerun cells and the near-boundary Gaussian arm, regenerates those data sets, and extends the manifest (`tables/threshold-study-followup-design.csv`, `results/threshold-study/manifest-binary-rerun.csv`, `manifest-ggmnb-test.csv`, `manifest-ggmnb-reference.csv`). |
| `R/24-fig-decision-tree.R` | Simulation study, step 6: the diagnostic decision-tree figure with leaf numbers read from the rule tables of R/23 and R/25 (`figures/fig-diagnostic-decision-tree.pdf`). |
| `R/37-verdict-error-distances.R` | The error-distance envelope of the simulation appendix: how close to an evidence threshold every erroneous verdict sat (`tables/verdict-error-distances.csv`). |
| `R/threshold-study-helpers.R` | Shared machinery of the simulation study (truth generation, per-edge observables, checkpoint bookkeeping); sourced, not run. |
| `R/rhat-diagnostics-helpers.R` | Split-R-hat and indicator-diagnostic helpers; sourced by the study helpers, not run. |

## The simulation study

### What is archived

`results/threshold-study/` holds the complete state of the study as it was analyzed for the manuscript:

- `manifest.csv`: one row per fit (390 fits: the ordinal calibration arm, the binary and Gaussian transfer arms, and the near-boundary Gaussian follow-up arm), with the data set, the sampler setting, and the seeds of each fit. The `manifest-*.csv` files are the slices that were handed to workers.
- `datasets/`: the simulated data sets with their generating networks.
- `fits/`: one checkpoint per fit with the per-edge observables, the per-chain summaries, and the fit's timing.
- `calibration.rds`, `fit-adhd-calibration.rds`, `networks.csv`, `binary-field-gate-chosen.rds`, `inclusion-bf-correction.csv`: the calibration constants and the follow-up design.
- `logs/`: the workers' logs.

### Regenerating the appendix tables from the archive

This reproduces the manuscript's tables exactly; every step takes seconds except the binary field gate, which simulates and takes about 3 min:

```sh
Rscript R/23-threshold-study-analysis.R      # calibration arm: rule tables
Rscript R/25-threshold-study-transfer.R      # all arms: transfer table and the pooled pair rule
Rscript R/28-ggm-rule-search.R               # Gaussian-arm rule candidates
Rscript R/27-binary-field-gate.R             # the binary field gate
Rscript R/24-fig-decision-tree.R             # the decision-tree figure
Rscript R/37-verdict-error-distances.R       # the error-distance envelope
```

### Rerunning the study from scratch

Point the scripts at an empty directory first, so the archive is not overwritten:

```sh
export THRESHOLD_STUDY_DIR=results/threshold-study-rerun
```

Then, in this order:

1. `Rscript R/01-worked-example.R` (the design calibrates the ordinal arm on the Wenchuan fit) and `Rscript R/20-threshold-study-design.R`, which fits the ADHD calibration model, draws the generating networks, simulates the stage-1 and stage-2 data sets, and writes `manifest.csv`.
2. The timing pilot: `Rscript R/21-threshold-study-pilot.R --stage build` writes `manifest-pilot-test.csv` and `manifest-pilot-reference.csv`; run workers on both (step 3); then `Rscript R/21-threshold-study-pilot.R --stage collect` prints the projected cost of the grid.
3. The fits. A worker runs a slice of a manifest and checkpoints every fit:

   ```sh
   Rscript R/22-threshold-study-run.R --manifest $THRESHOLD_STUDY_DIR/manifest.csv --slot 0 --nslots 3 --cores 4
   ```

   Start `--nslots` workers with `--slot 0 ... --nslots-1`. Reference fits have 12 chains and run best as a single worker with `--cores 12`; split the manifest by its `kind` column (test versus reference) into two files to run them separately. Interrupted workers resume where they stopped, because fits with an existing checkpoint are skipped.
4. `Rscript R/26-recompute-inclusion-bf.R`, then `Rscript R/23-threshold-study-analysis.R` and `Rscript R/25-threshold-study-transfer.R`.
5. The follow-ups: `Rscript R/27-binary-field-gate.R` and `Rscript R/29-followup-design.R`, which regenerate the binary data sets with the gated coupling scale, add the near-boundary Gaussian arm, rewrite `manifest.csv`, and write `manifest-binary-rerun.csv`, `manifest-ggmnb-test.csv`, and `manifest-ggmnb-reference.csv`. Run workers on those three manifests (step 3), then repeat step 4, and finish with `R/28`, `R/24`, and `R/37`.

Environment variables: `THRESHOLD_STUDY_DIR` (checkpoint root, default `results/threshold-study`), `THRESHOLD_STUDY_OUTPUT` (where R/23 and R/24 write their tables, default `tables/`), and `THRESHOLD_STUDY_FIGURE` (the decision-tree figure path).

The full grid (390 fits) took roughly 15 CPU-hours, about a day of wall-clock time on a 15-core laptop with three workers; the pilot's projection is the tool for planning it elsewhere.

### A caveat on bit-for-bit reproduction

The archived checkpoints were produced during the development of bgms 0.2.0.0, with two development builds that straddled a change in the scale on which the package reports inclusion Bayes factors; `R/26` puts them on one scale, and the manuscript's appendix says so. A rerun with the released package is seeded fit by fit and will reproduce the design, but Markov chain Monte Carlo draws are bit-identical only for the same compiled binary, so the regenerated rule tables will agree statistically rather than digit for digit. The worked example and every other script in this repository were rerun on the released bgms 0.2.0.0 and reproduce the manuscript's numbers exactly.

## Reproducibility notes

- Every stochastic script sets a seed, typically `set.seed(123)`, and the model fits pass `seed = 123` to `bgms::bgm()`. The simulation study's design script uses its own fixed seeds, and every simulated data set and fit draws its seed from the manifest.
- `bgms::bgm()` parallelizes its chains over the available cores by default. Each chain is seeded separately, so the results do not depend on the number of cores; on Windows with RcppParallel >= 6, the package documentation notes that a seed reproduces only at the same core count.
- The bootstrap row (R/33) runs bootnet single-core on purpose, because parallel bootnet streams are not seed-reproducible.
- On 2026-10-16 every script except the multi-day fit grid was rerun from a clean clone with CRAN bgms 0.2.0.0 under R 4.6, and every CSV output was byte-identical to the file behind the manuscript; the simulation-study tables were regenerated byte-identically from the archived checkpoints.
- The bgms documentation site <https://bayesian-graphical-modelling-lab.github.io/bgms-docs/> mirrors the workflow of the tutorial.
