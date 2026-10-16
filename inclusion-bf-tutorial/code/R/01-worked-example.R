# 01-worked-example.R
# Worked example for the tutorial (manuscript Section: Worked example in R).
# Fits a Bayesian ordinal Markov random field to the Wenchuan PTSD data,
# classifies every edge as presence / absence / undecided, writes the edge
# tables, and produces the strength-versus-evidence figure.
#
# Requires bgms >= 0.2.0.0. Uses the tutorial's default prior distributions
# and a longer run than the default sampler settings (see the fit block);
# set smoke_test <- TRUE for a fast syntax check (NOT valid for inference).

smoke_test <- FALSE

set.seed(123)

# --- Data ---------------------------------------------------------------------

data("Wenchuan", package = "bgms")
wenchuan <- Wenchuan[stats::complete.cases(Wenchuan), ]

# --- Fit (prior settings stated explicitly for the reader) --------------------

# Longer than the Table-1 defaults (2,000 + 2,000): 10,000 iterations after
# 5,000 warmup in four chains, so that every parameter reaches an effective
# sample size in the hundreds on these data.
n_iter <- if (smoke_test) 100 else 1e4
n_warmup <- if (smoke_test) 100 else 5e3
n_chains <- if (smoke_test) 2 else 4

fit <- bgms::bgm(
  x = wenchuan,
  variable_type = "ordinal",
  iter = n_iter,
  warmup = n_warmup,
  edge_selection = TRUE,
  edge_prior = bgms::bernoulli_prior(0.5),
  interaction_prior = bgms::normal_prior(scale = 1),
  chains = n_chains,
  seed = 123
)

# Save the fitted object; every downstream script reads it from here.
saveRDS(fit, "tables/fit-wenchuan.rds")

# --- Convergence checks (report before any substantive result) -----------------

# One machinery for everything since bgms's Rao-Blackwellized summary: the
# edge indicators are summarized by their continuous RB draws, so split-Rhat,
# ESS, and MCSE apply to them like to any other parameter. bgms 0.2.0.0 drops
# the transition-based n_eff_mixt column; the inclusion table now carries the
# RB triple (mcse, n_eff, Rhat) beside the directional flip counts n0->1 and
# n1->0, which stay as the exploration record. Read the RB precision against
# those counts; a large RB ESS with zero flips is the boundary signature (a
# precise one-step estimate resting on no direct transition evidence), and
# cross-chain agreement decides. See the manuscript's convergence subsection.
rhat <- bgms::extract_rhat(fit)
ess <- bgms::extract_ess(fit)
diagnostics <- data.frame(
  component = names(rhat),
  max_rhat = vapply(rhat, function(v) max(unlist(v), na.rm = TRUE), 0),
  min_ess = vapply(ess[names(rhat)], function(v) min(unlist(v), na.rm = TRUE), 0)
)
print(diagnostics, row.names = FALSE)
utils::write.csv(diagnostics, file = "tables/diagnostics.csv", row.names = FALSE)

# Full indicator and weight summaries for the manuscript's numbers: the
# RB/transition ESS pair with directional flip counts, and the derived
# composite ESS with its inclusion share for the model-averaged weights.
fit_summary <- summary(fit)
utils::write.csv(cbind(edge = rownames(fit_summary$indicator),
                       as.data.frame(fit_summary$indicator)),
                 file = "tables/indicator-summary.csv", row.names = FALSE)
utils::write.csv(cbind(edge = rownames(fit_summary$pairwise),
                       as.data.frame(fit_summary$pairwise)),
                 file = "tables/weight-summary.csv", row.names = FALSE)

# --- Three-way edge classification ---------------------------------------------

# Inclusion Bayes factors straight from bgms, computed from the
# Rao-Blackwellized odds accumulators (prior inclusion odds divided out),
# which stay finite even for edges whose raw inclusion frequency saturates
# at 1 in a finite run. An edge can still be +Inf when its exclusion
# evidence falls below double precision entirely, and it is reported as
# such. The log argument (log = FALSE by default) returns the natural-log
# Bayes factor instead.
posterior_incl <- bgms::extract_posterior_inclusion_probabilities(fit)
bf_matrix <- bgms::extract_inclusion_bf(fit)
log_bf_matrix <- bgms::extract_inclusion_bf(fit, log = TRUE)

# Scale checks: the two extractors must agree (exp of the log Bayes factor
# equals the Bayes factor wherever it is finite), and any analysis with
# evidence of absence must contain Bayes factors below 1.
ok <- is.finite(log_bf_matrix)
stopifnot(isTRUE(all.equal(exp(log_bf_matrix[ok]), bf_matrix[ok])))
stopifnot(any(bf_matrix < 1))

# Enumerate edges in row-major pair order (the order bgms names edges in);
# upper.tri() would give column-major order and silently misalign edges.
node_names <- colnames(posterior_incl)
pairs <- t(utils::combn(seq_along(node_names), 2))
edge_names <- paste(node_names[pairs[, 1]], node_names[pairs[, 2]], sep = "-")

posterior_p <- posterior_incl[pairs]

# Natural log throughout, the unit bgms reports in (verdicts(), the sensitivity
# curve, and plot_edge_posterior() all speak it), so the evidence thresholds of
# 10 and 1/10 sit at +/- log(10) = 2.303.
log_bf <- log_bf_matrix[pairs]

# The cut behind bgms::verdicts(): thresholds on the Bayes factor scale,
# with the one infinite Bayes factor falling into the presence category via
# the Inf break. Kept explicit here so downstream tables do not depend on
# the verdicts() print format.
classify_edges <- function(bf, threshold = 10) {
  cut(bf, breaks = c(0, 1 / threshold, threshold, Inf),
      labels = c("absence", "undecided", "presence"))
}

edges <- data.frame(
  from = node_names[pairs[, 1]],
  to = node_names[pairs[, 2]],
  posterior_inclusion = posterior_p,
  inclusion_bf = bf_matrix[pairs],
  log_bf = log_bf,
  evidence = classify_edges(bf_matrix[pairs])
)
edges <- edges[order(edges$log_bf, decreasing = TRUE), ]

cat("\nEdge evidence classification:\n")
print(table(edges$evidence))

utils::write.csv(edges, file = "tables/edge-classification.csv", row.names = FALSE)

# --- Per-edge supplement table --------------------------------------------------

# The full table the reporting section's last item points at: one row per
# candidate edge with the inclusion Bayes factor, the verdict and its
# Monte-Carlo fragility flag (verdicts(), union of the two-state and
# Rao-Blackwellized standard-error channels), and the model-averaged weight
# with its 95% credible interval. This is the file to ship as the online
# supplement.
verdict_report <- bgms::verdicts(fit)
weight_draws <- bgms::extract_pairwise_interactions(fit)
weight_draws <- weight_draws[, match(verdict_report$parameter, colnames(weight_draws))]
supplement <- data.frame(
  edge = verdict_report$parameter,
  posterior_inclusion = verdict_report$pip,
  inclusion_bf = verdict_report$bf,
  log_bf = verdict_report$log_bf,
  verdict = verdict_report$verdict,
  fragile = verdict_report$fragile,
  weight_mean = colMeans(weight_draws),
  weight_lower = apply(weight_draws, 2, stats::quantile, 0.025),
  weight_upper = apply(weight_draws, 2, stats::quantile, 0.975)
)
supplement <- supplement[order(supplement$log_bf, decreasing = TRUE), ]
utils::write.csv(supplement, file = "tables/edge-supplement.csv", row.names = FALSE)
cat("\nwrote tables/edge-supplement.csv (", nrow(supplement), "edges )\n")

# --- Figure 1 -------------------------------------------------------------------

# Figure 1 (figures/fig1-edge-evidence.pdf) is produced by
# R/18-fig1-evidence-networks.R from tables/edge-classification.csv and
# tables/weight-summary.csv written here: the network split into three panels
# by evidence category on one shared layout. Rerun that script after this one.

# --- Figure 3: association strength versus evidence (quadrant view) ------------

# Evidence is not magnitude: model-averaged interaction strength against the
# inclusion Bayes factor, so strong-but-undecided and weak-but-certain edges
# are visible as distinct quadrants.

evidence_colors <- c(
  presence = "#12655F",
  undecided = "#B9B4A8",
  absence = "#A66A1E"
)

# Display cap: the largest log Bayes factors lie so far above the thresholds
# that an axis wide enough to show them hides the thresholds. Edges beyond a log
# Bayes factor of 25, including any +Inf edge, are drawn AT the cap with an
# open triangle; the caption reports the convention and the table carries the
# exact values.
bf_cap <- 25

# extract_pairwise_interactions() returns posterior samples (draws x edges);
# column means are the model-averaged strengths. Align columns by edge name.
interaction_samples <- bgms::extract_pairwise_interactions(fit)
strength <- colMeans(interaction_samples)
name_match <- match(edge_names, colnames(interaction_samples))
stopifnot(!anyNA(name_match))
strength <- strength[name_match]

edges_by_index <- data.frame(
  strength = strength,
  log_bf = log_bf,
  evidence = classify_edges(bf_matrix[pairs])
)
edges_by_index$bf_display <- pmin(edges_by_index$log_bf, bf_cap)
edges_by_index$capped <- edges_by_index$log_bf > bf_cap

fig3 <- ggplot2::ggplot(
  edges_by_index,
  ggplot2::aes(x = strength, y = bf_display, color = evidence,
               shape = capped)
) +
  ggplot2::geom_hline(yintercept = c(-log(10), log(10)), linetype = "dashed",
                      color = "grey50") +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_color_manual(values = evidence_colors, name = "Evidence") +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 2), guide = "none") +
  ggplot2::labs(
    x = "Model-averaged interaction strength",
    y = "Log inclusion Bayes factor"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

ggplot2::ggsave(
  filename = "figures/fig3-strength-vs-evidence.pdf",
  plot = fig3,
  width = 6.5,
  height = 4.5
)
