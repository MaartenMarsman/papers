# 22-run-edge-prior-check.R
# The edge-prior arm of the robustness section: prior_sensitivity_check()
# varies the interaction (slab) scale but not the prior on the structure, so
# the structure prior gets its own single refit. Settings are identical to
# R/01-worked-example.R except the edge prior: Beta-Bernoulli(1, 1), which
# treats the total number of edges as unknown instead of fixing each edge at
# fifty-fifty.
#
# Compares the resulting three-way verdicts with the reported analysis
# (tables/edge-classification.csv from R/01) and writes the cross-tab plus
# per-edge movers to tables/edge-prior-check.csv.

set.seed(123)

data("Wenchuan", package = "bgms")
wenchuan <- Wenchuan[stats::complete.cases(Wenchuan), ]

t0 <- Sys.time()
fit_bb <- bgms::bgm(
  x = wenchuan,
  variable_type = "ordinal",
  iter = 1e4,
  warmup = 5e3,
  edge_selection = TRUE,
  edge_prior = bgms::beta_bernoulli_prior(1, 1),
  interaction_prior = bgms::normal_prior(scale = 1),
  chains = 4,
  seed = 123
)
cat("fit:", round(as.numeric(Sys.time() - t0, units = "secs")), "s\n")

bf_bb <- bgms::extract_inclusion_bf(fit_bb)

classify_edges <- function(bf, threshold = 10) {
  cut(bf, breaks = c(0, 1 / threshold, threshold, Inf),
      labels = c("absence", "undecided", "presence"))
}

# Row-major pair order, as in R/01 (upper.tri would misalign edge names).
node_names <- colnames(bf_bb)
pairs <- t(utils::combn(seq_along(node_names), 2))
bb <- data.frame(
  from = node_names[pairs[, 1]],
  to = node_names[pairs[, 2]],
  inclusion_bf_bb = bf_bb[pairs],
  verdict_bb = classify_edges(bf_bb[pairs])
)

default <- utils::read.csv("tables/edge-classification.csv")
comparison <- merge(default, bb, by = c("from", "to"))
stopifnot(nrow(comparison) == nrow(default))

lv <- c("presence", "undecided", "absence")
comparison$evidence <- factor(comparison$evidence, levels = lv)
comparison$verdict_bb <- factor(comparison$verdict_bb, levels = lv)

cat("\nverdicts, Bernoulli(.5) default vs Beta-Bernoulli(1,1):\n")
print(table(default = comparison$evidence, beta_bernoulli = comparison$verdict_bb))

movers <- comparison[comparison$evidence != comparison$verdict_bb, ]
cat("\nedges changing category:", nrow(movers), "of", nrow(comparison), "\n")
cat("presence <-> absence crossings:",
    sum((movers$evidence == "presence" & movers$verdict_bb == "absence") |
          (movers$evidence == "absence" & movers$verdict_bb == "presence")), "\n")
if (nrow(movers) > 0) {
  print(movers[order(movers$evidence),
               c("from", "to", "inclusion_bf", "evidence",
                 "inclusion_bf_bb", "verdict_bb")])
}

utils::write.csv(comparison, "tables/edge-prior-check.csv", row.names = FALSE)
