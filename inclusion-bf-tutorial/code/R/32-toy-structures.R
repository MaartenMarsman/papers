# 32-toy-structures.R
# Toy tier for Section 2: the Bayesian machinery on the smallest example with
# real structure uncertainty. Three Wenchuan symptoms (intrusive thoughts,
# sleep problems, concentration problems) admit 2^3 = 8 candidate networks.
# The script fits the same model as the worked example to just these three
# symptoms, tabulates the posterior probability of every candidate structure,
# and computes the inclusion Bayes factor of one edge by the literal summation
# of the manuscript's Equation 1, so the text can show the sum being evaluated.
#
# The trio deliberately contains both focal edges of the paper (intrusion-
# concen, intrusion-sleep). The three-variable verdicts need not match the
# seventeen-variable ones, and that difference is itself the lesson the
# manuscript draws: an edge is conditional on what else is in the network.
#
# Requires bgms >= 0.2.0.0. Runtime: a few minutes (four million draws).
# Output: tables/toy-structures.csv and printed numbers. The fitted object
# itself is not saved (see the note after the fit).

set.seed(123)

data("Wenchuan", package = "bgms")
wenchuan <- Wenchuan[stats::complete.cases(Wenchuan), ]
trio <- c("intrusion", "sleep", "concen")

# Four million posterior draws: with 1e6 sampling iterations
# per chain the structure probabilities are stable to three decimals and the
# never-visited structures are convincingly zero.
fit3 <- bgms::bgm(
  x = wenchuan[, trio],
  iter = 1e6,
  warmup = 5e3,
  chains = 4,
  seed = 123
)
# The full fitted object is NOT saved: at four million draws it is ~1 GB, and
# every downstream consumer (R/35, the manuscript numbers) reads the small
# tables/toy-structures.csv written below. Rerun this script to regenerate.

# --- Posterior probability of every candidate structure ------------------------

indicators <- bgms::extract_indicators(fit3)
cat("indicator columns:", colnames(indicators), "\n")

# Fixed edge order for labels and the csv: within-trio pairs.
edge_names <- colnames(indicators)
n_edges <- length(edge_names)
stopifnot(n_edges == 3)

# One row per draw -> structure id "101" etc. (1 = edge present).
pattern <- apply(indicators, 1, paste0, collapse = "")
all_patterns <- apply(expand.grid(rep(list(0:1), n_edges))[, n_edges:1], 1,
                      paste0, collapse = "")
posterior_prob <- vapply(
  all_patterns,
  function(p) mean(pattern == p),
  numeric(1)
)

structures <- data.frame(
  structure = all_patterns,
  prior_prob = rep(1 / 2^n_edges, 2^n_edges),
  posterior_prob = posterior_prob
)
for (j in seq_len(n_edges)) {
  structures[[edge_names[j]]] <- as.integer(substr(structures$structure, j, j))
}
structures <- structures[order(-structures$posterior_prob), ]
print(structures, row.names = FALSE, digits = 3)

utils::write.csv(structures, "tables/toy-structures.csv", row.names = FALSE)

# --- The inclusion Bayes factor by literal summation ---------------------------

cat("\nEquation-1 summation per edge (posterior odds / prior odds):\n")
for (j in seq_len(n_edges)) {
  present <- structures[[edge_names[j]]] == 1
  post_incl <- sum(structures$posterior_prob[present])
  prior_incl <- sum(structures$prior_prob[present])
  bf_sum <- (post_incl / (1 - post_incl)) / (prior_incl / (1 - prior_incl))
  cat(sprintf("  %-20s P(incl | data) = %.3f   BF by summation = %.2f\n",
              edge_names[j], post_incl, bf_sum))
}

cat("\nextract_inclusion_bf (Rao-Blackwellized, for comparison):\n")
print(bgms::extract_inclusion_bf(fit3), digits = 3)

cat("\nPosterior inclusion probabilities (Rao-Blackwellized):\n")
print(bgms::extract_posterior_inclusion_probabilities(fit3), digits = 3)

cat("\nPosterior mean weights (model-averaged):\n")
print(colMeans(bgms::extract_pairwise_interactions(fit3)), digits = 3)

# --- Full-network comparison for the manuscript's closing sentence -------------

fit17 <- readRDS("tables/fit-wenchuan.rds")
bf17 <- bgms::extract_inclusion_bf(fit17)
cat("\nSame pairs in the seventeen-variable network:\n")
cat(sprintf("  intrusion-sleep  BF = %.3f\n", bf17["intrusion", "sleep"]))
cat(sprintf("  intrusion-concen BF = %.2f\n", bf17["intrusion", "concen"]))
cat(sprintf("  sleep-concen     BF = %.3f\n", bf17["sleep", "concen"]))
