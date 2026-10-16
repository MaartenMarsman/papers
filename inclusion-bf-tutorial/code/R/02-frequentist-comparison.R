# 02-frequentist-comparison.R
# The frequentist analysis from the introduction (EBICglasso on the Wenchuan
# PTSD data) and its confrontation with the Bayesian three-way edge
# classification from 01-worked-example.R.
#
# NOTE ON PROVENANCE: this reproduces the MODERN standard workflow (regularized
# EBICglasso, Epskamp & Fried 2018) applied to McNally et al.'s (2015) data.
# It is NOT McNally's own method: the 2015 paper estimated association
# (zero-order, r >= .30) and concentration (UNREGULARIZED partial correlation)
# networks plus a relative-importance network -- no graphical lasso.
#
# Produces the cross-tabulation (tables/glasso-vs-bayes.csv,.tex) and the
# stacked-bar comparison figure (figures/fig-glasso-vs-bayes.pdf).
#
# Requires tables/edge-classification.csv from R/01-worked-example.R.

set.seed(123)

# --- Data (identical preprocessing to 01-worked-example.R) ---------------------

data("Wenchuan", package = "bgms")
wenchuan <- Wenchuan[stats::complete.cases(Wenchuan), ]

# --- Frequentist workflow: polychoric correlations + EBICglasso ----------------

correlations <- qgraph::cor_auto(wenchuan)
glasso_network <- qgraph::EBICglasso(
  S = correlations,
  n = nrow(wenchuan),
  gamma = 0.5
)

node_names <- colnames(glasso_network)
pairs <- t(utils::combn(seq_along(node_names), 2))

frequentist <- data.frame(
  from = node_names[pairs[, 1]],
  to = node_names[pairs[, 2]],
  glasso_weight = glasso_network[pairs],
  glasso_selected = glasso_network[pairs] != 0
)

cat("EBICglasso selects", sum(frequentist$glasso_selected), "of",
    nrow(frequentist), "edges\n")

# --- Join with the Bayesian three-way classification ---------------------------

bayes <- utils::read.csv("tables/edge-classification.csv")
bayes$evidence <- factor(bayes$evidence,
                         levels = c("presence", "undecided", "absence"))

comparison <- merge(frequentist, bayes, by = c("from", "to"))
stopifnot(nrow(comparison) == nrow(frequentist))

crosstab <- table(
  glasso = ifelse(comparison$glasso_selected, "selected", "zero"),
  evidence = comparison$evidence
)
cat("\nEBICglasso decision versus Bayesian evidence:\n")
print(crosstab)

# The two headline numbers for the manuscript:
zero_edges <- comparison[!comparison$glasso_selected, ]
selected_edges <- comparison[comparison$glasso_selected, ]
cat("\nOf the", nrow(zero_edges), "edges glasso sets to zero:",
    sum(zero_edges$evidence == "absence"), "have evidence of absence,",
    sum(zero_edges$evidence == "undecided"), "are undecided\n")
cat("Of the", nrow(selected_edges), "edges glasso selects:",
    sum(selected_edges$evidence == "presence"), "have evidence of presence,",
    sum(selected_edges$evidence == "undecided"), "are undecided\n")

utils::write.csv(comparison, file = "tables/glasso-vs-bayes.csv",
                 row.names = FALSE)

# --- Comparison figure: Bayesian verdicts within each glasso decision ----------

evidence_colors <- c(
  presence = "#12655F",
  undecided = "#B9B4A8",
  absence = "#A66A1E"
)

plot_data <- as.data.frame(crosstab)
plot_data$glasso <- factor(
  plot_data$glasso,
  levels = c("selected", "zero"),
  labels = c("Edge selected by EBICglasso", "Edge set to zero by EBICglasso")
)

fig_comparison <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = glasso, y = Freq, fill = evidence)
) +
  ggplot2::geom_col(position = "stack", width = 0.6) +
  ggplot2::scale_fill_manual(values = evidence_colors, name = "Bayesian evidence") +
  ggplot2::labs(x = NULL, y = "Number of edges") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())

ggplot2::ggsave(
  filename = "figures/fig-glasso-vs-bayes.pdf",
  plot = fig_comparison,
  width = 6.5,
  height = 4
)

# --- Narrative anchors for the introduction ------------------------------------

# Named edges for the running example: strongest evidence of presence, the
# strongest glasso edge that is nevertheless undecided, and glasso zeros with
# genuine evidence of absence.
cat("\nTop 5 edges by inclusion Bayes factor:\n")
print(utils::head(
  comparison[order(comparison$inclusion_bf, decreasing = TRUE),
             c("from", "to", "glasso_weight", "inclusion_bf", "evidence")], 5))

undecided_selected <- selected_edges[selected_edges$evidence == "undecided", ]
cat("\nGlasso-selected but undecided (candidates for the intro example):\n")
print(undecided_selected[order(abs(undecided_selected$glasso_weight),
                               decreasing = TRUE),
                         c("from", "to", "glasso_weight", "inclusion_bf")])

# --- Independent audit of the confrontation --------------------------------

# The extreme disagreements in the cross-tab invite the suspicion that the
# comparison is broken. These checks involve neither bgms nor the lasso
# penalty: unregularized partial correlations (polychoric and raw Pearson),
# an ordinal regression on the raw categories, and a milder lasso tuning.
P_poly <- -stats::cov2cor(solve(correlations))
P_raw <- -stats::cov2cor(solve(stats::cor(wenchuan)))
P_spear <- -stats::cov2cor(solve(stats::cor(wenchuan, method = "spearman")))
n <- nrow(wenchuan)
k <- ncol(wenchuan) - 2
t_partial <- function(r) r * sqrt((n - 2 - k) / (1 - r^2))
cat("\nUnregularized partial correlations (audit):\n")
for (pair in list(c("intrusion", "concen"), c("intrusion", "sleep"))) {
  cat(sprintf(
    "%s-%s: polychoric %.3f (t = %.2f), raw Pearson %.3f (t = %.2f), Spearman %.3f\n",
    pair[1], pair[2],
    P_poly[pair[1], pair[2]], t_partial(P_poly[pair[1], pair[2]]),
    P_raw[pair[1], pair[2]], t_partial(P_raw[pair[1], pair[2]]),
    P_spear[pair[1], pair[2]]))
}
lm_fit <- stats::lm(concen ~ ., data = as.data.frame(wenchuan))
cat("Linear regression concen ~ rest, intrusion coefficient:\n")
print(round(summary(lm_fit)$coefficients["intrusion", ], 4))
polr_fit <- MASS::polr(factor(concen, ordered = TRUE) ~ .,
                       data = as.data.frame(wenchuan), Hess = TRUE)
cat("Ordinal regression concen ~ rest, intrusion coefficient:\n")
print(round(summary(polr_fit)$coefficients["intrusion", ], 3))
glasso_g0 <- qgraph::EBICglasso(S = correlations, n = n, gamma = 0)
cat("intrusion-concen weight under gamma = 0:",
    round(glasso_g0["intrusion", "concen"], 4),
    "(selected:", glasso_g0["intrusion", "concen"] != 0, ")\n")
cat("intrusion-sleep weight under gamma = 0.5:",
    round(glasso_network["intrusion", "sleep"], 4),
    "(selected:", glasso_network["intrusion", "sleep"] != 0, ")\n")

# Why the selection flips evidences (numbers quoted in the manuscript's
# confrontation subsection): the selection is systematically harsher on
# negative partials (suppressor-type edges whose contribution the positive
# manifold can absorb), and the selected-but-absence cells are the lasso's
# smallest surviving weights.
ut <- upper.tri(P_poly)
d_audit <- data.frame(absP = abs(P_poly[ut]),
                      sel = glasso_network[ut] != 0,
                      sgn = sign(P_poly[ut]))
cat("\nSign asymmetry of the selection (unregularized polychoric partials):\n")
print(table(sign = ifelse(d_audit$sgn < 0, "negative", "positive"),
            glasso = ifelse(d_audit$sel, "selected", "zero")))
cat("marginal polychoric correlation intrusion-concen:",
    round(correlations["intrusion", "concen"], 2),
    paste0("| negative marginal correlations among ", nrow(frequentist), ":"),
    sum(correlations[upper.tri(correlations)] < 0), "\n")
cat("mean |partial| among zeroed edges, negative vs positive:",
    round(mean(d_audit$absP[!d_audit$sel & d_audit$sgn < 0]), 3), "vs",
    round(mean(d_audit$absP[!d_audit$sel & d_audit$sgn > 0]), 3), "\n")
r_focal <- abs(P_poly["intrusion", "concen"])
cat(paste0("rank of intrusion-concen |partial| among ", nrow(d_audit),
           " (1 = largest):"),
    sum(d_audit$absP >= r_focal), "\n")
cat("selected edges with smaller |partial| than intrusion-concen:",
    sum(d_audit$sel & d_audit$absP < r_focal), "of", sum(d_audit$sel), "\n")
cat("median |glasso weight|, selected-but-absence edges:",
    round(stats::median(abs(selected_edges$glasso_weight[
      selected_edges$evidence == "absence"])), 3),
    "| selected-and-presence:",
    round(stats::median(abs(selected_edges$glasso_weight[
      selected_edges$evidence == "presence"])), 3), "\n")
