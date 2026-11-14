# Table for Section 3.6: realized diagonal mean and realized variance of the
# included reference entry k_12 relative to their specified values (both one),
# from R/results/section36_marginals_bgms.rds (section36_marginals_bgms.R).
# Standard errors are printed to the console and belong in the caption.
#
# Run:    R_PROFILE_USER=/dev/null Rscript R/scripts/section36_marginals_table.R
# Output: tables/section36_marginals.tex

PROJECT_ROOT <- Sys.getenv("SPIKESLAB_ROOT", unset = getwd())
res <- readRDS(file.path(PROJECT_ROOT, "R", "results", "section36_marginals_bgms.rds"))
s <- res$summary
deltas <- sort(unique(s$delta)); qs <- sort(unique(s$q))
cell <- function(q, d, col) s[s$q == q & s$delta == d, col]

lines <- c(
  sprintf("%% Source: R/results/section36_marginals_bgms.rds (bgms %s, %d graphs x %d chains x %d draws per cell).",
          res$bgms_version, res$design$n_graphs, res$design$n_chains, res$design$n_samples),
  sprintf("%% Standard errors over reference graphs: diagonal mean at most %.3f, variance at most %.4f.",
          max(s$diag_ratio_se), max(s$var_ratio_se)),
  sprintf("\\begin{tabular}{l%s}", paste(rep("rr", length(deltas)), collapse = "")),
  "\\toprule",
  paste0(" & ", paste(sprintf("\\multicolumn{2}{c}{$\\delta=%g$}", deltas), collapse = " & "), " \\\\"),
  paste(sprintf("\\cmidrule(lr){%d-%d}", 2 + 2 * (seq_along(deltas) - 1), 3 + 2 * (seq_along(deltas) - 1)), collapse = ""),
  paste0("$q$ & ", paste(rep("Mean $k_{ll}$ & Var $k_{12}$", length(deltas)), collapse = " & "), " \\\\"),
  "\\midrule")
for (q in qs) {
  vals <- unlist(lapply(deltas, function(d) c(sprintf("%.2f", cell(q, d, "diag_ratio")),
                                                sprintf("%.2f", cell(q, d, "var_ratio")))))
  lines <- c(lines, paste0(q, " & ", paste(vals, collapse = " & "), " \\\\"))
}
lines <- c(lines, "\\midrule",
           paste0("Specified & ", paste(rep("1 & 1", length(deltas)), collapse = " & "), " \\\\"),
           paste0("Lower bound $(\\delta+1)/\\beta$ & ",
                  paste(sprintf("%g & ", deltas + 1), collapse = " & "), " \\\\"),
           "\\bottomrule", "\\end{tabular}")
writeLines(lines, file.path(PROJECT_ROOT, "tables", "section36_marginals.tex"))
cat(paste(lines, collapse = "\n"), "\n")
cat(sprintf("\nSE ranges: diagonal %.4f to %.4f; variance %.4f to %.4f\n",
            min(s$diag_ratio_se), max(s$diag_ratio_se), min(s$var_ratio_se), max(s$var_ratio_se)))
