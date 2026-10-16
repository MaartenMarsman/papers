# 33-bootstrap-vs-posterior.R
# Executes one row of the translation table side by side, for the manuscript's
# practical section: the bootnet edge-weight bootstrap of the estimation
# workflow against the posterior credible intervals that the Bayesian fit
# carries natively. Reports both intervals for the paper's two focal edges
# (intrusion-concen, intrusion-sleep) and the wall-clock cost of each route.
#
# Requires bootnet (and the fitted object from R/01-worked-example.R).
# Output: tables/bootstrap-edge-cis.csv, printed numbers for the manuscript.

set.seed(123)

data("Wenchuan", package = "bgms")
wenchuan <- Wenchuan[stats::complete.cases(Wenchuan), ]

t_boot <- system.time({
  net <- bootnet::estimateNetwork(wenchuan, default = "EBICglasso")
  boots <- bootnet::bootnet(net, nBoots = 1000, nCores = 1,  # one core: parallel streams are not seed-reproducible
                            statistics = "edge")
})

boot_tab <- summary(boots, statistics = "edge")
boot_tab <- as.data.frame(boot_tab)
focal <- boot_tab[
  (boot_tab$node1 == "intrusion" & boot_tab$node2 %in% c("concen", "sleep")) |
  (boot_tab$node2 == "intrusion" & boot_tab$node1 %in% c("concen", "sleep")), ]
cols <- intersect(
  c("node1", "node2", "sample", "mean", "CIlower", "CIupper", "q2.5", "q97.5"),
  colnames(focal)
)
focal <- focal[, cols]
print(focal, row.names = FALSE, digits = 3)
cat(sprintf("\nbootstrap wall clock: %.0f s (1000 bootstraps, 1 core)\n",
            t_boot["elapsed"]))

utils::write.csv(focal, "tables/bootstrap-edge-cis.csv", row.names = FALSE)

# Posterior credible intervals for the same pairs, already in the fit.
fit <- readRDS("tables/fit-wenchuan.rds")
draws <- bgms::extract_pairwise_interactions(fit)
for (pair in c("intrusion-concen", "intrusion-sleep")) {
  ci <- stats::quantile(draws[, pair], probs = c(.025, .975))
  cat(sprintf("posterior %s: mean %.3f, 95%% CI [%.3f, %.3f]\n",
              pair, mean(draws[, pair]), ci[1], ci[2]))
}
