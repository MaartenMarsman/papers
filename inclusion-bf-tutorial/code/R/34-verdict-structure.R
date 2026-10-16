# 34-verdict-structure.R
# Reproduces every number in the manuscript's substantive reading of the
# verdict panels (Section 3.4): the DSM-IV criterion structure of the absence
# panel, the concentration of undecided verdicts inside the avoidance and
# numbing cluster, the strongest individual absences, and the strongest
# boundary-crossing (bridge-candidate) edges.
#
# Reads tables/edge-supplement.csv (from R/01-worked-example.R).
# The criterion grouping (B = intrusion, C = avoidance and numbing,
# D = hyperarousal) follows the DSM-IV criterion clusters as labeled in the
# data documentation.

edges <- utils::read.csv("tables/edge-supplement.csv")

criterion <- c(
  intrusion = "B", dreams = "B", flash = "B", upset = "B", physior = "B",
  avoidth = "C", avoidact = "C", amnesia = "C", lossint = "C",
  distant = "C", numb = "C", future = "C",
  sleep = "D", anger = "D", concen = "D", hyper = "D", startle = "D"
)

parts <- strsplit(edges$edge, "-", fixed = TRUE)
edges$from <- vapply(parts, `[`, character(1), 1)
edges$to <- vapply(parts, `[`, character(1), 2)
edges$pair <- paste0(
  pmin(criterion[edges$from], criterion[edges$to]),
  pmax(criterion[edges$from], criterion[edges$to])
)

cat("criterion-pair by verdict:\n")
print(table(edges$pair, edges$verdict))

within <- edges$pair %in% c("BB", "CC", "DD")
cat(sprintf(
  "\nabsences crossing criteria: %d of %d\n",
  sum(!within & edges$verdict == "absence"),
  sum(edges$verdict == "absence")
))

cat("\nstrongest absences (largest BF against inclusion):\n")
absence <- edges[edges$verdict == "absence", ]
absence <- absence[order(absence$log_bf), ]
absence$bf_01 <- exp(-absence$log_bf)
print(utils::head(absence[, c("edge", "log_bf", "bf_01")], 5),
      row.names = FALSE, digits = 3)

cat("\nboundary-crossing presence edges (bridge candidates):\n")
bridges <- edges[!within & edges$verdict == "presence", ]
bridges <- bridges[order(-bridges$log_bf), ]
print(utils::head(bridges[, c("edge", "pair", "log_bf", "weight_mean")], 5),
      row.names = FALSE, digits = 3)
