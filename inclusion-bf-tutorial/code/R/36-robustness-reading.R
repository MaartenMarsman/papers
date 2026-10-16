# 36-robustness-reading.R
# Reproduces every number in the manuscript's prior-sensitivity reading of
# the robustness checks: the split of the edges into those robust to the
# interaction scale and the genuine movers, the direction of the movers, and
# the substantive reading of which named findings survive the scale check.
#
# Reads tables/prior-sensitivity-edges.csv (from R/16-run-sensitivity-check.R).

ps <- utils::read.csv("tables/prior-sensitivity-edges.csv")

# A genuine mover changed beyond run-to-run noise AND was precise enough at
# the chosen scale for that comparison; the insufficient edges are excluded
# (they are the printed report's "not certifiable" category, a Monte Carlo
# matter, not a prior matter). The cross-table below counts the four classes.
ps$genuine <- ps$mover == "moved-beyond-wobble" & !ps$insufficient

cat("sensitivity classes (mover label x insufficient):\n")
print(table(ps$mover, ps$insufficient))
cat(sprintf("\nrobust to the scale (all but the genuine movers): %d of %d\n",
            sum(!ps$genuine), nrow(ps)))
cat(sprintf("genuine movers: %d\n", sum(ps$genuine)))

# No mover crosses between presence and absence.
movers <- ps[ps$genuine, ]
anchor_cols <- c("verdict_x0.4", "verdict_x0.63", "verdict_x1",
                 "verdict_x1.6", "verdict_x2.5")
crossings <- apply(movers[, anchor_cols], 1, function(v)
  any(v == "presence") && any(v == "absence"))
cat(sprintf("movers whose anchors include both presence and absence: %d\n",
            sum(crossings)))

# Named findings of the panels reading.
named <- c("intrusion-concen", "intrusion-future", "lossint-concen",
           "future-concen", "dreams-future", "flash-concen",
           "distant-sleep", "amnesia-startle", "intrusion-sleep")
cat("\nnamed findings:\n")
sel <- ps[ps$edge %in% named, c("edge", "chosen_scale_verdict", "genuine",
                                "insufficient", "stability_lower",
                                "stability_upper")]
print(sel[match(named, sel$edge), ], row.names = FALSE, digits = 3)

# Boundary absences (avoidance/numbing vs the other two criteria).
criterion <- c(
  intrusion = "B", dreams = "B", flash = "B", upset = "B", physior = "B",
  avoidth = "C", avoidact = "C", amnesia = "C", lossint = "C",
  distant = "C", numb = "C", future = "C",
  sleep = "D", anger = "D", concen = "D", hyper = "D", startle = "D"
)
parts <- strsplit(ps$edge, "-", fixed = TRUE)
p1 <- criterion[vapply(parts, `[`, character(1), 1)]
p2 <- criterion[vapply(parts, `[`, character(1), 2)]
pair <- paste0(pmin(p1, p2), pmax(p1, p2))
ba <- ps$chosen_scale_verdict == "absence" & pair %in% c("BC", "CD")
cat(sprintf("\nboundary absences (B-C, C-D): %d\n", sum(ba)))
cat(sprintf("  hold across the range (stable or within noise): %d\n",
            sum(ba & !ps$genuine & !ps$insufficient)))
cat(sprintf("  genuinely scale-dependent: %d\n", sum(ba & ps$genuine)))
cat(sprintf("  too noisy to assess: %d\n", sum(ba & ps$insufficient)))
