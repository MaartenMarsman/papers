# 18-fig1-evidence-networks.R
# Figure 1 of the manuscript: the worked-example network split into three
# panels by the evidence classification (presence / undecided / absence),
# drawn on one shared layout so the panels are directly comparable.
# The presence panel shows posterior mean interaction weights; the undecided
# and absence panels show the classified pairs at uniform width, because for
# those edges the classification, not a weight, is the result.
# The two symptoms of the motivating example (intrusion, concen) are drawn
# as filled squares, matching the introduction's EBICglasso figure; their
# edge sits in the presence panel.
#
# Reads tables/edge-classification.csv and tables/weight-summary.csv (both
# from R/01-worked-example.R) and tables/network-layout.csv (from
# R/17-fig-intro-network.R): the panels use the introduction figure's node
# positions so the two manuscript network figures are comparable node by
# node. Rerun R/01 and R/17 first if these files are absent.
#
# Output: figures/fig1-edge-evidence.pdf

set.seed(123)

edges <- utils::read.csv("tables/edge-classification.csv")
weights <- utils::read.csv("tables/weight-summary.csv")

wparts <- do.call(rbind, strsplit(weights$edge, "-", fixed = TRUE))
weights$from <- wparts[, 1]
weights$to <- wparts[, 2]
edges <- merge(edges, weights[, c("from", "to", "mean")],
               by = c("from", "to"), all.x = TRUE)
stopifnot(!anyNA(edges$mean))

node_names <- sort(unique(c(edges$from, edges$to)))
p <- length(node_names)
focal <- node_names %in% c("intrusion", "concen")

mat_for <- function(value, subset) {
  m <- matrix(0, p, p, dimnames = list(node_names, node_names))
  sel <- edges[subset, ]
  for (k in seq_len(nrow(sel))) {
    m[sel$from[k], sel$to[k]] <- value[subset][k]
    m[sel$to[k], sel$from[k]] <- value[subset][k]
  }
  m
}

w_presence <- mat_for(edges$mean, edges$evidence == "presence")
a_undecided <- mat_for(rep(1, nrow(edges)), edges$evidence == "undecided")
a_absence <- mat_for(rep(1, nrow(edges)), edges$evidence == "absence")

# One layout for all panels, shared with the introduction figure so that
# readers find every node in the same place in both manuscript figures.
layout_df <- utils::read.csv("tables/network-layout.csv")
stopifnot(setequal(layout_df$node, node_names))
layout_all <- as.matrix(
  layout_df[match(node_names, layout_df$node), c("x", "y")]
)

panel_style <- list(
  labels = node_names,
  label.cex = 1.2,
  label.scale.equal = TRUE,
  vsize = ifelse(focal, 12, 10),
  shape = ifelse(focal, "square", "circle"),
  color = ifelse(focal, "#CDE5E2", "white"),
  border.width = ifelse(focal, 3.5, 1),
  border.color = ifelse(focal, "#12655F", "grey30"),
  label.color = ifelse(focal, "#12655F", "black"),
  mar = c(4, 3, 6, 3)
)

# Triangle arrangement (presence and absence on top, undecided centered
# below, matching the manuscript's "pause on the undecided category"
# emphasis): each panel gets nearly half the page width instead of a
# third, so node labels stay legible at manuscript size.
grDevices::pdf("figures/fig1-edge-evidence.pdf", width = 9, height = 8.2)
graphics::layout(rbind(c(1, 1, 2, 2), c(0, 3, 3, 0)))

do.call(qgraph::qgraph, c(list(
  w_presence, layout = layout_all, theme = "colorblind",
  edge.width = 1.6,
  title = sprintf("evidence of presence (%d edges)",
                  sum(edges$evidence == "presence")),
  title.cex = 1.9
), panel_style))

do.call(qgraph::qgraph, c(list(
  a_absence, layout = layout_all,
  edge.color = "#A66A1E", lty = 2, edge.width = 0.9,
  title = sprintf("evidence of absence (%d edges)",
                  sum(edges$evidence == "absence")),
  title.cex = 1.9
), panel_style))

do.call(qgraph::qgraph, c(list(
  a_undecided, layout = layout_all,
  edge.color = "#8A8578", lty = 3, edge.width = 1.3,
  title = sprintf("undecided (%d edges)",
                  sum(edges$evidence == "undecided")),
  title.cex = 1.9
), panel_style))

grDevices::dev.off()

cat("panels:", sum(edges$evidence == "presence"), "presence,",
    sum(edges$evidence == "undecided"), "undecided,",
    sum(edges$evidence == "absence"), "absence\n")
