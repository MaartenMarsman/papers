# 17-fig-intro-network.R
# Introduction figure: the Wenchuan PTSD network estimated with the field's
# standard workflow (polychoric correlations + EBICglasso, as in
# R/02-frequentist-comparison.R), drawn with qgraph. The two nodes of the
# motivating example (intrusion, concen) are emphasized; the point of the
# figure is the edge that is NOT there between them.
#
# Output: figures/fig-intro-network.pdf

set.seed(123)

data("Wenchuan", package = "bgms")
wenchuan <- Wenchuan[stats::complete.cases(Wenchuan), ]

correlations <- qgraph::cor_auto(wenchuan)
glasso_network <- qgraph::EBICglasso(
  S = correlations,
  n = nrow(wenchuan),
  gamma = 0.5
)

node_names <- colnames(glasso_network)
focal <- node_names %in% c("intrusion", "concen")

cat("edges selected:", sum(glasso_network[upper.tri(glasso_network)] != 0),
    "of", sum(upper.tri(glasso_network)), "\n")
cat("intrusion-concen weight:",
    glasso_network["intrusion", "concen"], "\n")

# Compute the spring layout once and save it: the evidence figure
# (R/18-fig1-evidence-networks.R) reuses these node positions so that both
# manuscript network figures are directly comparable node by node.
layout_spring <- qgraph::qgraph(glasso_network, layout = "spring",
                                DoNotPlot = TRUE)$layout
utils::write.csv(
  data.frame(node = node_names, x = layout_spring[, 1], y = layout_spring[, 2]),
  file = "tables/network-layout.csv", row.names = FALSE
)

grDevices::pdf("figures/fig-intro-network.pdf", width = 6.5, height = 5)
qgraph::qgraph(
  glasso_network,
  layout = layout_spring,
  theme = "colorblind",
  labels = node_names,
  label.cex = 1.1,
  label.scale.equal = TRUE,
  vsize = ifelse(focal, 11, 9),
  shape = ifelse(focal, "square", "circle"),
  color = ifelse(focal, "#CDE5E2", "white"),
  border.width = ifelse(focal, 4, 1),
  border.color = ifelse(focal, "#12655F", "grey30"),
  label.color = ifelse(focal, "#12655F", "black"),
  edge.width = 1.4,
  mar = c(3, 3, 3, 3)
)
grDevices::dev.off()
