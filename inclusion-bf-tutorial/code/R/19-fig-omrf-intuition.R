# 19-fig-omrf-intuition.R
# Didactic illustration embedded in Box 1 (the ordinal Markov random field):
# a two-node network (focal item j, neighbor k, positive edge weight omega)
# followed by one panel per value of the neighbor's answer, each showing the
# model-implied distribution of the focal item's answers.
# The leftmost bar panel (neighbor at 0) is what the category intercepts alone
# encode; moving right, the edge weight shifts the distribution upward one
# step at a time. Parameter values are illustrative, not estimated.
#
# Output: figures/fig-omrf-intuition.pdf

set.seed(123)

# Category intercepts mu_c (baseline log preference of category c relative to
# category 0) chosen to give a right-skewed baseline typical of symptom data,
# and a positive edge weight omega on bgms's unified partial association
# scale: raising the neighbor one category shifts the adjacent-category log
# odds by 2 * omega.
mu <- c(0, -0.4, -1.1, -1.9, -2.7)
omega <- 0.2
categories <- 0:4
neighbor_values <- 0:4

# P(x_j = c | neighbor = s) proportional to exp(mu_c + 2 * omega * c * s)
cell_list <- vector("list", length(neighbor_values))
for (i in seq_along(neighbor_values)) {
  s <- neighbor_values[i]
  weights <- exp(mu + 2 * omega * categories * s)
  cell_list[[i]] <- data.frame(
    neighbor = s,
    category = categories,
    probability = weights / sum(weights)
  )
}
cells <- do.call(rbind, cell_list)

cat("distribution at neighbor = 0:",
    round(cells$probability[cells$neighbor == 0], 2), "\n")
cat("distribution at neighbor = 4:",
    round(cells$probability[cells$neighbor == 4], 2), "\n")

panel_labels <- paste("Neighbor at", neighbor_values)
names(panel_labels) <- neighbor_values

# Left panel: the two-node network these distributions belong to.
network_nodes <- data.frame(
  x = c(0, 0),
  y = c(1, 0),
  node = c("k", "j"),
  role = c("neighbor", "focal item")
)
net <- ggplot2::ggplot(network_nodes, ggplot2::aes(x = x, y = y)) +
  ggplot2::annotate(
    "segment",
    x = 0, xend = 0, y = 0.12, yend = 0.88,
    colour = "#12655F",
    linewidth = 1.5
  ) +
  ggplot2::geom_point(
    shape = 21,
    size = 10,
    fill = "#CDE5E2",
    colour = "#12655F",
    stroke = 0.8
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = node),
    colour = "#12655F",
    fontface = "italic",
    size = 3.6
  ) +
  ggplot2::geom_text(
    ggplot2::aes(y = y + c(0.52, -0.52), label = role),
    colour = "grey30",
    size = 3.1
  ) +
  ggplot2::annotate(
    "text",
    x = 0.3, y = 0.5,
    label = paste0("omega[jk] == ", omega),
    parse = TRUE,
    hjust = 0,
    colour = "#12655F",
    size = 3.4
  ) +
  ggplot2::coord_cartesian(xlim = c(-1.1, 2.0), ylim = c(-0.72, 1.72)) +
  ggplot2::theme_void()

fig <- ggplot2::ggplot(
  cells,
  ggplot2::aes(x = factor(category), y = probability)
) +
  ggplot2::geom_col(
    fill = "#CDE5E2",
    colour = "#12655F",
    linewidth = 0.4,
    width = 0.75
  ) +
  ggplot2::facet_wrap(
    ~neighbor,
    nrow = 1,
    labeller = ggplot2::labeller(neighbor = panel_labels)
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 0.62),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(
    x = "Answer on the focal item",
    y = "Probability"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 9.5),
    axis.title = ggplot2::element_text(size = 10.5)
  )

grDevices::pdf("figures/fig-omrf-intuition.pdf", width = 6.3, height = 1.48)
gridExtra::grid.arrange(net, fig, nrow = 1, widths = c(1.1, 5.2))
grDevices::dev.off()
