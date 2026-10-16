# 35-fig-toy-structures.R
# Concept figure for the three-symptom example (Section 2): prior-to-posterior
# updating. Two rows of four candidate
# structures with space between them. Top row: the probability bars sit ABOVE
# each network; bottom row: below. Each structure carries a pair of bars, its
# prior probability (light, 1/8 everywhere; the value .125 is printed once,
# above the empty network) and its posterior probability (dark). The four
# shaded cells hold the structures that contain the edge between intrusive
# thoughts (I) and concentration problems (C), drawn in the focal teal;
# summing their posterior bars gives the edge's posterior inclusion
# probability (the manuscript's Equation 1). No summary panel: the sum is
# stated in the caption.
#
# Style matches R/17 and R/18: focal teal #12655F, recessive greys, no axes.
# Reads tables/toy-structures.csv (from R/32-toy-structures.R).
# Output: figures/fig-toy-structures.pdf

library(ggplot2)

structures <- utils::read.csv("tables/toy-structures.csv",
                              colClasses = c(structure = "character"))

# Reading order: by number of edges, then enumeration.
order_wanted <- c("000", "100", "010", "001", "110", "101", "011", "111")
structures <- structures[match(order_wanted, structures$structure), ]

teal <- "#12655F"
grey_edge <- "grey40"
col_prior <- "grey80"
col_post <- "grey35"

pitch <- 1.5
bar_h <- 0.7

# Vertical anatomy. Bottom row: bars grow up from y = 0, network above them.
# Top row: network first, bars grow up from y = 3.35, so the bars sit on top.
bot_glyph <- data.frame(node = c("I", "S", "C"),
                        x = c(0.50, 0.13, 0.87), y = c(1.83, 1.35, 1.35))
top_glyph <- data.frame(node = c("I", "S", "C"),
                        x = c(0.50, 0.13, 0.87), y = c(3.03, 2.55, 2.55))
bot_base <- 0
top_base <- 3.35

edge_def <- data.frame(
  key = c("intrusion.sleep", "intrusion.concen", "sleep.concen"),
  from = c("I", "I", "S"),
  to = c("S", "C", "C"),
  focal = c(FALSE, TRUE, FALSE)
)

edges <- NULL
nodes <- NULL
bars <- NULL
labels <- NULL
shades <- NULL
baselines <- NULL

for (i in seq_len(nrow(structures))) {
  s <- structures[i, ]
  top <- i <= 4
  x0 <- ((i - 1) %% 4) * pitch
  glyph <- if (top) top_glyph else bot_glyph
  base <- if (top) top_base else bot_base
  n <- glyph
  n$x <- n$x + x0
  nodes <- rbind(nodes, n)
  for (j in seq_len(nrow(edge_def))) {
    if (s[[edge_def$key[j]]] == 1) {
      p1 <- n[n$node == edge_def$from[j], ]
      p2 <- n[n$node == edge_def$to[j], ]
      edges <- rbind(edges, data.frame(
        x = p1$x, y = p1$y, xend = p2$x, yend = p2$y,
        focal = edge_def$focal[j]
      ))
    }
  }
  if (s$intrusion.concen == 1) {
    shades <- rbind(shades, data.frame(
      xmin = x0 - 0.06, xmax = x0 + 1.06,
      ymin = if (top) 2.28 else -0.14,
      ymax = if (top) 3.62 else 2.08
    ))
  }
  bars <- rbind(bars,
    data.frame(xmin = x0 + 0.16, xmax = x0 + 0.44,
               ymin = base, ymax = base + bar_h * s$prior_prob),
    data.frame(xmin = x0 + 0.56, xmax = x0 + 0.84,
               ymin = base, ymax = base + bar_h * s$posterior_prob)
  )
  baselines <- rbind(baselines,
                     data.frame(x = x0 + 0.06, xend = x0 + 0.94, y = base))
  if (i == 1) {
    # The prior is uniform over the 2^3 = 8 structures, so every prior bar is
    # 1/8 = .125; the value is printed once, above the empty network.
    labels <- rbind(labels, data.frame(
      x = x0 + 0.30, y = base + bar_h * s$prior_prob + 0.12, text = ".125"
    ))
  }
  if (s$posterior_prob >= 0.0005) {
    labels <- rbind(labels, data.frame(
      x = x0 + 0.70,
      y = base + bar_h * s$posterior_prob + 0.12,
      text = sub("^0", "", sprintf("%.3f", s$posterior_prob))
    ))
  }
}

legend_df <- data.frame(
  xmin = c(0.02, 2.30), ymin = c(3.82, 3.82),
  fill = c(col_prior, col_post),
  text = c("prior probability", "posterior probability")
)

p <- ggplot() +
  geom_rect(data = shades,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = teal, alpha = 0.08) +
  geom_segment(data = baselines,
               aes(x = x, xend = xend, y = y, yend = y),
               linewidth = 0.3, color = "grey60") +
  geom_rect(data = bars,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = rep(c(col_prior, col_post), nrow(structures)),
            color = "grey45", linewidth = 0.25) +
  geom_segment(data = edges,
               aes(x = x, y = y, xend = xend, yend = yend),
               color = ifelse(edges$focal, teal, grey_edge),
               linewidth = ifelse(edges$focal, 1.0, 0.55)) +
  geom_point(data = nodes, aes(x = x, y = y),
             shape = 21, size = 7.0, fill = "white", color = "grey30",
             stroke = 0.45) +
  geom_text(data = nodes, aes(x = x, y = y, label = node),
            size = 2.3, color = "black") +
  geom_text(data = labels, aes(x = x, y = y, label = text),
            size = 2.8, color = "black") +
  geom_rect(data = legend_df,
            aes(xmin = xmin, xmax = xmin + 0.22, ymin = ymin,
                ymax = ymin + 0.16),
            fill = legend_df$fill, color = "grey45", linewidth = 0.25) +
  geom_text(data = legend_df,
            aes(x = xmin + 0.30, y = ymin + 0.08, label = text),
            size = 2.7, hjust = 0, color = "black") +
  coord_fixed(ratio = 1, clip = "off") +
  theme_void() +
  theme(plot.margin = margin(4, 4, 4, 4))

ggsave("figures/fig-toy-structures.pdf", p, width = 6.3, height = 4.6)

incl <- structures$intrusion.concen == 1
cat(sprintf(
  "written figures/fig-toy-structures.pdf (shaded prior sum %.3f -> posterior sum %.3f)\n",
  sum(structures$prior_prob[incl]), sum(structures$posterior_prob[incl])
))
