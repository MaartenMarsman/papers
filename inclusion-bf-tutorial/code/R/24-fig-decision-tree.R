# Known-truth calibration study, step 5: the decision-tree figure.
# Design: the simulation appendix of the manuscript.
#
# One tree per parameter class, in the vocabulary of the manuscript's "Did the
# sampler do its job?" subsection: the ESS pair, flip counts, chain agreement,
# settled edges. The TOPOLOGY is declared here because the shape of the advice
# is a scientific choice; the NUMBERS on every leaf are read back from
# tables/threshold-study-rules-*.csv so the figure cannot drift from the
# calibration behind it.
#
# Colors: the manuscript's own accent pair (R/01-worked-example.R), checked
# for colorblind safety: worst-pair CVD delta-E 11.9 and normal-vision
# delta-E 20.3, both above their thresholds. The pair is flagged for chroma because #12655F is a deliberately
# muted teal used across every figure in this manuscript; a more saturated
# substitute scores WORSE on colorblind separation (delta-E 7.3), so the house
# color is also the more legible one. Every node carries a text label, so
# identity never rests on color alone.
#
# Output: figures/fig-diagnostic-decision-tree.pdf

# Overridable so the renderer can be exercised against smoke-test tables.
input_prefix <- Sys.getenv("THRESHOLD_STUDY_OUTPUT", "tables/")
figure_path <- Sys.getenv("THRESHOLD_STUDY_FIGURE",
                          "figures/fig-diagnostic-decision-tree.pdf")

certify_color <- "#12655F"
escalate_color <- "#A66A1E"
node_fill <- "#F2F0EC"
node_line <- "#8C877D"
ink <- "#2B2B2B"

rules_indicator <- utils::read.csv(paste0(input_prefix, "threshold-study-pair-rule.csv"),
                                   stringsAsFactors = FALSE)
rules_weight <- utils::read.csv(paste0(input_prefix, "threshold-study-rules-weight.csv"),
                                stringsAsFactors = FALSE)
rules_continuous <- utils::read.csv(
  paste0(input_prefix, "threshold-study-rules-continuous.csv"),
  stringsAsFactors = FALSE)

# Leaf annotation: the pass rate and conditional error rate the calibration
# measured for the rule this branch implements.
annotate_rule <- function(table, rule_name) {
  row <- table[table$rule == rule_name, ]
  if (!nrow(row)) {
    stop(sprintf("rule \"%s\" not found in the rules table (available: %s)",
                 rule_name, paste(table$rule, collapse = "; ")),
         call. = FALSE)
  }
  sprintf("passes %.1f%%; %s in %s\n(95%% upper bound %.2f%%)",
          100 * row$pass_rate[1],
          if (row$errors[1] == 0) "no errors" else
            sprintf("%d errors", row$errors[1]),
          format(row$n_pass[1], big.mark = ","),
          100 * row$error_upper95[1])
}

# --- Flowchart renderer --------------------------------------------------------

# The layout works in INCHES, with coord_fixed() pinning one data unit to one
# inch and the device sized to match. Sizing boxes in arbitrary data units
# does not survive: text is laid out in points, so a box tall enough at one
# figure height clips its own label at another. In inches the arithmetic is
# the same arithmetic the text renderer uses.
#
# All three trees share this one coordinate system rather than being stitched
# together as separate plots, which would give each its own scale and render
# identical boxes at different physical sizes.

figure_width <- 6.5
figure_height <- 8.7
text_size <- 2.5              # ggplot2 text size is in mm
line_inch <- 0.135            # line advance at that size, with padding
# Column spacing has to leave a real connector between boxes: at box_width
# 1.75 on a 1.9-inch pitch the arrows were 0.15 inches long and the branch
# labels were clipped by the boxes at either end.
box_width <- 1.55

line_count <- function(label) vapply(strsplit(label, "\n"), length, 0L)
box_height <- function(label) 0.16 + line_inch * line_count(label)

tree_layers <- function(nodes, edges, y_offset) {
  nodes$y <- nodes$y + y_offset
  nodes$height <- box_height(nodes$label)
  nodes$fill <- ifelse(nodes$type == "certify", certify_color,
                       ifelse(nodes$type == "escalate", escalate_color,
                              node_fill))
  nodes$text_color <- ifelse(nodes$type == "check", ink, "white")

  segments <- merge(edges, nodes[, c("id", "x", "y")], by.x = "from",
                    by.y = "id")
  names(segments)[names(segments) == "x"] <- "x_from"
  names(segments)[names(segments) == "y"] <- "y_from"
  segments <- merge(segments, nodes[, c("id", "x", "y")], by.x = "to",
                    by.y = "id")
  names(segments)[names(segments) == "x"] <- "x_to"
  names(segments)[names(segments) == "y"] <- "y_to"
  segments$x_start <- segments$x_from + box_width / 2
  segments$x_end <- segments$x_to - box_width / 2
  # Branch labels sit at the midpoint of the connector and just above it:
  # placed nearer either end they collide with the box they are leaving or
  # entering.
  segments$x_label <- (segments$x_start + segments$x_end) / 2
  segments$y_label <- (segments$y_from + segments$y_to) / 2 + 0.09

  # The title of a tree is placed from the tree, not guessed: hard-coded
  # heights land on top of whichever box happens to be tallest.
  list(nodes = nodes, segments = segments,
       top = max(nodes$y + nodes$height / 2))
}

draw_forest <- function(trees) {
  nodes <- do.call(rbind, lapply(trees, function(t) t$nodes))
  segments <- do.call(rbind, lapply(trees, function(t) t$segments))
  titles <- do.call(rbind, lapply(trees, function(t) t$title))

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = segments,
      ggplot2::aes(x = x_start, y = y_from, xend = x_end, yend = y_to),
      colour = node_line, linewidth = 0.4,
      arrow = ggplot2::arrow(length = ggplot2::unit(0.06, "inches"),
                             type = "closed")
    ) +
    ggplot2::geom_text(
      data = segments,
      ggplot2::aes(x = x_label, y = y_label, label = label),
      size = 2.3, colour = node_line, fontface = "italic"
    ) +
    ggplot2::geom_tile(
      data = nodes,
      ggplot2::aes(x = x, y = y, height = height, fill = I(fill)),
      width = box_width, colour = node_line, linewidth = 0.3
    ) +
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(x = x, y = y, label = label, colour = I(text_color)),
      size = text_size, lineheight = 1.0
    ) +
    ggplot2::geom_text(
      data = titles,
      ggplot2::aes(x = x, y = y, label = label), hjust = 0, size = 3,
      fontface = "bold", colour = ink
    ) +
    ggplot2::coord_fixed(ratio = 1, xlim = c(0, figure_width),
                         ylim = c(0, figure_height), expand = FALSE) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(plot.margin = ggplot2::margin(2, 2, 2, 2))
}

# --- The three trees -------------------------------------------------------------

# The calibrated rules. Each name must match the `rule` column of the table it
# is looked up in, and the drawn tree must implement the named rule so the
# leaf annotations describe what is drawn. The indicator rule is chosen on
# error catching rather than on a conditional error target, and it requires
# both channels: the two-state flag reads the INDICATOR chain's exploration,
# the RB rule reads the PRECISION of the Rao-Blackwellized draws, and the two
# have complementary blind spots.
indicator_rule_name <- "PAIR: both channels pass"
weight_rule_name <- "weight n_eff >= 1000"
continuous_rule_name <- "rhat < 1.010"

indicator_nodes <- data.frame(
  id = c("exploration", "precision", "certify", "escalate"),
  x = c(0.85, 3.05, 5.25, 5.25),
  y = c(1.50, 2.10, 2.75, 0.85),
  label = c(
    "Exploration channel:\nis every verdict boundary\nmore than 2 SE from the\nsmoothed logit inclusion\nprobability?",
    "Precision channel:\nis it also more than 2 RB\nstandard errors away?\n(settled edge: do all\nchains agree instead?)",
    "Read the verdict",
    "Run longer if the indicator\nflipped; if it is settled,\ndisperse the starts;\nlonger runs cannot help"
  ),
  type = c("check", "check", "certify", "escalate"),
  stringsAsFactors = FALSE
)
indicator_edges <- data.frame(
  from = c("exploration", "exploration", "precision", "precision"),
  to = c("precision", "escalate", "certify", "escalate"),
  label = c("yes", "no", "yes", "no"),
  stringsAsFactors = FALSE
)

weight_nodes <- data.frame(
  id = c("root", "certify", "escalate"),
  x = c(0.85, 3.05, 3.05), y = c(0.55, 1.05, 0.05),
  label = c(
    "Model-averaged weight:\nis the derived composite\nESS at least 1,000?",
    "Report the weight\nwith its interval", "Run longer"
  ),
  type = c("check", "certify", "escalate"), stringsAsFactors = FALSE
)
weight_edges <- data.frame(from = c("root", "root"),
                           to = c("certify", "escalate"),
                           label = c("yes", "no"), stringsAsFactors = FALSE)

continuous_nodes <- data.frame(
  id = c("root", "certify", "escalate"),
  x = c(0.85, 3.05, 3.05), y = c(0.55, 1.05, 0.05),
  label = c(
    "Category thresholds and\nprecisions: is split-R-hat\nbelow 1.01?",
    "Read the estimate", "Run longer"
  ),
  type = c("check", "certify", "escalate"), stringsAsFactors = FALSE
)
continuous_edges <- weight_edges

# Attach the measured pass and error rates to the certifying leaf of each tree.
indicator_nodes$label[indicator_nodes$id == "certify"] <- paste(
  indicator_nodes$label[indicator_nodes$id == "certify"],
  annotate_rule(rules_indicator, indicator_rule_name), sep = "\n")
# The indicator rule is read from the table pooled over all arms, written by
# R/25-threshold-study-transfer.R; the other two classes are read from the
# rule tables written by R/23-threshold-study-analysis.R.

weight_nodes$label[weight_nodes$id == "certify"] <- paste(
  weight_nodes$label[weight_nodes$id == "certify"],
  annotate_rule(rules_weight, weight_rule_name), sep = "\n")
continuous_nodes$label[continuous_nodes$id == "certify"] <- paste(
  continuous_nodes$label[continuous_nodes$id == "certify"],
  annotate_rule(rules_continuous, continuous_rule_name), sep = "\n")

# Vertical offsets stack the three trees in one column, tallest first.
with_title <- function(tree, label) {
  tree$title <- data.frame(x = 0.05, y = tree$top + 0.22, label = label,
                           stringsAsFactors = FALSE)
  tree
}

indicator_tree <- with_title(
  tree_layers(indicator_nodes, indicator_edges, y_offset = 4.95),
  "Edge indicators: the three-way verdict")
weight_tree <- with_title(
  tree_layers(weight_nodes, weight_edges, y_offset = 2.45),
  "Model-averaged weights (discrete arms; misses about half on GGMs)")
continuous_tree <- with_title(
  tree_layers(continuous_nodes, continuous_edges, y_offset = 0.35),
  "Category thresholds and precisions")

figure <- draw_forest(list(indicator_tree, weight_tree, continuous_tree))

ggplot2::ggsave(figure_path, plot = figure, width = figure_width,
                height = figure_height)
cat("wrote", figure_path, "\n")
