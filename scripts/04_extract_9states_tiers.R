# =============================================================================
# NS3 Protease — Extract Top 3 States from Each Population Tier
# Blue  (rare)         → top 3 least populated neurons
# Yellow (intermediate) → top 3 mid-range populated neurons
# Red   (frequent)     → top 3 most populated neurons
#
# Run in the SAME R session as sommd_analysis_v6.R
# Requires: som_model, dist_mat, pdbs, SOM_X, SOM_Y, RANDOM_SEED
# =============================================================================

library(bio3d)
library(igraph)

cat("\n=======================================================\n")
cat("  NS3 Protease — 9-State Conformational Tier Analysis\n")
cat("=======================================================\n")

# =============================================================================
# STEP 1 — Classify all neurons into tiers by population
# =============================================================================

assignments <- som_model$unit.classif
n_neurons   <- SOM_X * SOM_Y
pop         <- tabulate(assignments, nbins = n_neurons)

# Only consider visited neurons
visited_idx <- which(pop > 0)
visited_pop <- pop[visited_idx]

# Divide into terciles: bottom 33% = blue/rare, middle = yellow/intermediate, top 33% = red/frequent
q33 <- quantile(visited_pop, 0.33)
q66 <- quantile(visited_pop, 0.66)

blue_neurons   <- visited_idx[visited_pop <= q33]
yellow_neurons <- visited_idx[visited_pop > q33 & visited_pop <= q66]
red_neurons    <- visited_idx[visited_pop > q66]

cat(sprintf("\n  Population thresholds:\n"))
cat(sprintf("  Blue  (rare,         pop ≤ %d): %d neurons\n", round(q33), length(blue_neurons)))
cat(sprintf("  Yellow (intermediate, %d < pop ≤ %d): %d neurons\n", round(q33), round(q66), length(yellow_neurons)))
cat(sprintf("  Red   (frequent,     pop > %d): %d neurons\n", round(q66), length(red_neurons)))

# Top 3 within each tier (by population)
get_top3 <- function(neuron_ids, pop_vec, tier_label) {
  tier_pop <- pop_vec[neuron_ids]
  ranked   <- neuron_ids[order(tier_pop, decreasing = TRUE)]
  top3_ids <- head(ranked, 3)
  pcts     <- round(pop_vec[top3_ids] / sum(pop_vec) * 100, 2)
  cat(sprintf("\n  %s tier — Top 3 neurons:\n", tier_label))
  for (r in 1:3) {
    cat(sprintf("    %s%d: Neuron %3d | %d frames | %.1f%% of trajectory\n",
                substr(tier_label, 1, 1), r,
                top3_ids[r], pop_vec[top3_ids[r]], pcts[r]))
  }
  return(top3_ids)
}

red_top3    <- get_top3(red_neurons,    pop, "Red   (Frequent)")
yellow_top3 <- get_top3(yellow_neurons, pop, "Yellow (Intermediate)")
blue_top3   <- get_top3(blue_neurons,   pop, "Blue  (Rare)")

all_top9 <- c(red_top3, yellow_top3, blue_top3)

# =============================================================================
# STEP 2 — Representative frames for all 9 states
# =============================================================================

cat("\n--- Computing representative frames ---\n")

codebook <- som_model$codes[[1]]

get_rep_frame <- function(neuron_id) {
  diffs <- sweep(dist_mat, 2, codebook[neuron_id, ])
  which.min(rowSums(diffs^2))
}

rep_frames <- sapply(all_top9, get_rep_frame)
tier_labels <- c("R1","R2","R3","Y1","Y2","Y3","B1","B2","B3")
tier_names  <- c("Red_1","Red_2","Red_3",
                 "Yellow_1","Yellow_2","Yellow_3",
                 "Blue_1","Blue_2","Blue_3")

# =============================================================================
# STEP 3 — Save representative PDB structures
# =============================================================================

cat("\n--- Saving PDB structures for all 9 states ---\n")

for (i in seq_along(all_top9)) {
  outfile <- paste0("results/structures/holo_state_", tier_names[i], ".pdb")
  write.pdb(pdb  = pdbs,
            xyz  = pdbs$xyz[rep_frames[i], ],
            file = outfile)
  cat(sprintf("  Saved: %s  (neuron %d, frame %d, pop=%d)\n",
              outfile, all_top9[i], rep_frames[i], pop[all_top9[i]]))
}

# =============================================================================
# STEP 4 — Transition network with all 9 states highlighted
# =============================================================================

cat("\n--- Building enhanced transition network ---\n")

active_neurons <- which(pop > 0)
active_pop_vec <- pop[active_neurons]
n_active       <- length(active_neurons)

# Build transition matrix
trans <- matrix(0L, nrow = n_neurons, ncol = n_neurons)
n_frames <- nrow(dist_mat)
for (t in seq_len(n_frames - 1)) {
  f  <- assignments[t]
  to <- assignments[t + 1]
  if (f != to) trans[f, to] <- trans[f, to] + 1L
}
rs         <- rowSums(trans)
trans_prob <- trans / ifelse(rs == 0, 1, rs)

net <- graph_from_adjacency_matrix(trans_prob, mode = "directed", weighted = TRUE)
net <- delete.edges(net, which(E(net)$weight < 1e-6))
net <- delete.vertices(net, which(pop == 0))

# Base node colours from population gradient
pop_norm   <- (active_pop_vec - min(active_pop_vec)) /
              (max(active_pop_vec) - min(active_pop_vec) + 1)
base_pal   <- colorRampPalette(c("#2166ac", "#fee090", "#d73027"))(100)
node_color <- base_pal[pmax(1, ceiling(pop_norm * 100))]
node_size  <- 2 + sqrt(active_pop_vec) * 0.9

# Highlight colours for each tier
red_cols    <- c("#FF0000", "#CC0000", "#990000")       # bright → dark red
yellow_cols <- c("#FFD700", "#FFA500", "#FF8C00")       # gold shades
blue_cols   <- c("#00BFFF", "#1E90FF", "#0000CD")       # bright → dark blue

# Border colours to make highlighted nodes pop
border_col  <- rep("white", n_active)
border_width <- rep(0, n_active)

label_vec   <- rep(NA, n_active)

assign_highlight <- function(neuron_ids, colors, prefix) {
  for (r in seq_along(neuron_ids)) {
    net_idx <- match(neuron_ids[r], active_neurons)
    if (!is.na(net_idx)) {
      node_color[net_idx]  <<- colors[r]
      node_size[net_idx]   <<- node_size[net_idx] * 1.6   # make highlighted bigger
      border_col[net_idx]  <<- "black"
      border_width[net_idx] <<- 2
      label_vec[net_idx]   <<- paste0(prefix, r)
    }
  }
}

assign_highlight(red_top3,    red_cols,    "R")
assign_highlight(yellow_top3, yellow_cols, "Y")
assign_highlight(blue_top3,   blue_cols,   "B")

png("results/figures/holo_9states_network.png", width = 1200, height = 1100)
set.seed(RANDOM_SEED)
layout_net <- layout_with_fr(net)
plot(net,
     layout             = layout_net,
     vertex.size        = node_size,
     vertex.color       = node_color,
     vertex.frame.color = border_col,
     vertex.frame.width = border_width,
     vertex.label       = label_vec,
     vertex.label.color = "black",
     vertex.label.font  = 2,
     vertex.label.cex   = 1.1,
     edge.arrow.size    = 0.12,
     edge.width         = E(net)$weight * 6,
     edge.color         = adjustcolor("grey40", alpha.f = 0.5),
     main = "NS3 Protease Conformational Transition Network\nTop 3 States per Population Tier")

legend("bottomleft",
       legend = c("R1–R3: Frequent (dominant) states",
                  "Y1–Y3: Intermediate states",
                  "B1–B3: Rare (transient) states",
                  "── background nodes"),
       col    = c("#FF0000", "#FFD700", "#00BFFF", "#aaaaaa"),
       pch    = 19, pt.cex = 1.8, bty = "n", cex = 1.05,
       title  = "State Tiers")
dev.off()
cat("  Saved: plot_9states_network.png\n")

# =============================================================================
# STEP 5 — Population bar chart for all 9 states
# =============================================================================

cat("--- Plotting population bar chart ---\n")

state_pops <- pop[all_top9]
state_pcts <- round(state_pops / sum(pop) * 100, 1)
bar_cols   <- c(red_cols, yellow_cols, blue_cols)

png("results/figures/holo_9states_barchart.png", width = 900, height = 600)
par(mar = c(6, 5, 4, 2))
bp <- barplot(state_pcts,
              names.arg = tier_labels,
              col       = bar_cols,
              border    = "black",
              ylim      = c(0, max(state_pcts) * 1.25),
              ylab      = "% of Trajectory Frames",
              xlab      = "",
              main      = "NS3 Protease — Population of Top 9 Conformational States",
              cex.names = 1.2,
              cex.axis  = 1.1,
              las       = 1)

# Value labels on bars
text(bp, state_pcts + max(state_pcts) * 0.03,
     paste0(state_pcts, "%"), cex = 1.0, font = 2)

# X-axis tier labels
mtext("State Label", side = 1, line = 4, cex = 1.1)

# Tier brackets
abline(v = c(mean(bp[3:4])), lty = 2, col = "grey50")
abline(v = c(mean(bp[6:7])), lty = 2, col = "grey50")

legend("topright",
       legend = c("Frequent (Red)", "Intermediate (Yellow)", "Rare (Blue)"),
       fill   = c("#FF0000", "#FFD700", "#1E90FF"),
       bty    = "n", cex = 1.1)
dev.off()
cat("  Saved: plot_9states_barchart.png\n")

# =============================================================================
# FINAL SUMMARY TABLE
# =============================================================================

cat("\n=======================================================\n")
cat("  NS3 Protease — 9 Key Conformational States\n")
cat("=======================================================\n")
cat(sprintf("  %-6s %-8s %-10s %-8s %-10s %-14s\n",
            "Label", "Tier", "Neuron", "Frames", "% Traj", "Rep. Frame"))
cat(sprintf("  %s\n", paste(rep("-", 60), collapse="")))

tier_type <- c(rep("Frequent", 3), rep("Interm.", 3), rep("Rare", 3))
for (i in 1:9) {
  cat(sprintf("  %-6s %-8s %-10d %-8d %-10.1f %-14d\n",
              tier_labels[i],
              tier_type[i],
              all_top9[i],
              pop[all_top9[i]],
              round(pop[all_top9[i]] / sum(pop) * 100, 1),
              rep_frames[i]))
}
cat("-------------------------------------------------------\n")
cat("  PDB files saved: state_Red_1/2/3.pdb\n")
cat("                   state_Yellow_1/2/3.pdb\n")
cat("                   state_Blue_1/2/3.pdb\n")
cat("  Plots saved:     plot_9states_network.png\n")
cat("                   plot_9states_barchart.png\n")
cat("=======================================================\n")
