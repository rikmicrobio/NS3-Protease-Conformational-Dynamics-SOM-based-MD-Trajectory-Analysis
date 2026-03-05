# =============================================================================
# NS3 Protease — APO vs HOLO Conformational Comparison
# Run AFTER both sommd_analysis_v6.R AND apo_analysis.R
# Requires all objects from both sessions to be in memory
# =============================================================================

library(igraph)

cat("\n=======================================================\n")
cat("  NS3 Protease — APO vs HOLO Comparison\n")
cat("=======================================================\n")

# =============================================================================
# PLOT 1 — Side-by-side population maps
# =============================================================================

png("results/figures/comparison_population_maps.png", width = 1400, height = 700)
par(mfrow = c(1, 2), cex.main = 1.6, font.main = 2, cex.lab = 1.4, font.lab = 2, cex.axis = 1.2)
plot(som_model_apo, type = "counts",  main = "APO NS3 Protease\nNeuron Population")
plot(som_model,     type = "counts",  main = "HOLO NS3 Protease (+ Ligand)\nNeuron Population")
dev.off()
cat("  Saved: comparison_population_maps.png\n")

# =============================================================================
# PLOT 2 — Side-by-side U-matrices
# =============================================================================

png("results/figures/comparison_umatrix.png", width = 1400, height = 700)
par(mfrow = c(1, 2), cex.main = 1.6, font.main = 2, cex.lab = 1.4, font.lab = 2, cex.axis = 1.2)
plot(som_model_apo, type = "dist.neighbours", main = "APO NS3 Protease\nU-Matrix (State Boundaries)")
plot(som_model,     type = "dist.neighbours", main = "HOLO NS3 Protease\nU-Matrix (State Boundaries)")
dev.off()
cat("  Saved: comparison_umatrix.png\n")

# =============================================================================
# PLOT 3 — Side-by-side transition networks
# =============================================================================

build_net <- function(assignments, pop, n_neurons) {
  trans <- matrix(0L, n_neurons, n_neurons)
  for (t in seq_len(length(assignments) - 1)) {
    f <- assignments[t]; to <- assignments[t+1]
    if (f != to) trans[f, to] <- trans[f, to] + 1L
  }
  rs  <- rowSums(trans)
  tp  <- trans / ifelse(rs == 0, 1, rs)
  net <- graph_from_adjacency_matrix(tp, mode = "directed", weighted = TRUE)
  net <- delete.edges(net, which(E(net)$weight < 1e-6))
  net <- delete.vertices(net, which(pop == 0))
  net
}

n_neurons    <- SOM_X * SOM_Y
pop_holo     <- tabulate(som_model$unit.classif,     nbins = n_neurons)
pop_apo_full <- tabulate(som_model_apo$unit.classif, nbins = n_neurons)

net_holo <- build_net(som_model$unit.classif,     pop_holo,     n_neurons)
net_apo2 <- build_net(som_model_apo$unit.classif, pop_apo_full, n_neurons)

style_net <- function(net, pop_full) {
  ap  <- pop_full[pop_full > 0]
  pn  <- (ap - min(ap)) / (max(ap) - min(ap) + 1)
  pal <- colorRampPalette(c("#2166ac","#fee090","#d73027"))(100)
  list(color = pal[pmax(1, ceiling(pn*100))],
       size  = 2 + sqrt(ap) * 0.9)
}

s_apo  <- style_net(net_apo2, pop_apo_full)
s_holo <- style_net(net_holo, pop_holo)

png("results/figures/comparison_networks.png", width = 1600, height = 900)
par(mfrow = c(1, 2), mar = c(2, 2, 5, 2),
    cex.main = 1.8, font.main = 2)
set.seed(RANDOM_SEED)

plot(net_apo2,
     layout          = layout_with_fr(net_apo2),
     vertex.size     = s_apo$size,
     vertex.color    = s_apo$color,
     vertex.label    = NA,
     edge.arrow.size = 0.12,
     edge.width      = E(net_apo2)$weight * 6,
     edge.color      = adjustcolor("grey40", 0.5),
     main = "APO NS3 Protease\nTransition Network")

plot(net_holo,
     layout          = layout_with_fr(net_holo),
     vertex.size     = s_holo$size,
     vertex.color    = s_holo$color,
     vertex.label    = NA,
     edge.arrow.size = 0.12,
     edge.width      = E(net_holo)$weight * 6,
     edge.color      = adjustcolor("grey40", 0.5),
     main = "HOLO NS3 Protease (+ Ligand)\nTransition Network")

# Centered bottom legend spanning both panels
par(fig = c(0, 1, 0, 1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
plot.new()
legend(x         = "bottom",
       legend    = c("Rare states", "Intermediate states", "Frequent states"),
       col       = c("#2166ac", "#fee090", "#d73027"),
       pch       = 19,
       pt.cex    = 3.0,
       horiz     = TRUE,
       bty       = "o",
       box.col   = "#d73027",
       box.lwd   = 2.5,
       bg        = "white",
       cex       = 1.6,
       text.font = 2,
       x.intersp = 1.5,
       y.intersp = 1.2,
       xpd       = TRUE,
       inset     = c(0, 0.02))

dev.off()
cat("  Saved: comparison_networks.png\n")

# =============================================================================
# PLOT 4 — Convergence curves overlaid
# =============================================================================

png("results/figures/comparison_convergence.png", width = 900, height = 550)
par(mar = c(6, 6, 4, 2))

plot(som_model_apo$changes,
     type = "l", col = "#2166ac", lwd = 2,
     xlab = "",
     ylab = "",
     main = "SOM Convergence — APO vs HOLO NS3 Protease",
     cex.main = 1.4,
     cex.axis = 1.2,
     ylim = range(c(som_model_apo$changes, som_model$changes)))

# Bold axis labels drawn separately for font control
title(xlab = "Iteration",
      ylab = "Mean Distance to Closest Unit",
      cex.lab = 1.5,
      font.lab = 2)

lines(som_model$changes, col = "#d73027", lwd = 2)

legend("topright",
       legend = c("APO","HOLO (+ Ligand)"),
       col    = c("#2166ac","#d73027"),
       lwd = 2, bty = "n", cex = 1.3)

dev.off()
cat("  Saved: comparison_convergence.png\n")

# =============================================================================
# SUMMARY STATISTICS COMPARISON TABLE
# =============================================================================

apo_active  <- sum(pop_apo_full > 0)
holo_active <- sum(pop_holo > 0)

apo_edges   <- ecount(net_apo2)
holo_edges  <- ecount(net_holo)

apo_frames  <- nrow(dist_mat_apo)
holo_frames <- nrow(dist_mat)

cat("\n=======================================================\n")
cat("  APO vs HOLO — Conformational Landscape Comparison\n")
cat("=======================================================\n")
cat(sprintf("  %-30s %-12s %-12s\n", "Metric", "APO", "HOLO"))
cat(sprintf("  %s\n", paste(rep("-", 56), collapse="")))
cat(sprintf("  %-30s %-12d %-12d\n", "Frames analysed",       apo_frames,   holo_frames))
cat(sprintf("  %-30s %-12d %-12d\n", "Active SOM neurons",    apo_active,   holo_active))
cat(sprintf("  %-30s %-12d %-12d\n", "Transition edges",      apo_edges,    holo_edges))
cat(sprintf("  %-30s %-12.4f %-12.4f\n", "Final SOM error",
            tail(som_model_apo$changes, 1), tail(som_model$changes, 1)))
cat("=======================================================\n")
cat("\n  Interpretation guide:\n")
cat("  - More active neurons in APO → broader conformational sampling\n")
cat("  - More edges in HOLO → more interconnected / dynamic landscape\n")
cat("  - Higher APO error → more heterogeneous structural ensemble\n")
cat("  - Compare dominant states (R1-R3) PDB files in PyMOL:\n")
cat("    load apo_state_Red_1.pdb; load state_Red_1.pdb\n")
cat("    align apo_state_Red_1, state_Red_1\n")
cat("=======================================================\n")
