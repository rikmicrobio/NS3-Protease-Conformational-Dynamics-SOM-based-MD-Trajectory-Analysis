# =============================================================================
# NS3 Protease APO — Full SOM Transition Network Analysis
# Matches exactly the holo analysis for direct comparison
#
# PRE-REQUISITE: Run these commands in your terminal first:
#
#   Step 1 — Create C-alpha index group:
#     gmx make_ndx -f md_100.gro -o index_apo.ndx
#     → type: a CA
#     → type: q
#
#   Step 2 — Extract Cα-only PDB from apo trajectory (every 10th frame):
#     gmx trjconv -f traj_0_80ns.xtc -s md_100.tpr \
#                 -n index_apo.ndx -o apo_ca.pdb -skip 10
#     → select the "CA" group when prompted
#
#   NOTE: .tpr is used as topology here (-s md_100.tpr) because it contains
#         the full system topology for the apo simulation.
#         If trjconv complains, try: -s md_100.gro instead
# =============================================================================

library(bio3d)
library(kohonen)
library(igraph)

# =============================================================================
# PARAMETERS — kept identical to holo for fair comparison
# =============================================================================

APO_CA_PDB  <- "data/apo_ca.pdb"     # output from gmx trjconv above
NC_CUTOFF   <- 8.0
SOM_X       <- 12
SOM_Y       <- 12
SOM_ITER    <- 500
RANDOM_SEED <- 42

# =============================================================================
# STEP 1 — Read apo Cα trajectory
# =============================================================================

cat("\n[1/6] Reading apo Cα trajectory...\n")

apo_pdbs    <- read.pdb(APO_CA_PDB, multi = TRUE)
n_frames_apo <- dim(apo_pdbs$xyz)[1]
n_ca_apo     <- dim(apo_pdbs$xyz)[2] / 3

cat("  Frames loaded:", n_frames_apo, "\n")
cat("  Cα atoms     :", n_ca_apo, "\n")

xyz_apo <- apo_pdbs$xyz

# =============================================================================
# STEP 2 — Fit/align to frame 1
# =============================================================================

cat("\n[2/6] Fitting apo trajectory to frame 1...\n")

all_cols    <- seq_len(ncol(xyz_apo))
xyz_apo_fit <- fit.xyz(fixed       = xyz_apo[1, ],
                       mobile      = xyz_apo,
                       fixed.inds  = all_cols,
                       mobile.inds = all_cols)
rm(xyz_apo); gc()
cat("  Alignment done.\n")

# =============================================================================
# STEP 3 — Native contact distance matrix
# =============================================================================

cat("\n[3/6] Computing native contact distance matrix...\n")

ref_mat  <- matrix(xyz_apo_fit[1, ], nrow = n_ca_apo, ncol = 3, byrow = TRUE)
dist_ref <- as.matrix(dist(ref_mat))

idx    <- which(dist_ref < NC_CUTOFF & dist_ref > 0, arr.ind = TRUE)
idx    <- idx[idx[,1] < idx[,2], ]
idx    <- idx[abs(idx[,1] - idx[,2]) > 3, ]
n_cont_apo <- nrow(idx)
cat("  Native contact pairs:", n_cont_apo, "\n")

p1 <- (idx[,1] - 1) * 3 + 1
p2 <- (idx[,2] - 1) * 3 + 1

cat("  Computing distances across", n_frames_apo, "frames...\n")
dist_mat_apo <- matrix(0.0, nrow = n_frames_apo, ncol = n_cont_apo)

for (i in seq_len(n_frames_apo)) {
  fr               <- xyz_apo_fit[i, ]
  dx               <- fr[p1]   - fr[p2]
  dy               <- fr[p1+1] - fr[p2+1]
  dz               <- fr[p1+2] - fr[p2+2]
  dist_mat_apo[i,] <- sqrt(dx*dx + dy*dy + dz*dz)
  if (i %% 200 == 0) cat("    frame", i, "/", n_frames_apo, "\n")
}
rm(xyz_apo_fit); gc()
cat("  Distance matrix:", nrow(dist_mat_apo), "x", ncol(dist_mat_apo), "\n")

# =============================================================================
# STEP 4 — Train SOM
# =============================================================================

cat("\n[4/6] Training SOM (", SOM_X, "x", SOM_Y, ",", SOM_ITER, "iters)...\n")
set.seed(RANDOM_SEED)

som_grid_apo  <- somgrid(xdim = SOM_X, ydim = SOM_Y, topo = "hexagonal")
som_model_apo <- som(dist_mat_apo,
                     grid      = som_grid_apo,
                     rlen      = SOM_ITER,
                     alpha     = c(0.05, 0.01),
                     keep.data = TRUE)
cat("  Training complete.\n")

png("results/figures/apo_plot1_convergence.png", width = 800, height = 500)
plot(som_model_apo, type = "changes",
     main = "APO NS3 Protease — SOM Convergence")
dev.off(); cat("  Saved: apo_plot1_convergence.png\n")

png("results/figures/apo_plot2_umatrix.png", width = 800, height = 700)
plot(som_model_apo, type = "dist.neighbours",
     main = "APO NS3 Protease — U-Matrix\n(dark = state boundaries)")
dev.off(); cat("  Saved: apo_plot2_umatrix.png\n")

png("results/figures/apo_plot3_population.png", width = 800, height = 700)
plot(som_model_apo, type = "counts",
     main = "APO NS3 Protease — Neuron Population")
dev.off(); cat("  Saved: apo_plot3_population.png\n")

# =============================================================================
# STEP 5 — Transition network
# =============================================================================

cat("\n[5/6] Building apo transition network...\n")

assignments_apo <- som_model_apo$unit.classif
n_neurons       <- SOM_X * SOM_Y
pop_apo         <- tabulate(assignments_apo, nbins = n_neurons)

trans_apo <- matrix(0L, nrow = n_neurons, ncol = n_neurons)
for (t in seq_len(n_frames_apo - 1)) {
  f  <- assignments_apo[t]
  to <- assignments_apo[t + 1]
  if (f != to) trans_apo[f, to] <- trans_apo[f, to] + 1L
}
rs_apo        <- rowSums(trans_apo)
trans_prob_apo <- trans_apo / ifelse(rs_apo == 0, 1, rs_apo)

net_apo <- graph_from_adjacency_matrix(trans_prob_apo, mode = "directed", weighted = TRUE)
net_apo <- delete.edges(net_apo, which(E(net_apo)$weight < 1e-6))
net_apo <- delete.vertices(net_apo, which(pop_apo == 0))

active_pop_apo <- pop_apo[pop_apo > 0]
cat("  Active neurons:", length(active_pop_apo), "/", n_neurons, "\n")
cat("  Network edges :", ecount(net_apo), "\n")

pop_norm_apo   <- (active_pop_apo - min(active_pop_apo)) /
                  (max(active_pop_apo) - min(active_pop_apo) + 1)
pal            <- colorRampPalette(c("#2166ac", "#fee090", "#d73027"))(100)
node_color_apo <- pal[pmax(1, ceiling(pop_norm_apo * 100))]
node_size_apo  <- 2 + sqrt(active_pop_apo) * 0.9

png("results/figures/apo_plot4_transition_network.png", width = 1000, height = 950)
set.seed(RANDOM_SEED)
plot(net_apo,
     layout          = layout_with_fr(net_apo),
     vertex.size     = node_size_apo,
     vertex.color    = node_color_apo,
     vertex.label    = NA,
     edge.arrow.size = 0.15,
     edge.width      = E(net_apo)$weight * 6,
     main = "APO NS3 Protease — Transition Network\nBlue=rare  |  Red=frequent")
legend("bottomleft",
       legend = c("Rare state", "Intermediate", "Frequent state"),
       col    = c("#2166ac", "#fee090", "#d73027"),
       pch = 19, pt.cex = 1.8, bty = "n", cex = 1.1)
dev.off(); cat("  Saved: apo_plot4_transition_network.png\n")

# =============================================================================
# STEP 6 — Extract top 3 states per tier (Blue / Yellow / Red)
# =============================================================================

cat("\n[6/6] Extracting top 3 states per population tier...\n")

visited_idx_apo <- which(pop_apo > 0)
visited_pop_apo <- pop_apo[visited_idx_apo]

q33_apo <- quantile(visited_pop_apo, 0.33)
q66_apo <- quantile(visited_pop_apo, 0.66)

blue_neurons_apo   <- visited_idx_apo[visited_pop_apo <= q33_apo]
yellow_neurons_apo <- visited_idx_apo[visited_pop_apo > q33_apo & visited_pop_apo <= q66_apo]
red_neurons_apo    <- visited_idx_apo[visited_pop_apo > q66_apo]

get_top3_neurons <- function(neuron_ids, pop_vec) {
  head(neuron_ids[order(pop_vec[neuron_ids], decreasing = TRUE)], 3)
}

red_top3_apo    <- get_top3_neurons(red_neurons_apo,    pop_apo)
yellow_top3_apo <- get_top3_neurons(yellow_neurons_apo, pop_apo)
blue_top3_apo   <- get_top3_neurons(blue_neurons_apo,   pop_apo)
all_top9_apo    <- c(red_top3_apo, yellow_top3_apo, blue_top3_apo)

codebook_apo <- som_model_apo$codes[[1]]
get_rep_frame_apo <- function(neuron_id) {
  diffs <- sweep(dist_mat_apo, 2, codebook_apo[neuron_id, ])
  which.min(rowSums(diffs^2))
}
rep_frames_apo <- sapply(all_top9_apo, get_rep_frame_apo)

tier_names_apo <- c("Red_1","Red_2","Red_3",
                    "Yellow_1","Yellow_2","Yellow_3",
                    "Blue_1","Blue_2","Blue_3")
tier_labels_apo <- c("R1","R2","R3","Y1","Y2","Y3","B1","B2","B3")

# Save PDB structures
cat("  Saving PDB structures...\n")
for (i in seq_along(all_top9_apo)) {
  outfile <- paste0("results/structures/apo_state_", tier_names_apo[i], ".pdb")
  write.pdb(pdb  = apo_pdbs,
            xyz  = apo_pdbs$xyz[rep_frames_apo[i], ],
            file = outfile)
  cat(sprintf("    Saved: %s (neuron %d, frame %d, pop=%d)\n",
              outfile, all_top9_apo[i],
              rep_frames_apo[i], pop_apo[all_top9_apo[i]]))
}

# 9-state network plot
active_neurons_apo <- which(pop_apo > 0)
node_color_9  <- node_color_apo
node_size_9   <- node_size_apo
border_col_9  <- rep("white", length(active_pop_apo))
border_wid_9  <- rep(0,       length(active_pop_apo))
label_vec_9   <- rep(NA,      length(active_pop_apo))

red_cols    <- c("#FF0000","#CC0000","#990000")
yellow_cols <- c("#FFD700","#FFA500","#FF8C00")
blue_cols   <- c("#00BFFF","#1E90FF","#0000CD")

assign_hl <- function(neuron_ids, colors, prefix) {
  for (r in seq_along(neuron_ids)) {
    ni <- match(neuron_ids[r], active_neurons_apo)
    if (!is.na(ni)) {
      node_color_9[ni]  <<- colors[r]
      node_size_9[ni]   <<- node_size_9[ni] * 1.6
      border_col_9[ni]  <<- "black"
      border_wid_9[ni]  <<- 2
      label_vec_9[ni]   <<- paste0(prefix, r)
    }
  }
}
assign_hl(red_top3_apo,    red_cols,    "R")
assign_hl(yellow_top3_apo, yellow_cols, "Y")
assign_hl(blue_top3_apo,   blue_cols,   "B")

png("results/figures/apo_plot5_9states_network.png", width = 1200, height = 1100)
set.seed(RANDOM_SEED)
plot(net_apo,
     layout             = layout_with_fr(net_apo),
     vertex.size        = node_size_9,
     vertex.color       = node_color_9,
     vertex.frame.color = border_col_9,
     vertex.frame.width = border_wid_9,
     vertex.label       = label_vec_9,
     vertex.label.color = "black",
     vertex.label.font  = 2,
     vertex.label.cex   = 1.1,
     edge.arrow.size    = 0.12,
     edge.width         = E(net_apo)$weight * 6,
     edge.color         = adjustcolor("grey40", alpha.f = 0.5),
     main = "APO NS3 Protease — Conformational Transition Network\nTop 3 States per Population Tier")
legend("bottomleft",
       legend = c("R1–R3: Frequent", "Y1–Y3: Intermediate", "B1–B3: Rare"),
       col    = c("#FF0000","#FFD700","#00BFFF"),
       pch = 19, pt.cex = 1.8, bty = "n", cex = 1.05)
dev.off(); cat("  Saved: apo_plot5_9states_network.png\n")

# Population bar chart
state_pops_apo <- pop_apo[all_top9_apo]
state_pcts_apo <- round(state_pops_apo / sum(pop_apo) * 100, 1)
bar_cols_9     <- c(red_cols, yellow_cols, blue_cols)

png("results/figures/apo_plot6_9states_barchart.png", width = 900, height = 600)
par(mar = c(6, 5, 4, 2))
bp <- barplot(state_pcts_apo,
              names.arg = tier_labels_apo,
              col       = bar_cols_9,
              border    = "black",
              ylim      = c(0, max(state_pcts_apo) * 1.3),
              ylab      = "% of Trajectory Frames",
              main      = "APO NS3 Protease — Population of Top 9 Conformational States",
              cex.names = 1.2, cex.axis = 1.1, las = 1)
text(bp, state_pcts_apo + max(state_pcts_apo) * 0.03,
     paste0(state_pcts_apo, "%"), cex = 1.0, font = 2)
abline(v = c(mean(bp[3:4])), lty = 2, col = "grey50")
abline(v = c(mean(bp[6:7])), lty = 2, col = "grey50")
legend("topright",
       legend = c("Frequent (Red)","Intermediate (Yellow)","Rare (Blue)"),
       fill   = c("#FF0000","#FFD700","#1E90FF"),
       bty = "n", cex = 1.1)
dev.off(); cat("  Saved: apo_plot6_9states_barchart.png\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n=======================================================\n")
cat("  APO NS3 Protease — 9 Key Conformational States\n")
cat("=======================================================\n")
cat(sprintf("  %-6s %-10s %-10s %-8s %-10s %-12s\n",
            "Label","Tier","Neuron","Frames","% Traj","Rep.Frame"))
cat(sprintf("  %s\n", paste(rep("-",62), collapse="")))
tier_type_apo <- c(rep("Frequent",3), rep("Interm.",3), rep("Rare",3))
for (i in 1:9) {
  cat(sprintf("  %-6s %-10s %-10d %-8d %-10.1f %-12d\n",
              tier_labels_apo[i], tier_type_apo[i],
              all_top9_apo[i], pop_apo[all_top9_apo[i]],
              round(pop_apo[all_top9_apo[i]] / sum(pop_apo) * 100, 1),
              rep_frames_apo[i]))
}
cat("-------------------------------------------------------\n")
cat("  PDB files : apo_state_Red/Yellow/Blue_1/2/3.pdb\n")
cat("  Plots     : apo_plot1-6_*.png\n")
cat("=======================================================\n")
cat("\n  NEXT STEP: Run apo_vs_holo_comparison.R to compare\n")
cat("  the two conformational landscapes side by side.\n")
