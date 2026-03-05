# =============================================================================
# 01_holo_som_analysis.R
# NS3 Protease HOLO — SOM Transition Network Analysis
#
# PRE-REQUISITE (run ONCE in terminal before this script):
#
#   gmx make_ndx -f data/md_100.gro -o data/index_holo.ndx
#   # type: a CA  then: q
#
#   gmx trjconv -f data/md_100_center.xtc \
#               -s data/md_100.gro \
#               -n data/index_holo.ndx \
#               -o data/holo_ca.pdb \
#               -skip 10
#   # Select CA group when prompted
#
# Then run this script from the NS3-SOM-MD/ root directory:
#   source("scripts/01_holo_som_analysis.R")
# =============================================================================

library(SOMMD)
library(bio3d)
library(kohonen)
library(igraph)

# =============================================================================
# PARAMETERS — adjust here if needed
# =============================================================================

HOLO_CA_PDB <- "data/holo_ca.pdb"
NC_CUTOFF   <- 8.0    # native contact distance cutoff (Angstrom)
SOM_X       <- 12     # SOM grid width  (12x12 = 144 neurons)
SOM_Y       <- 12     # SOM grid height
SOM_ITER    <- 500    # training iterations (increase to 1000 for better convergence)
RANDOM_SEED <- 42

# Output directories
dir.create("results/figures",    showWarnings = FALSE, recursive = TRUE)
dir.create("results/structures", showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# STEP 1 — Read Cα trajectory
# =============================================================================

cat("\n[1/6] Reading HOLO Cα trajectory...\n")

pdbs        <- read.pdb(HOLO_CA_PDB, multi = TRUE)
n_frames    <- dim(pdbs$xyz)[1]
n_ca        <- dim(pdbs$xyz)[2] / 3
xyz_all     <- pdbs$xyz

cat("  Frames loaded:", n_frames, "\n")
cat("  Cα atoms     :", n_ca, "\n")

# =============================================================================
# STEP 2 — Fit/align all frames to frame 1
# =============================================================================

cat("\n[2/6] Fitting trajectory to frame 1...\n")

all_cols <- seq_len(ncol(xyz_all))
xyz_fit  <- fit.xyz(fixed       = xyz_all[1, ],
                    mobile      = xyz_all,
                    fixed.inds  = all_cols,
                    mobile.inds = all_cols)
rm(xyz_all); gc()
cat("  Alignment done.\n")

# =============================================================================
# STEP 3 — Native contact distance matrix
# =============================================================================

cat("\n[3/6] Computing native contact distance matrix...\n")

ref_mat  <- matrix(xyz_fit[1, ], nrow = n_ca, ncol = 3, byrow = TRUE)
dist_ref <- as.matrix(dist(ref_mat))

idx    <- which(dist_ref < NC_CUTOFF & dist_ref > 0, arr.ind = TRUE)
idx    <- idx[idx[,1] < idx[,2], ]
idx    <- idx[abs(idx[,1] - idx[,2]) > 3, ]
n_cont <- nrow(idx)
cat("  Native contact pairs:", n_cont, "\n")

p1 <- (idx[,1] - 1) * 3 + 1
p2 <- (idx[,2] - 1) * 3 + 1

cat("  Computing distances across", n_frames, "frames...\n")
dist_mat <- matrix(0.0, nrow = n_frames, ncol = n_cont)
for (i in seq_len(n_frames)) {
  fr           <- xyz_fit[i, ]
  dx           <- fr[p1]   - fr[p2]
  dy           <- fr[p1+1] - fr[p2+1]
  dz           <- fr[p1+2] - fr[p2+2]
  dist_mat[i,] <- sqrt(dx*dx + dy*dy + dz*dz)
  if (i %% 200 == 0) cat("    frame", i, "/", n_frames, "\n")
}
rm(xyz_fit); gc()
cat("  Distance matrix:", nrow(dist_mat), "x", ncol(dist_mat), "\n")

# =============================================================================
# STEP 4 — Train SOM
# =============================================================================

cat("\n[4/6] Training HOLO SOM (", SOM_X, "x", SOM_Y, ",", SOM_ITER, "iters)...\n")
set.seed(RANDOM_SEED)

som_grid  <- somgrid(xdim = SOM_X, ydim = SOM_Y, topo = "hexagonal")
som_model <- som(dist_mat,
                 grid      = som_grid,
                 rlen      = SOM_ITER,
                 alpha     = c(0.05, 0.01),
                 keep.data = TRUE)
cat("  Training complete.\n")

# --- Plot 1: Convergence ---
png("results/figures/holo_plot1_convergence.png", width = 800, height = 500)
par(cex.main = 1.4, font.main = 2, cex.lab = 1.3, font.lab = 2, cex.axis = 1.2)
plot(som_model, type = "changes",
     main = "HOLO NS3 Protease — SOM Convergence")
dev.off(); cat("  Saved: results/figures/holo_plot1_convergence.png\n")

# --- Plot 2: U-matrix ---
png("results/figures/holo_plot2_umatrix.png", width = 800, height = 700)
par(cex.main = 1.4, font.main = 2)
plot(som_model, type = "dist.neighbours",
     main = "HOLO NS3 Protease — U-Matrix\n(dark = state boundaries)")
dev.off(); cat("  Saved: results/figures/holo_plot2_umatrix.png\n")

# --- Plot 3: Population ---
png("results/figures/holo_plot3_population.png", width = 800, height = 700)
par(cex.main = 1.4, font.main = 2)
plot(som_model, type = "counts",
     main = "HOLO NS3 Protease — Neuron Population")
dev.off(); cat("  Saved: results/figures/holo_plot3_population.png\n")

# =============================================================================
# STEP 5 — Transition network
# =============================================================================

cat("\n[5/6] Building HOLO transition network...\n")

assignments <- som_model$unit.classif
n_neurons   <- SOM_X * SOM_Y
pop         <- tabulate(assignments, nbins = n_neurons)

trans <- matrix(0L, nrow = n_neurons, ncol = n_neurons)
for (t in seq_len(n_frames - 1)) {
  f  <- assignments[t]; to <- assignments[t + 1]
  if (f != to) trans[f, to] <- trans[f, to] + 1L
}
rs         <- rowSums(trans)
trans_prob <- trans / ifelse(rs == 0, 1, rs)

net <- graph_from_adjacency_matrix(trans_prob, mode = "directed", weighted = TRUE)
net <- delete.edges(net, which(E(net)$weight < 1e-6))
net <- delete.vertices(net, which(pop == 0))

active_pop <- pop[pop > 0]
cat("  Active neurons:", length(active_pop), "/", n_neurons, "\n")
cat("  Network edges :", ecount(net), "\n")

pop_norm   <- (active_pop - min(active_pop)) / (max(active_pop) - min(active_pop) + 1)
pal        <- colorRampPalette(c("#2166ac", "#fee090", "#d73027"))(100)
node_color <- pal[pmax(1, ceiling(pop_norm * 100))]
node_size  <- 2 + sqrt(active_pop) * 0.9

png("results/figures/holo_plot4_transition_network.png", width = 1000, height = 950)
set.seed(RANDOM_SEED)
plot(net,
     layout          = layout_with_fr(net),
     vertex.size     = node_size,
     vertex.color    = node_color,
     vertex.label    = NA,
     edge.arrow.size = 0.15,
     edge.width      = E(net)$weight * 6,
     main = "HOLO NS3 Protease — Transition Network")
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
plot.new()
legend(x="bottom", legend=c("Rare states","Intermediate states","Frequent states"),
       col=c("#2166ac","#fee090","#d73027"), pch=19, pt.cex=3.0,
       horiz=TRUE, bty="o", box.col="#d73027", box.lwd=2.5,
       bg="white", cex=1.6, text.font=2, x.intersp=1.5, xpd=TRUE, inset=c(0,0.02))
dev.off(); cat("  Saved: results/figures/holo_plot4_transition_network.png\n")

# =============================================================================
# STEP 6 — Extract top 3 states per population tier
# =============================================================================

cat("\n[6/6] Extracting top 9 conformational states (3 per tier)...\n")

visited_idx <- which(pop > 0)
visited_pop <- pop[visited_idx]
q33 <- quantile(visited_pop, 0.33)
q66 <- quantile(visited_pop, 0.66)

blue_neurons   <- visited_idx[visited_pop <= q33]
yellow_neurons <- visited_idx[visited_pop > q33 & visited_pop <= q66]
red_neurons    <- visited_idx[visited_pop > q66]

get_top3 <- function(ids, p) head(ids[order(p[ids], decreasing = TRUE)], 3)

red_top3    <- get_top3(red_neurons,    pop)
yellow_top3 <- get_top3(yellow_neurons, pop)
blue_top3   <- get_top3(blue_neurons,   pop)
all_top9    <- c(red_top3, yellow_top3, blue_top3)

codebook <- som_model$codes[[1]]
get_rep  <- function(nid) {
  diffs <- sweep(dist_mat, 2, codebook[nid, ])
  which.min(rowSums(diffs^2))
}
rep_frames  <- sapply(all_top9, get_rep)
tier_names  <- c("Red_1","Red_2","Red_3","Yellow_1","Yellow_2","Yellow_3","Blue_1","Blue_2","Blue_3")
tier_labels <- c("R1","R2","R3","Y1","Y2","Y3","B1","B2","B3")

for (i in seq_along(all_top9)) {
  outfile <- paste0("results/structures/holo_state_", tier_names[i], ".pdb")
  write.pdb(pdb=pdbs, xyz=pdbs$xyz[rep_frames[i],], file=outfile)
  cat("  Saved:", outfile, "\n")
}

# 9-state annotated network
active_neurons <- which(pop > 0)
nc9  <- node_color; ns9 <- node_size
bc9  <- rep("white", length(active_pop)); bw9 <- rep(0, length(active_pop))
lv9  <- rep(NA, length(active_pop))
rc   <- c("#FF0000","#CC0000","#990000")
yc   <- c("#FFD700","#FFA500","#FF8C00")
bc_  <- c("#00BFFF","#1E90FF","#0000CD")
for (r in 1:3) {
  for (tier_info in list(list(red_top3,rc,"R"), list(yellow_top3,yc,"Y"), list(blue_top3,bc_,"B"))) {
    ni <- match(tier_info[[1]][r], active_neurons)
    if (!is.na(ni)) {
      nc9[ni] <- tier_info[[2]][r]; ns9[ni] <- ns9[ni]*1.6
      bc9[ni] <- "black"; bw9[ni] <- 2
      lv9[ni] <- paste0(tier_info[[3]], r)
    }
  }
}

png("results/figures/holo_plot5_9states_network.png", width = 1200, height = 1100)
set.seed(RANDOM_SEED)
plot(net, layout=layout_with_fr(net), vertex.size=ns9, vertex.color=nc9,
     vertex.frame.color=bc9, vertex.frame.width=bw9, vertex.label=lv9,
     vertex.label.color="black", vertex.label.font=2, vertex.label.cex=1.1,
     edge.arrow.size=0.12, edge.width=E(net)$weight*6,
     edge.color=adjustcolor("grey40",0.5),
     main="HOLO NS3 Protease — Top 9 Conformational States")
legend("bottomleft", legend=c("R1-R3: Frequent","Y1-Y3: Intermediate","B1-B3: Rare"),
       col=c("#FF0000","#FFD700","#00BFFF"), pch=19, pt.cex=1.8, bty="n", cex=1.1)
dev.off(); cat("  Saved: results/figures/holo_plot5_9states_network.png\n")

# Save population CSV
write.csv(data.frame(neuron=seq_len(n_neurons), population=pop),
          "results/structures/holo_representative_frames.csv", row.names=FALSE)

cat("\n=== HOLO ANALYSIS COMPLETE ===\n")
cat("  Frames      :", n_frames, "\n")
cat("  Contacts    :", n_cont, "\n")
cat("  Active nodes:", length(active_pop), "/", n_neurons, "\n")
cat("  Edges       :", ecount(net), "\n")
