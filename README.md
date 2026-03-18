# Conformational Dynamics — SOM-based MD Trajectory Analysis

[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue?logo=r)](https://www.r-project.org/)
[![SOMMD](https://img.shields.io/badge/SOMMD-CRAN-brightgreen)](https://cran.r-project.org/web/packages/SOMMD/index.html)
[![GROMACS](https://img.shields.io/badge/GROMACS-2020%2B-orange)](https://www.gromacs.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> Self-Organizing Map (SOM) analysis of NS3 protease molecular dynamics trajectories — APO vs HOLO (ligand-bound) conformational landscape comparison.

---

![Image](https://github.com/user-attachments/assets/a42bdd21-0a75-498a-8005-60e00db94a7b)

## Table of Contents

- [Overview](#overview)
- [Scientific Background](#scientific-background)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Input Files](#input-files)
- [Workflow](#workflow)
  - [Step 0 — GROMACS Preprocessing](#step-0--gromacs-preprocessing)
  - [Step 1 — HOLO Analysis](#step-1--holo-analysis)
  - [Step 2 — APO Analysis](#step-2--apo-analysis)
  - [Step 3 — APO vs HOLO Comparison](#step-3--apo-vs-holo-comparison)
  - [Step 4 — Dominant State Extraction](#step-4--dominant-state-extraction)
- [Output Files](#output-files)
- [Interpreting the Results](#interpreting-the-results)
- [Citation](#citation)

---

## Overview

This repository contains a complete, reproducible pipeline for analysing the conformational dynamics of **NS3 protease** from molecular dynamics (MD) simulations using Self-Organizing Maps (SOM), as implemented in the [SOMMD](https://github.com/alepandini/SOMMD) R package.

The analysis covers:
- **Conformational clustering** of MD trajectory frames via SOM
- **Transition network construction** to map pathways between conformational states
- **Tier-based state extraction** — identifying rare, intermediate, and dominant conformations
- **APO vs HOLO comparison** — quantifying the effect of ligand binding on the conformational landscape

---

## Scientific Background

Self-Organizing Maps are unsupervised neural networks that project high-dimensional structural data (Cα pairwise distances) onto a 2D grid. Each node (neuron) on the grid represents a cluster of structurally similar frames. The resulting **transition network** captures the probability of moving between conformational states during the simulation, providing a graph-based view of the energy landscape.

**Key metrics produced:**
| Metric | Biological meaning |
|--------|-------------------|
| Neuron population | How much time the protein spends in each conformation |
| U-matrix | Boundaries between distinct conformational states |
| Transition edges | Kinetic connectivity between states |
| Representative frames | Actual MD snapshots for each state — loadable in PyMOL/VMD |

---

## Repository Structure

```
NS3-SOM-MD/
│
├── README.md                          # This file
│
├── scripts/
│   ├── 01_holo_som_analysis.R         # HOLO SOM + transition network
│   ├── 02_apo_som_analysis.R          # APO SOM + transition network
│   ├── 03_apo_vs_holo_comparison.R    # Side-by-side comparison plots
│   └── 04_extract_9states_tiers.R     # Top 3 states per population tier
│
├── data/
│   ├── md_100.gro                     # HOLO topology (GROMACS GRO)
│   ├── md_100_center.xtc              # HOLO trajectory (GROMACS XTC)
│   ├── md_100.tpr                     # APO topology (GROMACS TPR)
│   └── traj_0_80ns.xtc                # APO trajectory (GROMACS XTC)
│
├── results/
│   ├── figures/                       # All generated plots (PNG)
│   └── structures/                    # Representative PDB files per state
│
└── LICENSE
```

> **Note:** Trajectory files (`.xtc`, `.tpr`, `.gro`) are not tracked by git due to size.  
> Add `*.xtc`, `*.tpr`, `*.dcd` to your `.gitignore`.

---

## Requirements

### R packages

```r
install.packages(c("SOMMD", "bio3d", "kohonen", "igraph"), dependencies = TRUE)
```

| Package | Version tested | Purpose |
|---------|---------------|---------|
| `SOMMD` | ≥ 1.0 | SOM training, trajectory handling |
| `bio3d` | ≥ 2.4 | Structure/trajectory I/O, alignment |
| `kohonen` | ≥ 3.0 | SOM engine |
| `igraph` | ≥ 1.3 | Transition network construction |

### External tools

- **GROMACS** ≥ 2020 — required for trajectory preprocessing only
- **PyMOL** or **VMD** — recommended for visualising output PDB structures

---

## Input Files

| File | System | Description |
|------|--------|-------------|
| `md_100.gro` | HOLO | Equilibrated structure with ligand (LIG), solvent, ions |
| `md_100_center.xtc` | HOLO | 100 ns trajectory, centred and PBC-corrected |
| `md_100.tpr` | APO | GROMACS run input (topology + parameters) |
| `traj_0_80ns.xtc` | APO | 80 ns apo trajectory |

**System composition (HOLO):**
- 199 protein residues (3025 atoms)
- 1 ligand molecule (LIG)
- Solvent (SOL) + ions (NA, CL)
- Total: 32,988 atoms

---

## Workflow

### Step 0 — GROMACS Preprocessing

Convert trajectories to Cα-only multi-model PDB format. This step is done **once in the terminal** and drastically reduces memory usage in R.

```bash
# --- HOLO ---
# Create Cα index group
gmx make_ndx -f md_100.gro -o index_holo.ndx
# type: a CA  →  q

# Extract Cα frames (every 10th = ~880 frames)
gmx trjconv -f md_100_center.xtc \
            -s md_100.gro \
            -n index_holo.ndx \
            -o holo_ca.pdb \
            -skip 10
# Select CA group when prompted

# --- APO ---
gmx make_ndx -f md_100.gro -o index_apo.ndx
# type: a CA  →  q

gmx trjconv -f traj_0_80ns.xtc \
            -s md_100.tpr \
            -n index_apo.ndx \
            -o apo_ca.pdb \
            -skip 10
# Select CA group when prompted
```

> **Why PDB?** bio3d's `read.xtc` is not available in standard installations. The multi-model PDB approach is fully portable and requires no additional libraries.

---

### Step 1 — HOLO Analysis

```r
source("scripts/01_holo_som_analysis.R")
```

**What it does:**
1. Reads `holo_ca.pdb` (Cα frames)
2. Fits all frames to frame 1 (Kabsch alignment via `bio3d::fit.xyz`)
3. Computes native contact distance matrix (Cα pairs within 8 Å, non-local)
4. Trains a 12×12 hexagonal SOM (500 iterations)
5. Builds directed transition network
6. Extracts top 3 dominant states per population tier (Red/Yellow/Blue)

**Key parameters (adjustable at top of script):**

```r
NC_CUTOFF   <- 8.0    # native contact distance threshold (Angstrom)
SOM_X       <- 12     # SOM grid width  (increase for larger datasets)
SOM_Y       <- 12     # SOM grid height
SOM_ITER    <- 500    # training iterations (increase to 1000 for better convergence)
RANDOM_SEED <- 42     # reproducibility seed
```

---

### Step 2 — APO Analysis

```r
source("scripts/02_apo_som_analysis.R")
```

Identical pipeline to Step 1, applied to the apo trajectory. All parameters are kept identical for a fair comparison.

---

### Step 3 — APO vs HOLO Comparison

```r
# Run in the SAME R session as Steps 1 and 2
source("scripts/03_apo_vs_holo_comparison.R")
```

Generates side-by-side comparison figures and a summary statistics table.

---

### Step 4 — Dominant State Extraction

```r
# Run in the SAME R session as Steps 1 and 2
source("scripts/04_extract_9states_tiers.R")
```

Extracts the top 3 conformational states from each population tier (9 states total per system), saves representative PDB structures, and generates annotated network plots.

---

## Output Files

### Figures

| File | Description |
|------|-------------|
| `plot1_convergence.png` | SOM training convergence curve |
| `plot2_umatrix.png` | U-matrix — conformational state boundaries |
| `plot3_population.png` | Neuron population heatmap |
| `plot4_transition_network.png` | Full transition network |
| `plot_9states_network.png` | Network with top 9 states highlighted |
| `plot_9states_barchart.png` | Population bar chart for top 9 states |
| `comparison_convergence.png` | APO vs HOLO convergence overlaid |
| `comparison_population_maps.png` | Side-by-side population maps |
| `comparison_umatrix.png` | Side-by-side U-matrices |
| `comparison_networks.png` | Side-by-side transition networks |

### Structures

| File | Description |
|------|-------------|
| `state_Red_1/2/3.pdb` | HOLO dominant (frequent) conformations |
| `state_Yellow_1/2/3.pdb` | HOLO intermediate conformations |
| `state_Blue_1/2/3.pdb` | HOLO rare/transient conformations |
| `apo_state_Red_1/2/3.pdb` | APO dominant conformations |
| `apo_state_Yellow_1/2/3.pdb` | APO intermediate conformations |
| `apo_state_Blue_1/2/3.pdb` | APO rare conformations |
| `representative_frames.csv` | Frame indices and populations for all neurons |

### Visualising structures in PyMOL

```python
# Load and compare APO vs HOLO dominant states
load apo_state_Red_1.pdb
load state_Red_1.pdb
align apo_state_Red_1, state_Red_1
show cartoon
color blue, apo_state_Red_1
color red,  state_Red_1
```

---

## Interpreting the Results

### SOM Convergence (plot1)
The mean distance to the closest unit should **decrease and flatten** by the end of training. If the curve has not plateaued, increase `SOM_ITER` to 1000.

### U-Matrix (plot2)
- **Light/yellow regions** → neurons with similar codebook vectors = within a conformational state
- **Dark/red regions** → large distances between neighbouring neurons = state boundaries

### Population Map (plot3)
- **Red neurons** → frequently visited conformations (dominant states)
- **Grey/white neurons** → rarely or never visited

### Transition Network (plot4)
- **Node size** → population (time spent in state)
- **Edge width** → transition probability
- **Densely connected core** → dominant conformational ensemble
- **Long branches** → transition pathways to rare/excited states

### APO vs HOLO Comparison
| Observation | Interpretation |
|-------------|---------------|
| Fewer active neurons in HOLO | Ligand restricts conformational space |
| More transition edges in HOLO | Ligand promotes dynamic interconversion |
| Different dominant state positions | Ligand shifts the preferred conformation |
| Lower final SOM error in APO | APO ensemble is more homogeneous |

---

## Citation

If you use this pipeline, please cite:

**SOMMD:**
```
Motta, S., Callea, L., Bonati, L., Pandini, A. (2022).
PathDetect-SOM: A Neural Network Approach for the Identification of
Pathways in Ligand Binding Simulations.
Journal of Chemical Theory and Computation, 18(3), 1957–1968.
https://doi.org/10.1021/acs.jctc.1c01163
```

**bio3d:**
```
Grant, B.J., Rodrigues, A.P.C., ElSawy, K.M., McCammon, J.A., Caves, L.S.D. (2006).
Bio3d: An R package for the comparative analysis of protein structures.
Bioinformatics, 22(21), 2695–2696.
https://doi.org/10.1093/bioinformatics/btl461
```

---

## .gitignore

Add this to your `.gitignore` to avoid committing large binary files:

```gitignore
# Trajectory and topology files
*.xtc
*.trr
*.dcd
*.tpr
*.edr
*.cpt

# Large PDB trajectory files
*_ca.pdb

# R workspace
.RData
.Rhistory
*.Rproj.user/
```

---

*Pipeline developed for NS3 protease conformational dynamics analysis.*  
*Contributions and issues welcome.*
