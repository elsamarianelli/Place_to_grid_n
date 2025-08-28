This repository turns place-cell population activity into grid-like spatial maps using dimensionality-reduction. There are two main analysis routes:

1) Simulate exploration, build a neuron×neuron covariance matrix from place-cell firing, run PCA, and project components back to 2D space to obtain grid-like rate maps.
  
3) Simulate a foraging trajectory, learn a Successor Representation (SR) matrix from place-cell activity with temporal-difference updates, take leading eigenvectors, and project them back into space to yield grid-like maps.
Whilst the current write up includes the primary results using the Covariance PCA path without non-negativity constraints, this repo also includes several experimental options (below) that were explored but not used in the paper’s main analyses.

_**The 2 main scripts to run in this repo are:**_ 

**1) Generate data (covariance PCA or SR)- 'Algorithmic_PCA/SR_and_Covar_Param_sweeps_master.m'**

Inside the script you can toggle:
- Environment (square / trapezoid)
- Place-cell population type (uniform, boundary_modulated, size_only, density_only)
- Dimensionality reduction mode (covariance vs SR)
- PCA version (standard PCA or NNPCA/FISTA variants; see below)
- Trajectory type (for SR you can use uniform or the foraging model)

(pipeline visualisation...)
<img width="500" alt="image" src="https://github.com/user-attachments/assets/ce5eb587-69f1-4f87-8e5c-ed1303355b99" />

_**Pipelines**_

**1) Covariance pipeline**
- Simulate place cells in a 2D environment (square or trapezoid). Boundary-modulated variants control field size and/or centre distribution as a function of distance to walls; uniform variants keep these constant.
- Uniform trajectory: record each neuron’s time series, zero-mean, and build the covariance matrix (neurons × neurons).
- PCA with MATLAB pca(...,'Algorithm','eig','Centered',false). Select top components (e.g., 250).
- Project back to space using the PC weights over the original 2D place fields to obtain grid-like maps.
- Compute spatial metrics per bin (3×3 in squares; halves in trapezoid) and aggregate by region; save CSVs and figures. 

**2) SR pipeline**

- Foraging model (velocity-based): at 50 Hz for 360,000 steps (~2 h), heading noise ~ 𝒩(0, 0.2 rad), speed ~ Rayleigh(b = 16), capped at 50 units/step; turn away and slow within 7 units of walls to ensure realistic sampling.
- Record place-cell firing over the trajectory.
- TD learning of SR matrix M using place-cell firing as successor features (α=1, discount γ):
- Eigen-decompose 𝑀. Because 𝑀 is near-symmetric but slightly asymmetric from TD updates, retain real parts of eigenvectors.
- Project eigenvectors to space (reweighted input maps) to obtain grid-like outputs. (No non-negativity post-constraint is applied in this SR path, so “square” gridness can exceed hex gridness.)

**NNPCA / FISTA options** - was playing around with this and did manage to generate some grid like representations with the merged version, requires many place cells however (see below)
Set the PCA_type toggle inside the master script to pick one of the non-negative versions:
- sharp_assymptotics — Implements the “sharp asymptotics” approach (slower; benefits from more place cells and the “space invader” boundaries).
- FISTA_DOG — Follows the Dordek et al. use-case: FISTA with a Difference-of-Gaussians penalty; tweak Gaussian sigmas to shift resulting grid scales. (see examples generated below)
- 
  <img width="500" alt="image" src="https://github.com/user-attachments/assets/385f2cf2-2de4-4a03-bd99-3d841702ee74" />

- FISTA_merged — Merges Wang & Wang’s cost function with the DOG penalty; includes a faster FISTA implementation adapted from a public GitHub. Consider increasing the number of place cells for best results. (see examples generated below with the following parameters...)
_use_SR  =  false;            use_traj = 'uniform' or 'generate';                                    trap_add = 0;                In.shape = 'trapezoid';      NonNegative = true;       PCA_type ='FISTA_merged';                pc_density = false  ;        pc_size    =  false;        mean_firing_match = true;  In.pf_width_cntrl = 2;      n_iterations = 1;In.n_cells = 500.*5;       In.n_steps = 360000;      In.dim_x = 351;             In.dim_y = 252;In.n_polys = 1;In.NumberOfPC = 250; In.bound_ctrl = 2_

<img width="500"  alt="image" src="https://github.com/user-attachments/assets/10d814e4-af80-4f7c-93c7-ad46cd021b05" />

_Extra notes on these variations:_
- Covariance-FISTA (Wang & Wang-style)
Works on the neuron × neuron covariance built from place-cell activity.
Finds a non-negative unit vector v that maximises variance (then deflates the rank-1 contribution and repeats.)
- Spatial DoG-FISTA (Dordek-style)
Optimises directly in image space (a 2D component map)
Minimises Difference-of-Gaussians (DoG) energy, which acts as a band-pass prior that encourages periodic, grid-like patterns at controllable scales. Also uses non-negativity + unit-norm projection, plus deflation to diversify components. Intuition: explicitly favor the spatial frequency content that looks grid-like.


**2) Generate Plots - 'final_figures_master.m'**

- Bar plots (for square arena: Corners / Edges / Centre; for trapezoid: Left / Right). Bars are means of iteration means; error bars are SD across iterations.
- Heatmaps (square): 3×3 bin maps averaged across PCs and iterations for each metric.
 
Stats
- Square arena: paired Student’s t-tests on iteration means for each pair: Corners–Edges, Corners–Centre, Edges–Centre. Two-tailed α = .05. Reports t, df, two-sided p, 95% CI, and paired Cohen’s dₓ (dz).
- Trapezoid: paired t-test on iteration means Left vs Right, with the same outputs.

(example figures for square case)

<img width="500" alt="image" src="https://github.com/user-attachments/assets/b3d1d04b-4a8c-4be3-b0f7-8cf677e45f72" />
<img width="500" alt="image" src="https://github.com/user-attachments/assets/32fb2c47-4ee7-47f1-8334-125221125dde" />

_**Extras**_

**Real data analysis of recorded grids**
real_data_analysis/real_data_analysis.m
For running the same grid-metric pipeline on recorded datasets 
- Load & preprocess real sessions (positions + spikes → smoothed rate maps per neuron).
- Compute spatial autocorrelograms (SACs) and metrics
<img width="500"  alt="image" src="https://github.com/user-attachments/assets/bbcdd573-54e8-4591-aff3-4615e2db7bff" />

**Neural-network PCA (Hebbian / Oja’s rule)**
An online, Hebbian network that learns principal-component–like filters from place-cell inputs using a variant of Oja’s rule with an optional non-negativity constraint on the feedforward weights. 
The model can also include lateral connections in the grid-cell layer to encourage competition/structure. 
Learned weights are projected back onto the place-cell maps to visualise emergent grid-like patterns.
(Network.m contains code to emualte nn PCA using a hebbian network with non-negative weights. Simplified from Dordek version, set up to use Tanni warped PCs and Hasselmo trajectory instead
_Note_ I was playing around with this but never managed to get it to generate good grids, as noted, dordek uses very specific parameters for place cell input which may have caused problems here)
