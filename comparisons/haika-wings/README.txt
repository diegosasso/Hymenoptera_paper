HAIKA and ontophylo wing rates on the 77-tip tree

Completed 20 November 2026. Read analysis-report.txt for the interpretation and
results/primary-summary.csv for the compact numerical comparison. The original
77-tip family-proxy analysis has weak primary branch concordance (rho=0.068);
temporal correspondence varies with ANN initialization, sampling and smoothing.

This analysis keeps the study's original dated 77-tip tree and ontophylo map
samples. HAIKA image labels can be matched to all 48 historical morphology
families represented by those tips. The resulting 77-tip HAIKA dataset is a
family-proxy alignment, not 77 measured species. The user authorized these
proxies and requested an alignment as close to the original tree as possible.

Source
Peng et al. (2026), Science Advances, doi:10.1126/sciadv.aef0605.
Data archive: https://zenodo.org/records/22715169
Ontophylo repository: https://github.com/diegosasso/Hymenoptera_paper
See audit-haika.txt, audit-ontophylo.txt, and provenance/ for checksums,
source discrepancies, anatomical definitions, and computational verification.

Data choices
The original image-labelled 64-dimensional matrices from five archived ANN seeds are
used separately. Seed 3407 is the prespecified primary example. The released
28-dimensional image matrix is a sensitivity. Family image means receive equal weight
in PCA; retain at least two PCs and at least 90% of family-level variance.
These are new PCA axes, not the publication's shape/venation axes.

The released Case 3 species matrix has unresolved cross-taxon profile matches
and does not reproduce the published PC scores. The S4 genus matrix repeats
family templates. Both are retained for audit, not used as verified species
measurements. Source-image family assignments are explicit in
data/haika-family-proxy-crosswalk.csv and haika-family-proxy-source-rows.csv.

Models and comparison
Continuous: exact phytools 2.5.2 multirateBM penalized likelihood, evaluated with
a verified cached objective and analytic gradient. Three bounded starts per
fit; lambda 0.1, 1, and 10 for the primary seed. Node rates are converted to branch
averages using the logarithmic mean. Each PC is fitted separately; summed
diffusivities describe total retained variance, not an independently fitted
joint model. All fits require convergence and independent likelihood checks.

Discrete: existing 1,000 stochastic-map draws, with branch rate equal to summed
character transition counts divided by branch duration. Fore-wing (9 characters)
is primary, fore-wing sclerites (8 characters) an anatomical sensitivity; all 19 wing characters
and hind-wing (10 characters) are secondary. The eight external fore-wing
characters describe seven vein/stigma features and one nygmata feature; the
"sclerite" label is the repository's broad character-type classification.
They overlap substantially with HAIKA's wing-outline and venation information.
The ninth fore-wing character is a muscle. The representations still differ:
discrete anatomical states versus continuous image geometry and latent features.

Normalize each rate realization by its tree-length-weighted mean. Use identical
branches and branch-duration exposure to calculate 10-Myr marginal averages.
The narrow residual bin above 280 Ma is exported but excluded from correlations.
Map intervals are conditional on fitted discrete models; network-seed and
penalty variation are sensitivities, not HAIKA confidence or credible bands.

Controls
Two 48-family representative subsets assess repeated family-profile assignments
and representative choice. Full 77-tip interfamily branches are reported separately.
The 199-replicate null simulates homogeneous multivariate BM on the 48-family
tree, retaining fitted cross-PC covariance, then copies the simulated family
profiles to the same 77 positions. It calibrates correlations against this
specific tree/proxy null, conditional on the fixed PCA and fitted covariance.
It does not test general independence of two heterogeneous evolutionary histories.

Execution from this directory
python3 scripts/01-extract-haika.py
python3 scripts/01-extract-haika-family-crosswalk.py
Rscript --vanilla scripts/02-extract-ontophylo.R
Rscript --vanilla scripts/03-prepare-comparison.R
Rscript --vanilla scripts/verify-bm-objective.R
Run sbatch slurm/batch.slurm on Roihu: scripts/run-batch.py dispatches all
51 planned fits and 199 null replicates through scripts/04-fit-bm.R and
scripts/06-null-comparison.R. Run sbatch slurm/bounds.slurm for 18 additional
numerical-bound and strong-smoothing fits (scripts/09-bound-sensitivity.R).
Rscript --vanilla scripts/05-compare-rates.R --final
Rscript --vanilla scripts/07-collect-null.R --final
Rscript --vanilla scripts/08-plot-comparison.R
python3 scripts/10-write-report.py

Roihu
Root: /scratch/project_2004996/haika-wings
Module: r-env/461; project: project_2004996; partition: interactive.
Completed job 3285454: 51 fits and 199 null replicates, eight CPUs, 62 seconds,
approximately 2.5 GB maximum resident memory, exit code 0.
Completed job 3285525: 18 diagnostic fits, one CPU, four seconds,
approximately 0.4 GB maximum resident memory, exit code 0.
All 69 fits and 199 replicates passed their required checks. CSC permits batch
jobs in this partition: https://docs.csc.fi/supercomputing/running/interactive-usage/.
Earlier task-specific jobs 3285349, 3285385, 3285386 and 3285424 were cancelled
while pending and superseded by the completed batch allocation. Existing
unrelated jobs were not changed. See logs/, provenance/roihu-accounting.txt,
results/fit-status.csv, results/null-status.csv and results/bound-axis-diagnostics.csv.

The additional smoothing diagnostic uses lambda equal to the original root age
(about 280 Myr), with lower log-rate offsets 5 and 15 checked separately at
lambda 1. The small rates and detailed branch ranks are sensitive to these
choices. Their results are separate from the original fit manifest and the
lambda-1 null calibration. All family assignments remain explicit proxies.
