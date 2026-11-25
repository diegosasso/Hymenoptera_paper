# Phylogenetic sampling simulations

This analysis tests how reducing the 765-tip Hymenoptera tree of Blaimer et al.
(2023) to the 77 study taxa affects reconstructed phenomic rates. It compares
the full and study trees using the same simulated character origins.

## Design

Each selected internal node defines one binary character: its descendants have
state 1 and all other taxa have state 0. Origins are restricted to crown
Apocrita, excluding its incoming stem. Nodes are sampled uniformly without
replacement within each epoch. There are 100 trials per condition, giving
1,500 paired trials.

| Epoch (Ma) | Eligible nodes | Characters at 5%, 10%, 30%, 50%, 80% |
| --- | ---: | --- |
| 0–80 | 540 | 27, 54, 162, 270, 432 |
| 80–160 | 177 | 9, 18, 53, 89, 142 |
| 160–280 | 17 | 1, 2, 5, 9, 14 |

Intervals include their younger bound and exclude their older bound. Counts
are rounded to the nearest integer, with a minimum of one character.

Both trees use an equal-rates Mk model, the `maddfitz` root treatment, five
optimization starts, and 100 stochastic maps per variable character. Each
unique node pattern is fitted once per tree; trials receive fresh mapping
draws with saved seeds. Characters made constant by pruning are omitted.

Mapped changes are summed across characters and divided by branch duration.
Marginal rates average these combined branch rates over the lineages present
at each age. Rates are not divided by the number of characters.

`RESULTS.txt` gives the numerical findings and limitations. The main figures
compare all simulated characters on the full tree with the characters that
remain variable in the study sample. The matched figures use the same
surviving characters on both trees.

## Contents

- `inputs/`: reference tree, study tree and taxon correspondence.
- `design/`: eligible nodes, character patterns, fitted models and trial manifest.
- `scripts/`: simulation, inference, plotting and validation code.
- `slurm/`: job scripts for Roihu.
- `summary/`: numerical summaries and analysis checks.
- `figures/`: final PDF figures and PNG previews.
- `manuscript/`: LaTeX methods, captions and references.
- `provenance/`: frozen ontophylo source, package versions and checksums.

The development version of ontophylo used for the analysis is included in
`provenance/ontophylo`. Its source revision is recorded in
`provenance/ontophylo-git-head.txt`. Package author and licence information
is retained with the source.

Per-trial results, simulated selections, fit checkpoints and the large
`summary/marginal-trials.rds` object are local outputs excluded from Git.
Existing local copies are retained. A fresh checkout can regenerate them
with the scripts below; no download of archived raw results is configured.

## Inspect the results

Open `phylo-sampling.Rproj` in RStudio, or run commands from this directory.
The saved summaries and design objects are sufficient to regenerate the rate
figures:

```sh
Rscript scripts/06-plot-results.R
```

Plotting requires `ggplot2`, `ape` and their dependencies. With local per-trial
results available, choose a task number in `scripts/view-trial.R` and source
the script to inspect its characters and reconstructed rates.

The manuscript files use `natbib` citations. `materials-and-methods.tex` and
`figure-captions.tex` are fragments for the paper; `methods-preview.tex` is a
standalone wrapper for checking the methods layout.

## Run the analysis

Install the dependencies listed in `provenance/ontophylo/DESCRIPTION` first.
The installer below installs the included package source into `Rlib/`.
Set `PHYLO_SAMPLING_RLIB` to use another library, and optionally set
`PHYLO_SAMPLING_DEPENDENCY_RLIB` for dependencies installed elsewhere.

```sh
Rscript scripts/install-package.R
Rscript scripts/00-prepare-design.R
for chunk in $(seq 1 19); do
  Rscript scripts/01-fit-node-bank.R "$chunk"
done
Rscript scripts/02-collect-node-bank.R
for task in $(seq 1 1500); do
  Rscript scripts/03-run-trial.R "$task"
done
Rscript scripts/07-audit-results.R
Rscript scripts/05-collect-results.R --final
Rscript scripts/06-plot-results.R
Rscript scripts/plot-design.R
Rscript scripts/08-report-results.R
Rscript scripts/09-check-summary.R
Rscript scripts/10-check-epoch-metrics.R
Rscript scripts/11-check-recent-signal.R
```

The serial commands show the execution order. The full analysis is intended
for a cluster. Existing trials are checked against their settings and reused.
The seeds and model settings are recorded in the scripts and saved design.

### Roihu

Submit from this directory after creating `logs/`. Supply your allocation with
`sbatch --account=YOUR_ALLOCATION`; no account or login details are stored here.
Jobs use the submission directory, or `PHYLO_SAMPLING_DIR` if set. The default
module is `r-env/461`, overridable with `PHYLO_SAMPLING_R_MODULE`.

Submit the stages below with `--dependency=afterok:JOB_ID` to enforce the
listed order. Check the six pilot trials before submitting production jobs.

| Stage | Script | Array |
| --- | --- | --- |
| Install package and prepare design | `slurm/setup.slurm` | — |
| Fit unique characters | `slurm/fit-bank.slurm` | `--array=1-19` |
| Collect fitted models | `slurm/collect-bank.slurm` | — |
| Map pilot trials | `slurm/pilot.slurm` | `--array=1-6` |
| Check pilot results | `slurm/check-pilot.slurm` | — |
| Map all trials | `slurm/batch.slurm` | `--array=1-100%100` |
| Audit and collect results | `slurm/collect-results.slurm` | — |

Each production array element runs the 15 conditions for one replicate.
After collection, run the plotting, reporting and summary checks listed above.

## Interpretation

These simulations impose single-origin clade patterns. They do not simulate
a continuous-time process with a known generating rate matrix. Selected node
ages define the focal epoch; averaging changes over branches, and merging
branches during pruning, can spread the reconstructed signal beyond it.

Character loss, temporal resolution and reconstruction error are distinct.
The two trees also have different lineage-time denominators. Empirical 95%
ranges across trials describe replicate variation; within-trial mapping
intervals condition on the fitted model and fixed tree. Recent-increase
recovery counts are descriptive, without a constant-rate null or calibrated
shift test. The experiment concerns the fixed study subset and does not
evaluate missing fossil lineages.
