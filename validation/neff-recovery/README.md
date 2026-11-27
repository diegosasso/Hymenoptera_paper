# Neff+ recovery simulations

This experiment evaluates reconstruction of a phenomic signal associated with
a known diversification history. It contains 100 trials, each with an independently
simulated tree and 200 binary characters. The paired independent-character
control is described in [null/README.md](null/README.md); publication figures for
both scenarios are in [publication/README.md](publication/README.md).

## Design and inference

Trees have a crown age of 300 Ma and 120–250 tips. Speciation rates are 0.010,
0.028 and 0.010 per Myr at 300–180, 180–80 and 80–0 Ma, respectively; extinction
is zero. Ten overlapping focal branches define ten anatomical units. Each unit
contains 20 binary characters, originating in state 1 on its focal branch and
subsequently evolving under slow ER evolution (q01 = q10 = 0.0001 per Myr).
`scripts/01-simulate.R` gives the geometry and conditioning rules.

Each variable character is fitted with ER and `maddfitz` root treatment, using
three optimization starts and 100 conditional stochastic maps. Constant
observations are recorded and omitted. Full and retained-character histories
are saved separately. Tree, character, fitting and mapping seeds are fixed;
failed inference does not trigger replacement of a biological replicate.

Unit rates sum character changes and divide by branch duration. The primary
Neff+ definition uses branch rates strictly above each unit's 95th percentile,
separately within each map, and combines the joint acceleration patterns across
units. The inclusive cutoff is retained as a sensitivity analysis.

Diversification is known from the generator. Pearson tests use 15 matched
20-Myr bins; their p-values are nominal and do not correct temporal
autocorrelation. Recovery frequencies in this scenario describe the imposed
association. The paired control evaluates false associations under its
independent-character generator. See `RESULTS.txt` for numerical results.

## Files

- `design/`: the 100-trial seed manifest.
- `scripts/` and `slurm/`: simulation, reconstruction, checks and cluster jobs.
- `summary/`: correlations, recovery measures, denominators and metric provenance.
- `figures/`: batch summaries, all-trial diagnostic atlases and trial 001 examples.
- `manuscript/`: LaTeX methods, captions, references and a compiled preview.
- `provenance/`: frozen ontophylo source, versions, checksums and numerical audits.
- `null/`: the paired independent-character control.
- `publication/`: combined figures and their source data.
- `empirical-branch-shuffle/`: a separate randomization of the observed branch profiles.

The local `simulated/` and `results/` directories are excluded from Git.
Trial 001 is the prespecified detailed example. Within-trial bands describe
conditional mapping variation; batch bands describe variation among trial means.

## Run

Open `neff-recovery.Rproj`, or run commands from this directory. Install the
dependencies in `provenance/ontophylo/DESCRIPTION`, plus `diversitree`.
Plotting uses `ggplot2`, `gridExtra`, Cairo graphics, Arial and Poppler.

```sh
Rscript --vanilla scripts/install-package.R
Rscript --vanilla scripts/01-simulate.R 1
Rscript --vanilla scripts/02-reconstruct.R 1
Rscript --vanilla scripts/03-summarize.R 1
```

Replace `1` with trial IDs 1–100. Alternatively, use `run-one-trial.R` in
RStudio. The frozen package installs into `Rlib/`, or `NEFF_RLIB` if set.
`NEFF_DEPENDENCY_RLIB` accepts additional libraries separated by the platform
path separator. Reconstruction uses `SLURM_CPUS_PER_TASK` workers, or one
worker locally. Matching checkpoints resume with the saved settings and seeds.
Use a separate output directory when changing the design or inference settings.

On Roihu, create `logs/` and supply your allocation with `sbatch --account`.
Jobs use the submission directory, or `NEFF_DIR`, and the `r-env/461` module,
or `NEFF_R_MODULE`. Run setup and simulation, reconstruct trial 001, and check
it with `slurm/check-pilot.slurm` before the remaining trials. Package setup
uses 32 GB; reconstruction uses eight CPUs and 16 GB per trial.

After all trials finish:

```sh
Rscript --vanilla scripts/check-completion.R
Rscript --vanilla scripts/04-audit.R
Rscript --vanilla scripts/05-collect.R
Rscript --vanilla scripts/06-plot-batch.R
Rscript --vanilla scripts/07-plot-trials.R
Rscript --vanilla scripts/08-report.R
python3 scripts/09-render-previews.py
Rscript --vanilla scripts/11-additional-summaries.R
latexmk -cd -pdf manuscript/methods-preview.tex
python3 scripts/10-verify-delivery.py
```

The `*-standalone.tex` manuscript fragments contain expanded numerical results;
other fragments use `results-macros.tex`. Include the accompanying references
and update figure cross-references when incorporating them into the paper.
