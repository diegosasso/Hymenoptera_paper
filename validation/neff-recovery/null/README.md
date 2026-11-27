# Independent-character control

This control pairs each of the 100 associated trials with independent character
evolution on the same tree and known diversification trajectory. Run commands
from this directory, or open `null-recovery.Rproj`. The parent directory provides
the frozen ontophylo source and the paired trees.

## Design

Ten anatomical units each contain 20 independent binary characters. All
characters start in state 0 at the root and evolve on every branch under the
same ER model, q01 = q10 = 0.003 per Myr. There are no focal character origins,
rate pulses, tip-balance filters or character redraws. The rate was fixed before
examining these outcomes. The paired trees retain the conditioning used to
generate the associated scenario.

Independent character, fitting and mapping seeds are recorded in
`design/task-manifest.csv`, along with hashes of the paired inputs. Inference
uses ER, `maddfitz`, three optimization starts and 100 maps per variable
character. Only the tree and observed tip states enter reconstruction.
Constants are recorded and omitted; a unit with no variable characters has
unavailable reconstruction. Numerical failures do not cause resimulation.

Neff+ and the 15-bin Pearson comparisons follow the associated experiment.
The strict 95th-percentile cutoff is primary; inclusive ties are a sensitivity
analysis. P-values remain nominal. Two-sided significance and significance
with a positive correlation are reported separately. The resulting frequencies
measure false association under this generator, including tree geometry and
reconstruction effects.

## Outputs

`RESULTS.txt` summarizes the control and paired comparison. Important tables are:

- `summary/comparison-percentages.csv`: both scenarios, directions, thresholds,
  cutoff conventions, intervals and denominators.
- `summary/comparison-correlations.csv`: every trial-level correlation.
- `summary/metric-sources.csv`: simulated and reconstructed metric definitions.

Figures include scenario comparisons, a complete diagnostic atlas and all ten
units of trial 001. The diagnostic character-rate reference is the generating
expectation, 200q = 0.6; realized rates and reconstructed Neff+ can vary over time.
Combined publication figures are in `../publication/`.

`manuscript/` contains the paired methods and results, complementary supplementary
detail, captions and the compiled `methods-preview.pdf`. The preservation
manifest in `provenance/associated-inputs-sha256.txt` contains original hashes
for the retained associated data, code, manuscript and figure files.

## Run

Use the parent's installed library (`../Rlib`), or set `NEFF_RLIB`.
`NEFF_DEPENDENCY_RLIB` supplies additional dependency libraries if needed.

```sh
Rscript --vanilla scripts/01-simulate.R 1
Rscript --vanilla scripts/02-reconstruct.R 1
Rscript --vanilla scripts/03-summarize.R 1
Rscript --vanilla scripts/check-completion.R 1
Rscript --vanilla scripts/04-audit.R 1
```

Replace `1` with IDs 1–100, or use `run-one-trial.R` in RStudio. On Roihu,
submit from this directory after creating `logs/`, supplying your allocation
with `sbatch --account`. Trial 001 is the technical pilot and remains part of
the final 100 trials. Run the remaining inference jobs after its checks pass.

After all trials finish, run `scripts/check-completion.R`, `04-audit.R` and
`05-collect.R`, then:

```sh
Rscript --vanilla scripts/06-plot-batch.R
Rscript --vanilla scripts/07-plot-trials.R
Rscript --vanilla scripts/08-compare.R
python3 scripts/09-render-previews.py
Rscript --vanilla scripts/10-report.R
Rscript --vanilla scripts/11-diagnostics.R
latexmk -cd -pdf -interaction=nonstopmode -halt-on-error manuscript/methods-preview.tex
python3 scripts/12-verify-delivery.py
```

Cluster collection records checksums of portable outputs in
`provenance/portable-output-sha256.json`. Transfer the final `maps.rds` objects
and other per-trial outputs; redundant fit and map checkpoints are not required.
Raw `simulated/` and `results/` files remain local and are excluded from Git.
