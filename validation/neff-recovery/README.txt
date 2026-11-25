NEFF RECOVERY SIMULATIONS

This folder contains the 100-trial recovery experiment based on the final
forced-origin / slow-evolution pilot. Run all R commands from this folder,
or open neff-recovery.Rproj in RStudio and use run-one-trial.R.
See RESULTS.txt for numerical results.

DESIGN

100 independently simulated crown trees, age 300 Ma, 120-250 tips.
Speciation: 0.010 / 0.028 / 0.010 per Myr at 300-180 / 180-80 / 80-0 Ma.
Extinction is zero. Elevated speciation applies to every lineage.
Ten overlapping focal branches define ten anatomical units. Each unit has
20 binary characters, each born in state 1 on its focal branch, followed by
slow ER evolution (q01=q10=0.0001 per Myr) in its descendants.
The complete geometry and conditioning rules are in scripts/01-simulate.R.
Tree and character seeds were fixed before examining any batch results.

INFERENCE AND INTERPRETATION

Each variable character is fitted with ER and maddfitz root treatment, using
three random optimization restarts and 100 conditional stochastic maps.
Constant observations are recorded and omitted from reconstruction; they
are not assigned invented zero-change histories or redrawn. Full and
retained-character truth are both saved. Failed reconstructions are not
replaced by a new biological replicate.

The primary Neff+ definition uses branch rates strictly above each unit's
95th percentile, separately within each map. The original inclusive cutoff
is retained as a sensitivity analysis; constant trajectories have undefined
correlations. Unit rates are summed character changes / branch duration,
without dividing by character count. Neff uses joint patterns across units.

Diversification is KNOWN from the generator. It is not estimated from the
tree, and make.bisse.td() is a likelihood check, not a forward simulator.
Pearson tests use 15 complete matched 20-Myr bins. P-values are nominal and
do not correct temporal autocorrelation. Percentages under this imposed
association are not null-calibrated statistical power or evidence of causality.
There is no independent-character null scenario in this batch.

CONTENTS

  design/       Prespecified 100-trial seed manifest.
  scripts/      Explicit R analysis, validation, plotting and reporting scripts.
  slurm/        Portable Roihu setup, simulation, reconstruction and audit jobs.
  simulated/    Trees, 200-character matrices, exact histories, event times and
                selection diagnostics for every trial (local, excluded from Git).
  results/      Fitted models, portable per-character branch-count maps, Neff
                samples, rates and status tables (local, excluded from Git).
  summary/      Complete tables, p-value percentages, recovery measures and
                metric-source information. Percentages include all denominators.
  figures/      Vector PDFs and PNG previews, including every trial.
  manuscript/   Main-text and complementary supplementary LaTeX, captions,
                references, numerical inserts and compiled preview.
  provenance/   Frozen installable OntoPhylo source, checksums, audit results,
                version/session information and job record.
  logs/         Execution logs (local, excluded from Git).

MAIN OUTPUTS

  figures/recovery-summary.pdf
  figures/neff-through-time.pdf
  figures/p-value-percentages.pdf
  figures/p-value-distribution.pdf
  figures/all-trial-neff-curves.pdf       (all 100 trials, ten per page)
  figures/all-trial-diagnostics.pdf       (one full diagnostic page per trial)
  figures/trial-001-overview.pdf
  figures/trial-001-unit-acceleration.pdf (all ten units of trial 001)
  summary/p-value-percentages.csv
  manuscript/methods-preview.pdf

Trial 001 was selected in advance for detailed unit plots. No illustrative
trial was chosen because of its p-value. Within-trial bands describe mapping
variation; batch bands describe variation among trial means.

REPRODUCING THE ANALYSIS

1. Provide R and the dependencies listed in provenance/ontophylo/DESCRIPTION,
   plus diversitree. For plotting, ggplot2, gridExtra, Cairo/Arial and Poppler pdftoppm are used.
   Set NEFF_DEPENDENCY_RLIB to optional dependency-library directories
   separated by the platform path separator. The frozen OntoPhylo itself
   is installed into Rlib, or the directory given by NEFF_RLIB.

   Rscript --vanilla scripts/install-package.R

2. Run each trial once (replace 1 with IDs 1-100):

   Rscript --vanilla scripts/01-simulate.R 1
   Rscript --vanilla scripts/02-reconstruct.R 1
   Rscript --vanilla scripts/03-summarize.R 1

   Reconstruction uses SLURM_CPUS_PER_TASK workers, or one worker locally.
   Existing simulation and inference files are retained. Matching fit/map
   checkpoints can resume an interrupted execution with identical seeds.
   Do not change settings in an existing output directory. To repeat the
   experiment from scratch, copy the scripts, design and frozen source to
   a new directory. Never redraw trials until significant.

3. On Roihu, submit setup and simulation jobs, then reconstruct trial 001,
   check it with slurm/check-pilot.slurm, and submit trials 002-100 only if
   the technical audit succeeds. Pass your account with sbatch --account,
   and set dependency-library paths via --export. No account or SSH key is
   embedded in the scripts. Production used eight CPUs and 16 GB per trial.
   Package setup needs 32 GB for the bundled OntoPhylo datasets.

4. After all 100 trials finish:

   Rscript --vanilla scripts/check-completion.R
   Rscript --vanilla scripts/04-audit.R
   Rscript --vanilla scripts/05-collect.R
   Rscript --vanilla scripts/06-plot-batch.R
   Rscript --vanilla scripts/07-plot-trials.R
   Rscript --vanilla scripts/08-report.R
   python3 scripts/09-render-previews.py
   Rscript --vanilla scripts/11-additional-summaries.R
   # After compiling the methods preview:
   python3 scripts/10-verify-delivery.py

5. Compile manuscript/methods-preview.tex with latexmk -pdf from manuscript/.
   The *-standalone.tex snippets contain literal results for easy Overleaf
   insertion; the source snippets use results-macros.tex, generated from CSVs.
   Include references.bib, captions and figures, and update figure numbering
   through LaTeX cross-references in the final manuscript.

The original pilot and the Taxon-sampling-sim / validation/phylo-sampling
workflow were consulted. Their files were not modified. This experiment
uses a frozen development source snapshot, not an unrecorded live checkout.
