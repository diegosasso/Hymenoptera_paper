# Neff+ branch-shuffling publication figure

This directory contains one publication figure using **all 9,999 unrestricted
branch-label randomizations**. No length-constrained results are included. The
earlier exploratory analysis is preserved in its original directory.

The supported conclusion is: **random placement of accelerated anatomical units,
conditional on the empirical tree and per-unit counts, does not reproduce the
observed positive Neff+–diversification association.** This does not imply that
shuffling creates no correlation, or that the empirical association is causal.

## Figure

`figures/neff-shuffling.pdf` is the vector master (183 x 170 mm, embedded Arial).
`figures/neff-shuffling.png` is its 450-dpi preview.

- **a:** all 9,999 shuffled Neff+ trajectories, their arithmetic mean, and the
  observed trajectory on the full original time grid. Every trajectory first
  averages Neff+ across the 1,000 reconstruction samples. All shuffles are drawn;
  no subsampling, smoothing or filtering is used. Redundant vertices on constant
  segments are removed solely to keep the vector file smaller: this preserves the
  line geometry of the original 1,000 time points exactly.
- **b:** the distribution of Pearson r across the thirteen complete 20-Myr bins
  from 260 to 0 Ma, compared with the observed r. The annotation gives the
  positive-tail empirical permutation p value, including the plus-one correction.
- **c:** the distribution of nominal two-sided Pearson p values, shown as
  -log10(p), with the p = 0.05 threshold and observed nominal p marked. Counts of
  significant negative and positive correlations are distinguished explicitly.

`figure-caption.tex` and `methods.tex` are LaTeX snippets for the manuscript.
`materials-and-methods.tex` is the main-text section covering recovery, threshold
sensitivity, component analyses and randomization, explaining the purpose and
interpretation of each test. The file includes a standalone preview wrapper;
copy the subsection into the manuscript and update the Supplementary Note numbers.

`supplementary-worked-example.tex` explains how to calculate Neff+ and Pc from
the four scenarios in Extended Data Figure 4. It includes a small comparison table
and one modified example with an unaccelerated lineage. Run
`Rscript validation/neff-shuffling/scripts/check-worked-example.R` from the
repository root to reproduce all values and verify them against the existing
metric functions; results are in `source-data/worked-example-values.csv`.

## Results

The observed correlation is r = 0.950066, with nominal two-sided Pearson
p = 6.70693e-7. The randomized correlations have median -0.474269 and central 95%
interval [-0.571086, -0.296410]; their maximum is -0.064637. None reaches the
observed positive correlation, giving empirical permutation p = 0.0001.

There are 658/9,999 (6.5807%) nominal p values below 0.05. All correspond to
**negative** correlations. There are zero significant positive correlations.
Nominal p values assume independent bin observations and are included as a
descriptive comparison; the empirical test compares complete shuffled trajectories
and does not rely on independence of the thirteen bins.

## Recreate the figure

Open `make-figure.R` in RStudio, set the working directory to this folder, and
source the script. It reads only data and style files inside this directory and
does not rerun simulations. Alternatively, from the repository root:

```sh
Rscript --vanilla validation/neff-shuffling/make-figure.R
```

R packages required for plotting are ggplot2 and gridExtra; grid is part of R.
Arial is the project font. The figure's source data and statistics are saved in
`source-data/`; all computed correlations and p values can be inspected directly.

## Recreate the randomizations

`run-shuffling.R` is a readable, unrestricted-only R script using the empirical
inputs already in this repository. The optional C++ calculation helper is included
in `scripts/neff-kernel.cpp` and is checked against the plain-R definition before
simulation. Run from this folder or the repository root:

```sh
Rscript --vanilla validation/neff-shuffling/run-shuffling.R
Rscript --vanilla validation/neff-shuffling/make-figure.R
```

The defaults are 9,999 randomizations and seed 20261123; a full rerun takes about
five minutes on the machine used here. This additionally requires ape, Rcpp and a
C++ compiler. A small verification run can use `NEFF_SHUFFLES=20` and a separate
`NEFF_OUTPUT_DIR`; this prevents replacing the publication data with a pilot.

The supplied publication data were imported from the existing completed
unrestricted run using `prepare-source-data.R`, without repeating or selecting
randomizations. That import script is retained for provenance; it is not needed
to remake the figure or to rerun the analysis. Import checksums are recorded in
`provenance/`.

## Data and calculation

`source-data/figure-data.rds` contains every shuffled curve in a compact format:

- `slice_age`, `slice_group`: the 1,000-point age grid and its 75 branch-membership
  groups. Slices with identical intersecting branches have identical Neff+.
- `shuffle_neff`: a 75 x 9,999 matrix of mean Neff+, one column per randomization.
  The full 1,000-point curve for shuffle j is
  `shuffle_neff[slice_group, j]`.
- `observed_neff`: the observed mean curve on the same 1,000-point grid.
- `observed_bins`, `shuffle_bins`: empirical values and the 13 x 9,999 matrix
  used in the correlations. `bin_weights` reproduces temporal averaging.
- unit names, number of maps, tip/branch counts, number of shuffles and seed.

The other source files give every correlation and nominal p value, summary
statistics, the observed bin data, the full-grid mean trajectories and original
empirical-input checksums. `provenance/figure-validation.csv` records that every
shuffle is drawn and that vertex compression preserves the curves.

Accelerated assignments use the original **rate >= 95th percentile** rule,
separately for each unit and reconstruction sample. The complete 1,000-sample
profile of a branch is permuted as a block within a unit; units receive independent
permutations. Counts are therefore preserved within every sample and unit, and
each assignment remains fixed along its branch across time slices. Neff+ is
recalculated from the joint 15-unit binary patterns, not from mean branch rates.
The family-level PyRate mean net diversification trajectory (lambda minus mu)
remains fixed throughout.

For reproducibility, the original threshold convention is retained even when the
95th percentile is zero: this occurs in 27 female-genitalia samples and labels all
152 branches as accelerated in those samples. It affects observed and randomized
calculations equally. No original input or earlier analysis has been edited.
