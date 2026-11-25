# Empirical branch-label randomization

This test addresses whether the placement of accelerated anatomical units on the
empirical tree produces a stronger Neff+–diversification association than expected
from random branch assignments. It uses the original 77-tip, 152-branch tree and
15 adult anatomical units, excluding the larval unit. Existing inputs and earlier
simulation studies are unchanged.

The observed raw correlation was **r = 0.9501**. Neither null produced a correlation
this large in 9,999 permutations (empirical upper-tail P = 0.0001 for each). However,
the length-constrained null had median r = 0.6805 and 86.13% of its permutations
had nominal Pearson P < 0.05. Thus the observed alignment exceeds both specified
nulls, while a positive association is common when the empirical relationship
between branch length and acceleration is retained. Full numbers are in
`results/RESULTS.txt` and `results/permutation-summary.csv`.

## Design

For each unit and each of its 1,000 saved reconstruction samples, a branch is
accelerated when its rate is **greater than or equal to** that sample's 95th
percentile. This reproduces the empirical manuscript definition, including ties.
There is no new ancestral-state reconstruction, character simulation, tree
simulation or PyRate fitting.

Each randomization independently permutes the 152 branch labels within each unit.
The same permutation moves that unit's complete 1,000-sample branch profile, so
every sample retains its original accelerated-branch count and the reconstruction
uncertainty within each unit is preserved. The resulting assignments remain fixed
along entire branches and across all time slices. Units receive independent
permutations; their original cross-unit co-occurrence is not held fixed.

Two null models were specified before the full run:

1. **Unrestricted:** exchange profiles among all branches, as requested.
2. **Within length quartiles:** exchange profiles only among branches in the same
   quartile of the empirical branch-length distribution. The four groups contain
   38, 39, 37 and 38 branches, respectively, because tied lengths are kept together.
   This sensitivity analysis retains the association between branch
   length and the frequency of accelerated assignments at this coarse resolution.

There are 9,999 permutations per model, using seed 20261123. They are branch-label
randomizations of an existing reconstruction, not independent simulated trees.

## Calculation

At each slice, Neff+ is `f_plus * exp(H_plus)`, where `f_plus` is the fraction of
intersecting branches with at least one accelerated unit and `H_plus` is Shannon
entropy of the frequencies of distinct nonzero **whole 15-unit patterns** among
these branches. It is zero when no branch is accelerated. Neff+ is calculated
separately for each reconstruction sample and then averaged; the nonlinear metric
is not calculated from a mean-rate mask.

All calculations retain the original 1,000-slice grid, branch-intersection rule,
endpoint convention and age assignment. Slices sharing the same set of branches
are evaluated once for speed. The primary statistic is Pearson's r between mean
Neff+ and the fixed family-level PyRate mean **net** diversification trajectory
(lambda minus mu), averaged into the same thirteen complete 20-Myr bins from
260 to 0 Ma as the manuscript. PyRate values are not shuffled or refitted.

For the prespecified positive association, the empirical permutation P value is
`(1 + number of shuffled r >= observed r) / 10000`. The table additionally gives
the tail based on absolute r, and secondary correlations after linear detrending
and first differencing. Nominal two-sided Pearson P values are reported separately;
the residual-correlation convention matches the manuscript. These nominal values
are not the basis for the randomization inference. The empirical test compares
whole trajectories and does not treat the thirteen bins as independent draws.

## Validation and interpretation

The reconstructed mean Neff+ curve must match the existing saved curve to within
1e-10 before any permutation runs. The compiled entropy calculation is also
checked against a plain-R calculation. Branch-count preservation is checked for
every sample and unit in the first permutation of each model, including within
length quartiles for the constrained model. Input checksums, R session details,
the observed bin values and every null correlation are saved.

The input audit retains an existing threshold edge case: in 27 of the 1,000
female-genitalia samples, the 95th percentile is zero and the inclusive threshold
labels all 152 branches as accelerated. This is reproduced in observed and shuffled
data rather than silently changing the empirical statistic.

The unrestricted null asks whether random placement with the observed per-unit
counts and fixed tree can reproduce the association. It does not preserve the
strong dependence of accelerated assignments on branch length. The constrained
null addresses that dependence approximately. Neither null preserves all possible
clade structure, within-unit phylogenetic clustering or character-model detection
biases. A significant result demonstrates unusually strong temporal alignment under
the stated randomization, and does not establish a causal effect of phenomic
diversification on lineage diversification. Any positive correlation common under
the constrained null must also be reported; rejection of that null does not imply
that branch geometry contributes nothing to the observed association.

## Run

Open `run-shuffles.R` in RStudio and set the working directory to this folder (or
the repository root), then source it. After it finishes, source `make-plots.R`.
From the repository root, the equivalent commands are:

```sh
Rscript --vanilla validation/neff-recovery/empirical-branch-shuffle/run-shuffles.R
Rscript --vanilla validation/neff-recovery/empirical-branch-shuffle/make-plots.R
```

Required R packages are ape, Rcpp, ggplot2 and gridExtra. Rcpp needs a working C++
compiler. `neff-kernel.cpp` only accelerates entropy evaluation; all randomization
and analysis choices are visible in R. The code also includes the plain-R metric
calculation used to verify the helper.

An optional quick timing check can use `NEFF_SHUFFLES=20` and a separate
`NEFF_OUTPUT_DIR`. The full results must use the default 9,999 shuffles.

## Outputs

- `results/permutation-summary.csv`: observed r, null intervals, exceedance counts,
  empirical P values and percentages of nominal Pearson P values below 0.05.
- `results/null-correlations.csv`: every permutation, both nulls and all statistics.
- `results/observed-bins.csv`: the thirteen paired empirical observations and
  sampled-tree branch counts.
- `results/analysis.rds`: reproducible numeric results, including all randomized
  group-wise and binned mean Neff+ curves.
- `results/accelerated-branch-counts.csv` and `results/branch-length-audit.csv`:
  threshold/count and branch-duration diagnostics.
- `results/figures/branch-shuffle-test.pdf`: randomized Neff+ curves and primary
  correlation distributions; PNG preview and source data are also provided.
- `results/figures/empirical-trajectories.pdf`: Neff+, net diversification and
  intersecting-branch counts on separate, aligned axes.
- `results/figure-captions.txt`: captions defining intervals, statistics and axes.
- `results/RESULTS.txt`: numerical findings and their interpretation.

Existing manuscript LaTeX and publication figures are not edited by these scripts.
