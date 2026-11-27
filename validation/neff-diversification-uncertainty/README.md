# Neff+–diversification reconstruction uncertainty

This analysis propagates the variability in saved reconstructions into
correlations between Neff+ and diversification. PyRate comparisons use the
original trajectories; BAMM comparisons linearly detrend both trajectories.

## Reproduce

The script is at the repository root:

```sh
Rscript --vanilla make-neff-diversification-uncertainty.R
```

Set its initial working-directory line to your checkout location. The script
can also be sourced in RStudio. It requires `ape` and `BAMMtools`; the figure
font is Arial. The first run extracts the saved BAMM posterior. Later runs
reuse `binned-trajectories.rds` after checking input hashes. Set
`rebuild_inputs <- TRUE` if changing the metric definition or time bins.

## Outputs

- `extended-data-neff-diversification.pdf`: vector figure, 180 × 145 mm.
- Matching PNG and RGB, LZW-compressed TIFF files: 300-dpi artwork.
- `figure-caption.tex`: LaTeX caption; update the figure number in the manuscript.
- `correlation-summary.csv`: correlation quantiles and nominal significance frequencies.
- `figure-source-data.csv`: counts and percentages for every histogram bin.
- `binned-trajectories.rds`: input curves, draw metadata and input hashes.
- `pyrate-draw-index.csv`: source file, row and iteration of each PyRate draw.
- `validation.csv`: mean-curve reproduction, correlation and regression checks.
- `input-checksums.csv` and `session-info.txt`: input and software provenance.

## Calculation

Neff+ uses 1,000 saved stochastic-map replicates and 15 adult anatomical units.
Acceleration follows the empirical inclusive 95th-percentile rule within each
unit and map. Neff+ is the accelerated-branch frequency multiplied by the
exponential of the entropy of complete nonzero anatomical acceleration patterns.
The original 1,000 time slices and endpoint convention are retained, then each
replicate is averaged into 13 complete 20-Myr bins from 260 to 0 Ma.

PyRate uses the last 200 rows of each of the ten family-level
`BDCS-fam-ssing-20` logs. BAMM uses all 18,001 configurations in the saved
`bamm_data.RDA` on the 765-tip tree. Its 280-point rate grid matches
`make-BAMM/make-BAMM-data.R` and `data/bamm.rds`. BAMM is not pruned to the
phenomic sample. For both methods, extinction is subtracted from speciation
within the same posterior draw.

For BAMM, each Neff+ and diversification trajectory is regressed separately on
the 13 bin-midpoint ages. Their residual correlation is a partial Pearson
correlation controlling linear time, with two-sided p-values using 10 degrees
of freedom. The mean-curve reference lines use the same calculation. PyRate
uses undetrended trajectories and 11 degrees of freedom.

Each Neff+ trajectory is paired with every diversification trajectory, retaining
the full sequence of bins. This gives 2,000,000 PyRate and 18,001,000 BAMM
combinations with equal weights. Histograms are normalized separately to 100%.
Undefined correlations stop the script. The saved trajectories allow regeneration
of the pairwise results; full correlation matrices are not stored.

## Interpretation

The analysis propagates uncertainty represented by the saved reconstructions,
conditional on their trees, models, acceleration threshold and retained samples.
Pairing separately fitted distributions does not fit a joint model, and the
combinations are not independent observations or biological replicates.
Central 95% ranges describe reconstruction variability. The fraction of nominal
p-values below 0.05 is a sensitivity summary, not an integrated p-value or the
probability that the association is true. Linear detrending does not guarantee
removal of temporal autocorrelation, and this analysis does not establish causation.
