Neff+ correlations with family-level PyRate and full-tree BAMM

Run the separate script at the repository root:
  Rscript --vanilla make-neff-diversification-uncertainty.R
Or open it in RStudio and click Source. The first line sets the repository
working directory; edit that line if the checkout is elsewhere.
Required packages: ape and BAMMtools. Figure font: Arial.

The first run extracts the existing BAMM posterior and takes a few minutes.
Subsequent runs reuse binned-trajectories.rds after checking input hashes.
Set rebuild_inputs <- TRUE after changing the metric definition or time bins.
No existing input files, models, stochastic maps or other figures are modified.

Deliverables
  extended-data-neff-diversification.pdf: editable vector figure, 180 x 145 mm.
  extended-data-neff-diversification.png: 300-dpi preview.
  extended-data-neff-diversification.tiff: 300-dpi RGB, LZW-compressed artwork.
  figure-caption.tex: LaTeX caption; replace X with the assigned figure number.
  correlation-summary.csv: exact correlation quantiles and p-value percentages.
  figure-source-data.csv: counts and percentages for every displayed histogram bin.
  binned-trajectories.rds: all 13-bin input curves, draw metadata and input hashes.
  pyrate-draw-index.csv: source file, row and iteration for every PyRate draw.
  validation.csv: reproduction of original means and direct cor.test checks.
  input-checksums.csv and session-info.txt: input and software provenance.

Analysis
  Neff+: 1000 saved stochastic-map replicates, 15 adult anatomical units.
  Acceleration uses the empirical >= 95th-percentile rule, independently within
  each unit and map. Neff+ is frequency times exp(Shannon entropy) of complete
  nonzero anatomical acceleration patterns among branches present in a slice.
  The original 1000 time slices and endpoint convention are retained, then
  each replicate is averaged within 13 complete 20-Myr bins (260-0 Ma).

  PyRate: family-level BDCS-fam-ssing-20, final 200 rows from each of 10 logs.
  BAMM: all 18001 configurations in the saved bamm_data.RDA, on the full 765-tip
  tree. Its 280-point rate grid matches make-BAMM/make-BAMM-data.R and the
  existing data/bamm.rds summary. BAMM is not pruned to the phenomic sample.
  For both methods lambda and mu are subtracted within the same posterior draw.

  Each Neff trajectory is paired with every diversification trajectory.
  Complete time series are kept intact; bins are never independently resampled.
  There are 2000000 PyRate and 18001000 BAMM combinations, all equally weighted.
  No outcome filtering, multiple-testing correction or significance selection.
  Every finite pair contributes to the histograms, normalized separately to 100%.
  The script stops if any correlation is undefined rather than silently dropping it.
  Full pairwise matrices can be regenerated from the saved trajectories; only
  input curves and histogram summaries are stored to keep this repository small.

Interpretation
  This propagates the uncertainty represented by the saved reconstructions,
  conditional on the trees, models, acceleration threshold and retained samples.
  Pairing two separately fitted distributions does not fit a joint model.
  The combinations are not independent observations or biological experiments.
  The central 95% ranges describe reconstruction variability, not confidence
  intervals from cor.test. The fraction p<0.05 is a sensitivity summary, not
  an integrated p-value or the probability that an association is true.
  Pearson p-values use n=13 and df=11 and remain nominal; this analysis does
  not adjust for temporal autocorrelation or establish causation.
