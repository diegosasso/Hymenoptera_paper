PUBLICATION FIGURES: PAIRED NEFF RECOVERY

Plot labels are Correlated and Independent. These correspond to Associated
and Unassociated in the original analysis and source-data tables. Correlated
describes a design intended to generate a positive Neff-diversification
association; it does not guarantee a significant correlation in every trial.

Figure 1: example tree, known lambda schedule, and all 100 reconstructed
Neff trajectories per scenario. Thin lines are trial means over conditional
maps; the thick line is the mean across trials. The two Neff panels share
the same vertical scale. The example tree is the prespecified trial 001.

Figure 2: all 100 inferred correlations and p-values per scenario. Horizontal
jitter separates points without changing r or p. Boxes summarize the
distribution across trials, not uncertainty in any individual estimate.
The p-value axis is logarithmic. The strict acceleration cutoff is used
throughout, matching the primary comparison. No simulated-history results
or alternative cutoffs are mixed into these reconstructed-result figures.

FILES
figures/figure-1-neff-and-diversification.pdf        183 x 170 mm
figures/figure-2-correlations-and-p-values.pdf       183 x 112 mm
Matching PNG previews are rendered from the PDFs at 450 dpi.
figure-captions.tex provides captions and inclusion code for LaTeX.
source-data/ retains the plotted values, trial IDs, generating lambda,
example-tree coordinates, input hashes and plotting session information.

The PDFs contain embedded Arial text and vector artwork. The project uses
readable 9-10 pt labels and 11 pt panel letters; adapt these to a target
journal's specific requirements at final submission.

REPRODUCE
From validation/neff-recovery, or its existing RStudio project:
  source("publication/make-publication-figures.R")
Or:
  Rscript --vanilla publication/make-publication-figures.R

Requirements: ape, ggplot2, gridExtra, Cairo graphics, Arial and pdftoppm.
The script reads existing results and writes only publication/ outputs.
No simulation, model fitting, stochastic mapping or correlation testing
is rerun. The original scientific and figure files are preserved.
