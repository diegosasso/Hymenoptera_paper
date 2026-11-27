# Paired Neff+ recovery figures

The plot labels **Correlated** and **Independent** correspond to Associated and
Unassociated in the source tables. Correlated describes the generating design;
it does not imply a significant correlation in every trial.

## Figures

- **Figure 1:** the example tree, known speciation schedule and all 100
  reconstructed Neff+ trajectories per scenario. Thin lines are trial means
  over conditional maps; the thick line averages trials. Both Neff+ panels use
  the same vertical scale. The example is prespecified trial 001.
- **Figure 2:** all 100 inferred correlations and nominal p-values per scenario.
  Horizontal jitter separates points without changing their values. Boxes
  summarize variation across trials. The p-value axis is logarithmic.

Both figures use reconstructed results and the primary strict acceleration
cutoff. PDFs contain vector artwork and embedded Arial text; PNG previews are
rendered at 450 dpi. `figure-captions.tex` contains the LaTeX captions and
inclusion code. `source-data/` contains plotted values, trial IDs, generating
rates, tree coordinates, input hashes and software information.

## Reproduce

From `validation/neff-recovery`, run:

```sh
Rscript --vanilla publication/make-publication-figures.R
```

The script reads existing results and writes the publication figures and source
data. It requires `ape`, `ggplot2`, `gridExtra`, Cairo graphics, Arial and
Poppler's `pdftoppm`. No inference is rerun.
