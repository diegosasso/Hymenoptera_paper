HOMOGENEOUS ROOT ER CONTROL

Paired no-association scenario for all 100 completed associated trials.
Run commands from this null/ directory; the parent holds the preserved
associated inputs, frozen OntoPhylo source and installed R library.

PRESPECIFIED DESIGN
Use exactly the corresponding associated tree and its known diversification
trajectory. Ten anatomical units each contain 20 independent binary characters.
All characters start in state 0 at the root and evolve on every branch under
the same ER Q, q01=q10=0.003 per Myr. This is the rate in the earlier root ER
pilot. It was fixed before inspecting these 100 control outcomes.
There are no focal branches, imposed origins, pulse in Q, tip-balance filters,
or character redraws. Retaining the same trees also retains their original
size and focal-geometry conditioning; this is a null conditional on that tree
ensemble, not an unconditional sample of all birth-process trees.

The 100 independent character, fitting and mapping seed triples are recorded
in design/task-manifest.csv alongside hashes of the paired source inputs.
Inference uses ER, maddfitz, 3 restarts and 100 maps per variable character.
Only tree and terminal observations enter ASR. The known Q, true histories,
root state and diversification schedule are not supplied to inference.
Constants are recorded and omitted, without fabricated zero-change maps.
A unit without any variable character makes reconstruction unavailable.
Numerical failures are reported, never replaced by new simulated data.

Neff calculations and 15 matched 20-Myr Pearson tests are the same as in the
associated experiment. Strict 95th-percentile acceleration is primary;
inclusive ties are a sensitivity analysis. P-values remain nominal. Report
two-sided significance of either sign and positive-r significance separately.
The resulting null frequencies quantify false association under this particular
independent-character generator, including tree geometry and reconstruction.

RUNNING
Open null-recovery.Rproj in RStudio or work in this directory.
Use the parent's isolated frozen library (default ../Rlib), or NEFF_RLIB.
Set NEFF_DEPENDENCY_RLIB for additional dependency libraries if required.

Rscript --vanilla scripts/01-simulate.R 1
Rscript --vanilla scripts/02-reconstruct.R 1
Rscript --vanilla scripts/03-summarize.R 1
Rscript --vanilla scripts/check-completion.R 1
Rscript --vanilla scripts/04-audit.R 1

The Roihu scripts use the same modules/resources as the associated workflow.
Trial 001 is a technical preflight included in the final 100. All simulations
are drawn before their significance is inspected. The remaining 99 inference
jobs follow a successful preflight. A technical check, never a p-value, gates
production. After all 100 finish, run scripts/check-completion.R,
04-audit.R, 05-collect.R, then the plotting and comparison scripts.

All control outputs stay here. Parent scientific files are preserved and
verified against provenance/associated-before-sha256.txt.

FINAL OUTPUTS
RESULTS.txt summarizes control frequencies and the paired benchmark.
summary/comparison-percentages.csv contains both scenarios, two directions,
three p-value thresholds, both cutoffs, exact intervals and denominators.
summary/comparison-correlations.csv contains every trial-level result.
summary/metric-sources.csv distinguishes simulated and reconstructed metrics.

figures/scenario-p-value-percentages.pdf
figures/scenario-p-value-distribution.pdf
figures/scenario-neff-through-time.pdf
figures/all-trial-diagnostics.pdf          (all 100 control trials)
figures/trial-001-unit-acceleration.pdf    (all ten anatomical units)

The full diagnostic atlas shows trees coloured by generating lambda, known
marginal diversification, marginal character rates and Neff. The dotted
character-rate reference is the constant generating expectation 200q=0.6;
realized branch rates and reconstructed Neff need not be constant in time.

manuscript/main-text-methods-standalone.tex is the combined methods insert.
manuscript/main-text-results-standalone.tex reports the paired results.
manuscript/supplementary-methods-standalone.tex gives complementary detail.
manuscript/methods-preview.pdf contains the combined text, table and figures.
The original associated-only documents in the parent folder are preserved.

POSTPROCESSING (after all 100 finish)
Rscript --vanilla scripts/06-plot-batch.R
Rscript --vanilla scripts/07-plot-trials.R
Rscript --vanilla scripts/08-compare.R
python3 scripts/09-render-previews.py
Rscript --vanilla scripts/10-report.R
Rscript --vanilla scripts/11-diagnostics.R
latexmk -cd -pdf -interaction=nonstopmode -halt-on-error manuscript/methods-preview.tex

Cluster collection creates provenance/portable-output-sha256.json. Transfer
the outputs without redundant fit-checkpoints/ or map-checkpoints/ folders;
final maps.rds objects contain portable branch counts for every map.
Run python3 scripts/12-verify-delivery.py after transfer and rendering.
