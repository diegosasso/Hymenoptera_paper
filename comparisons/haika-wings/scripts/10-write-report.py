"""Write the numerical interpretation from verified comparison tables."""
import csv
from math import isclose
from pathlib import Path


def read_table(path):
    with open(path, newline="") as stream:
        return list(csv.DictReader(stream))


def select(rows, **conditions):
    selected = [row for row in rows if all(row[key] == value for key, value in conditions.items())]
    if len(selected) != 1:
        raise ValueError((conditions, len(selected)))
    return selected[0]


associations = read_table("results/association-summary.csv")
null = read_table("results/null-association-summary.csv")
bounds = read_table("results/bound-sensitivity.csv")
tips = read_table("data/tree-tip-metadata.csv")
case3_scores = read_table("data/haika-case3-published-pc-scores.csv")
fit_status = read_table("results/fit-status.csv")
null_status = read_table("results/null-status.csv")
assert len(fit_status) == 51 and all(row["comparison_complete"] == "TRUE" for row in fit_status)
assert len(null_status) == 199 and all(row["success"] == "TRUE" for row in null_status)
assert len(bounds) == 8 and all(row["success"] == "TRUE" for row in bounds)

summary_rows = []
for dataset in ["seed3407_full77", "seed3407_family48_first", "seed3407_family48_last"]:
    for scope in ["all_edges", "exclude_root_edges", "interfamily_edges", "time_full_width_bins"]:
        row = select(associations, dataset=dataset, axis="TOTAL_VARIANCE", **{
            "lambda": "1", "ontophylo_partition": "fore_wing", "scope": scope})
        summary = {key: row[key] for key in ["dataset", "scope", "n_units", "spearman_rho"]}
        summary["p_positive_concordance"] = ""
        summary["p_two_sided_tail"] = ""
        if dataset != "seed3407_family48_last":
            null_row = select(null, dataset=dataset, axis="TOTAL_VARIANCE", partition="fore_wing", scope=scope)
            assert isclose(float(row["spearman_rho"]), float(null_row["statistic"]), abs_tol=1e-12)
            summary["p_positive_concordance"] = null_row["p_positive_concordance"]
            summary["p_two_sided_tail"] = null_row["p_two_sided_tail"]
        summary_rows.append(summary)
with open("results/primary-summary.csv", "w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(summary_rows[0]))
    writer.writeheader()
    writer.writerows(summary_rows)

table = ["Alignment                  Branch rho    Interfamily rho    Time-bin rho"]
for dataset, label in [("seed3407_full77", "Original 77 tips"),
                       ("seed3407_family48_first", "48 families, first tip"),
                       ("seed3407_family48_last", "48 families, last tip")]:
    values = []
    for scope in ["all_edges", "interfamily_edges", "time_full_width_bins"]:
        row = select(summary_rows, dataset=dataset, scope=scope)
        values.append(float(row["spearman_rho"]))
    table.append(f"{label:27s}{values[0]:9.3f}{values[1]:19.3f}{values[2]:16.3f}")

case3_labels = {row["taxon_label"] for row in case3_scores}
case3_genera = {label.split("_")[0] for label in case3_labels}
morphology_genera = {row["Genus_Sharkey"] for row in tips}
molecular_genera = {row["Genus_Blaimer"] for row in tips}
exact_species = {row["molecular_species"].replace(" ", "_") for row in tips} & case3_labels
case3_overlap_text = (
    f"The published Case 3 scores have {len(exact_species)} exact molecular-species matches.\n"
    f"Genus overlap is {len(morphology_genera & case3_genera)} with the morphology labels, "
    f"{len(molecular_genera & case3_genera)} with the molecular labels,\n"
    f"and {len((morphology_genera | molecular_genera) & case3_genera)} distinct genera across "
    "their union. These are name overlaps, not verified\n"
    "independent images for the study exemplars. More fundamentally, some released\n"
)
strong77 = select(bounds, dataset="seed3407_full77", scenario="strong_smoothing_root_age")
strong48 = select(bounds, dataset="seed3407_family48_first", scenario="strong_smoothing_root_age")
smoothing_text = (
    "Because branch lengths are in Myr, penalty 1 allows extremely flexible log-rate\n"
    "variation over a roughly 280-Myr tree. An additional diagnostic penalty equal\n"
    "to the root age corresponds to approximately unit root-to-tip log-rate standard\n"
    "deviation on the penalty's scale. It removes all lower-bound hits, but branch-rank\n"
    "agreement with the primary model falls to "
    f"{float(strong77['rate_agreement_default_all_edges']):.3f} (77 tips) and "
    f"{float(strong48['rate_agreement_default_all_edges']):.3f} (48 tips).\n"
    "Forewing correlations across all "
    f"{strong77['n_edges']}/{strong48['n_edges']} branches become "
    f"{float(strong77['forewing_rho_all_edges']):.3f}/"
    f"{float(strong48['forewing_rho_all_edges']):.3f}. The\n"
    "sensitivity plots exclude the two basal branches from branch correlations;\n"
    "their corresponding values across "
    f"{int(strong77['n_edges']) - 2}/{int(strong48['n_edges']) - 2} branches are "
    f"{float(strong77['forewing_rho_exclude_root_edges']):.3f}/"
    f"{float(strong48['forewing_rho_exclude_root_edges']):.3f}. Both\n"
    "scopes use the same mean individually normalized ontophylo maps. Weak concordance\n"
    "persists, but detailed branch rankings are smoothing-sensitive. This diagnostic\n"
    "was prompted by boundary behavior; it was not selected to maximize agreement\n"
    "with ontophylo. The penalty-1 null probabilities do not apply to this diagnostic.\n"
)

report = """HAIKA and ontophylo wing evolutionary rates
Completed exploratory comparison, 20 November 2026

The original 77-tip dated tree can be retained by explicitly assigning HAIKA
family-level image profiles to its tips. Under this alignment, the primary
HAIKA and ontophylo forewing branch rates show weak concordance, and there is
no reproducible agreement in their marginal rate-through-time patterns across
all examined ANN initializations and sampling choices. This is a conditional
comparison of representations and reconstructions, not a species-matched
validation experiment or evidence that either method is intrinsically invalid.

The requested alignment
-----------------------
All 77 tip positions and all 152 branches of the original tree are retained.
The HAIKA assignment uses literal source-image family labels: 48 historical
Sharkey morphology families, represented by 112 source images. Forty-one
families have one image. The full crosswalk identifies both the molecular
species used to date/place each tip and its ontophylo morphology exemplar.
There are 22 differing morphology/molecular genus labels. Four tips also have
differing morphology/molecular family labels, which are explicitly recorded.

The 77 assignments contain 48 distinct family profiles. Nineteen families have
multiple study tips, producing 29 excess assignments; 50 branches lie wholly
within a single family with identical assigned HAIKA phenotypes. Low rates on
those branches may be imposed by proxy assignment rather than measured stasis.
Two 48-tip sensitivity trees retain one alphabetically first or last study
representative per family, conserving original ages and branch paths.
Pteromalidae, Cynipidae and Figitidae are nonmonophyletic here, making this a
representative-selection sensitivity rather than a collapse into family clades.

Data audit and HAIKA inputs
--------------------------
""" + case3_overlap_text + """Case 3 species vectors are numerically identical to differently labelled Case 2
image vectors: Tetragonula carbonaria matches a Trigonalidae row across all 28
columns, and Osmia/Heriades rows match an Osiris row. Standard PCA of that released
matrix does not reproduce the published species PC scores. The cause remains
unresolved. The S4 genus-labelled matrix also repeats identical vectors within
families. These products were retained for audit. The fitted comparison uses
the explicit family-profile assignment described above.

The fitted inputs are the original image-labelled 64-dimensional embeddings
from five archived ANN initializations. Seed 3407 is the prespecified primary
configuration; the other four seeds and the released 28-column image matrix
are representation sensitivities. Image vectors are averaged within each
matched family, then centered, unscaled PCA is fitted with each family weighted
once. At least two axes and at least 90% of variance are retained. Seed 3407
retains three axes explaining 96.94%. These newly calculated axes are not the
publication's named wing-shape and venation axes. All downloaded archive files
used in extraction passed their published checksums.

What the rates measure
----------------------
Ontophylo forewing rates sum mapped state changes in nine characters and divide
by branch duration, using 1,000 existing stochastic-map draws. Seven characters
describe vein/stigma features, one nygmata and one a muscle. Thus anatomical
overlap with HAIKA forewing outline and venation is substantial. The eight
nonmuscle characters form an additional comparison; the file label "sclerites"
is the source ontology's broad classification, not predominantly wing-base
structures. Hind-wing and all-wing partitions are secondary comparisons.

Each retained HAIKA PC is fitted separately with the exact phytools 2.5.2
multirateBM penalized likelihood. Three bounded optimization starts are used.
The cached analytic-gradient implementation was verified against phytools and
finite differences, and each final fit passes independent likelihood checks.
Node variance rates are converted to branch averages using their logarithmic
mean. TOTAL_VARIANCE is the sum of raw fitted PC diffusivities, an exploratory
trace of retained-feature variance; it is not a jointly fitted multivariate
posterior. The raw units are PC score squared/Myr versus discrete changes/Myr.

Both representations are normalized by their lineage-duration-weighted mean
rate before comparison. Marginal rates average the branch rates using their
duration of overlap with each time bin divided by total lineage exposure.
Both methods use identical branches, bins and weights. The original-segment
ontophylo marginal is retained as a pruning sensitivity. Twenty-eight bins of
at least 5 Myr enter temporal correlations; the tiny residual bin above 280 Ma
is exported but excluded. These are branch-averaged trajectories; plateaus and
bin maxima do not identify the dates of individual evolutionary events.

Primary numerical results
-------------------------
Spearman correlations of mean relative ontophylo rates against HAIKA retained
TOTAL_VARIANCE, seed 3407, smoothing penalty 1:

""" + "\n".join(table) + """

The 77-tip interfamily comparison contains 50 branches whose descendant clade
contains at least two families. Excluding only the two root branches gives
rho=0.084 across 150 branches. The primary forewing ontophylo curve peaks at
230-240 Ma, at 6.50 times its treewide mean; HAIKA is 0.65 times its own mean
in that bin. Its highest 77-tip bin is 240-250 Ma at 1.35 times its mean.
Removing the muscle character yields branch rho=0.157 and time rho=-0.455.

Across the five ANN seeds at penalty 1, full-tree branch correlations range
from 0.052 to 0.215, temporal correlations from -0.463 to +0.196, and interfamily
correlations from 0.083 to 0.482. Therefore not every configuration disagrees
or has negative temporal correspondence. The main finding is weak primary
concordance and sensitivity of stronger claims to representation and modeling.

Conditional null simulations
----------------------------
All 199 simulations completed. A homogeneous multivariate BM model is fitted
on the 48-family representative tree, retaining cross-PC evolutionary covariance.
Simulated family profiles are then copied to the same 77 positions and both
trees are refitted. This controls for shared tree geometry and proxy duplication
under a specific plug-in null, conditional on the fixed PCA and fitted covariance.
It is not a general test of independence between heterogeneous rate histories.

The primary all-branch and interfamily upper-tail probabilities for positive
concordance are 0.895 and 0.955, respectively. Thus those observed correlations
are not unusually positive under this null. The 77-tip temporal association has
an exploratory two-tail probability of 0.020, but the 48-family temporal result
has 0.640, and ANN seeds change the direction of temporal association. This is
not robust evidence of biological anticorrelation. All reported probabilities
are unadjusted exploratory values. Upper-tail values have Monte Carlo resolution
0.005; doubled two-tail values change in steps of 0.010. They apply only to the
matching penalty-1 configuration.

Numerical bounds and smoothing
------------------------------
All 51 planned fits and 18 additional diagnostic fits passed convergence and
likelihood checks. Nevertheless, the selected primary 77-tip fits place 184 of
459 node-axis rates at the lower numerical bound, compared with 11 of 285 for
the 48-tip tree. Positive objective derivatives at these boundaries show that
the objective would decrease locally toward smaller rates if the bound were relaxed.

Moving the lower log-rate offset from 10 to 5 or 15 leaves broad conclusions
similar: total branch-rate ranks correlate 0.971/0.975 with the primary 77-tip
rates, and their forewing correlations are 0.092/0.060. The 48-tip rank agreement
is at least 0.996. However, the wider offset 15 still has 174 boundary parameters
in the 77-tip analysis and one in the 48-tip analysis. Those small individual
rates are not uniquely resolved by these data.

""" + smoothing_text.rstrip() + """

Uncertainty and interpretation
------------------------------
Ontophylo ribbons and map-correlation quantiles describe variation among maps
conditional on fitted discrete-model parameters. They are not joint confidence
or credible intervals for the comparison. The correlation of mean rates is a
different nonlinear summary from quantiles of map-by-map correlations, so the
plotted point can lie outside those whiskers. ANN-seed and smoothing alternatives
are sensitivity scenarios, not HAIKA confidence intervals. Individual BM rates,
especially at poorly constrained basal and duplicated-family branches, should
not be interpreted as uniquely estimated biological values.

The comparison has been completed as far as the available, auditable HAIKA
assignments support: a full 77-tip family-proxy analysis with explicit sampling,
representation and numerical sensitivities. A direct species-level validation
would require a verified image-to-species crosswalk or images for the study
exemplars, followed by extraction in a common trained embedding space.

Files and reproducibility
-------------------------
data/haika-family-proxy-crosswalk.csv: every tip and source-family assignment.
data/haika-family-proxy-source-rows.csv: every selected source image.
data/haika-branch-proxy-flags.csv: branch-level duplication flags.
results/primary-summary.csv: compact numerical comparison and null probabilities.
results/branch-comparison.csv: all branch rates, partitions and map summaries.
results/marginal-comparison.csv: all marginal trajectories and exposure weights.
results/association-summary.csv: complete descriptive association sensitivities.
results/null-association-summary.csv: conditional null comparison.
results/bound-sensitivity.csv: numerical-bound and strong-smoothing checks.
figures/: four editable PDFs, PNG previews, captions and plotting source record.
scripts/: explicit extraction, preparation, fitting, collection and plotting code.
provenance/: source checksums, package versions, checks and frozen primary model.

Substantial fitting and simulations ran on Roihu, project_2004996, in
/scratch/project_2004996/haika-wings. Job 3285454 completed the main grid and null
simulations using eight CPUs; job 3285525 completed numerical sensitivities.
CSC explicitly permits batch jobs in the interactive partition used for these
small allocations. No existing analyses were modified or cancelled.

Sources: doi:10.1126/sciadv.aef0605; https://zenodo.org/records/22715169;
https://github.com/diegosasso/Hymenoptera_paper.
"""
Path("analysis-report.txt").write_text(report)
print("Wrote analysis-report.txt and results/primary-summary.csv")
