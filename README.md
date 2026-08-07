
Manuscript Title: Phenome-wide adaptive radiation drives lineage diversification in hyperdiverse Hymenoptera
Authors: Diego S. Porto, Lars Vilhelmsen, István Mikó, Sergei Tarasov
Last Updated: 2026-08-06

# Overview
This repository contains all data, scripts, and outputs associated with the manuscript.
The root folder contains all R scripts for reproducing the analyses of the manuscript. 
This main folder is organized into separated subfolders for raw data, processed data, output data, R functions, raw figures, etc. as detailed below.

# Folder structure
```
/
├── alluvial-plot/
├── Archive/
├── data/
├── data_mds/
├── data_nhpp/
├── data_out/
├── data_pyrate/
├── figures/
├── make-BAMM/
├── ml-morpho/
├── R/
├── stmaps/
├── stmaps_amalg/
├── corr-mol-vs-morph.Rmd
├── corr-rates-and-metrics.Rmd
├── make-adult-marginal-rate.Rmd
├── make-branch-specific-rates.Rmd
├── make-Neff-minus.Rmd
├── make-Neff-plus.Rmd
├── make-Neff-plus-sensitivity.Rmd
├── make-Pc.Rmd
├── make-pyrate-data.Rmd
├── miscellaneous-analyses.Rmd
├── ontophylo_base_analyses.Rmd
├── ontophylo_core_analyses.Rmd
├── organize_data.Rmd
├── pairwise-corr.Rmd
├── plot-enrich-over-16.Rmd
├── plot-pyrate-Neff.Rmd
├── plot-tip-changes.Rmd
├── plot-tree.Rmd
├── toy-example.Rmd
└── README
```

# Brief explanation of folders and main files

## Main subfolders

### alluvial-plot/

Contains R script to produce a subplot of figure 1 and its output raw figure file.

- anatomy_tree.pdf - raw pdf file of the hierarchical anatomy plot of figure 1.

- anatomy_tree.R - script to produce the hierarchical anatomy plot of figure 1.


### Archive/

Contains subfolders with files from earlier versions of the analyses of the manuscript. Not used in the final version of the manuscript. Please, ignore.


### data/

Contains all input data, raw or modified from previous publications, used in the analyses.

- adult_matrix.RDS - R object with the phylogenetic matrix of adult characters modified from [Sharkey et al. (2012)](https://doi.org/10.1111/j.1096-0031.2011.00366.x).

- bamm.RDS - R object with the BAMM diversification data modified from [Blaimer et al. (2023)](https://doi.org/10.1038/s41467-023-36868-4).

- char_annot.RDS - R object with the annotations to ontology terms associated with the adult characters.

- char_num.RDA - R objects with the counts of characters and states for each phenomic unit analyzed.

- HAO.obo - OBO file with a snapshot from the Hymenoptera Anatomy Ontology ([HAO](https://github.com/hymao/hao/blob/master/hao.obo)).

- hym_tree.RDS - R object with the dated phylogenetic tree modified rom topology C1 of [Blaimer et al. (2023)](https://doi.org/10.1038/s41467-023-36868-4).

- HYM-NUC-70%-SWSC.tre - file with a phylogenetic tree obtained from [Blaimer et al. (2023)](https://doi.org/10.1038/s41467-023-36868-4).

- larval_matrix.RDS - R object with the phylogenetic matrix of larval characters modified from [Ronquist et al. (2012)](https://doi.org/10.1093/sysbio/sys058).

- taxon_data.RDS - R object with the taxonomic information matches between the [Sharkey et al. (2012)](https://doi.org/10.1111/j.1096-0031.2011.00366.x) and [Blaimer et al. (2023)](https://doi.org/10.1038/s41467-023-36868-4) datasets.


### data_mds/

Contains output data from the morphospace reconstructions.

- phenome_full_MDS.RDS - R object with the morphospace reconstruction of the adult phenome.


### data_nhpp/

Contains output data from the ontophylo core analyses.

- adult_anat_ent_*.RDS - R objects with the continuously-varying branch-specific evolutionary rates inferred for all 15 elementary anatomical entities of the adult phenome.

- adult_body_reg_*.RDS - R objects with the continuously-varying branch-specific evolutionary rates inferred for all 5 major body regions of the adult phenome.

- adult_phenome_phenome_full.RDS - R object with the continuously-varying branch-specific evolutionary rates inferred for the full adult phenome.

- larval.RDS - R object with the continuously-varying branch-specific evolutionary rates inferred for the larval characters.


### data_nhpp/

Contains output data from several analyses, as detailed below.

br_rates_all.RDS - R object with the branch-specific rates for each of the 15 elementary anatomical regions of the adult phenome for all stochastic-mapping samples.

br_rates_all_ind.RDS - R object with the mean branch-specific rates for each of the 15 elementary anatomical regions of the adult phenome.

hym_dataset.RDA - R objects with processed data (phylogenetic tree, adult and larval character datasets, ontology annotations) for downstream analyses.

neff-92-5.RDS - R object with the output data from the Neff_plus sensitivity analysis with upper threshold of 92.5%.

neff-97-5.RDS - R object with the output data from the Neff_plus sensitivity analysis with upper threshold of 97.5%.

neff-per-replicate.RDS - R object with the output data from the Neff_plus calculations.

neff-under-per-replicate.RDS - R object with the output data from the Neff_minus calculations.

paramo_amalg_adult.RDA - R objects with data from ontophylo analyses (character amalgamation) of the adult characters.

paramo_stm_adult.RDA - R objects with raw data from ontophylo analyses (stochastic character mapping) of the adult characters.

paramo_stm_adult_final.RDA - R objects with processed data from base ontophylo analyses (stochastic character mapping) of the adult characters.

paramo_stm_larval.RDA - R objects with raw data from base ontophylo analyses (stochastic character mapping) of the larval characters.

paramo_stm_larval_final.RDA - R objects with processed data from base ontophylo analyses (stochastic character mapping) of the larval characters.

pc-per-replicate.RDS - R object with the output data from the phenomic coordination (Pc) calculations.

pyrate-dive.RDS - R object with the diversification data processed from [Jouault et al. (2025)](https://doi.org/10.1016/j.cub.2025.03.002).

pyrate-lambda.RDS - R object with the speciation/origination data processed from [Jouault et al. (2025)](https://doi.org/10.1016/j.cub.2025.03.002).

pyrate-mu.RDS - R object with the extinction data processed from [Jouault et al. (2025)](https://doi.org/10.1016/j.cub.2025.03.002).

quality_control_adult.RDA - R objects with the outputs from the quality control of stochastic-mapping analyses of the adult characters.

quality_control_larval.RDA - R objects with the outputs from the quality control of stochastic-mapping analyses of the larval characters.

rates-perm.RDS - R object with the output from the permutation tests of rate correlations among the 15 elementary anatomical regions.

rates-thru-time.RDS - R object with marginalized rate-through-time data for the adult characters.


### data_pyrate/BDCS/BDCS-fam-ssing-20/

Contains log files and raw data from the BDCS-fam-ssing-20 PyRate analyses of [Jouault et al. (2025)](https://doi.org/10.1016/j.cub.2025.03.002).


### figures/

Contains raw figures used to prepare the final plates of the manuscript.

- analyses_checks/
	- char_rates_final.png - Raw export of Supplementary Figure 1b.
	
	- char_rates_initial.png - Raw export of Supplementary Figure 1a.

- BAMM/BAMM_rates.pdf - Raw export of Supplementary Figure 8.

- cor/pairwise-corr.pdf - Raw export of Figure 3b.

- evo_rates/
	- adult_anat_ent_\*_edgeplot.* - Raw exports of Figure 3a.

	- adult_body_reg_\*_edgeplot.* - Raw exports of Supplementary Figure 2.
	
	- adult_phenome_phenome_full_edgeplot.* - Raw exports of Figure 2a.

	- larval_edgeplot.* - Raw exports of Figure 2b.

- metrics_thru_time/
	- marginal_rate.pdf - Raw export of Figure 4c.

	- neff_95.pdf - Raw export of Supplementary Figure 3.

	- neff_925.pdf - Raw export of Supplementary Figure 4.

	- neff_975.pdf - Raw export of Supplementary Figure 5.

	- neff-minus.pdf - Raw export of part of Figure 4d.
	
	- neff-pyrate-bin-20.pdf - Raw export of Figure 4b.
	
	- neff-pyrate-bin-50.pdf - Alternative raw export of Figure 4b with 50-million-year bins (not used in the manuscript).
	
	- neff-pyrate-steps.pdf - Raw export of Extended Data Figure 3b.
	
	- Pc-metric.pdf - Raw export of Extended Data Figure 2 and part of Figure 4d.
	
	- pyrate_vs_neff.pdf - Raw export of Extended Data Figure 3a.

- morphospaces/phenome_full_mds.* - Raw exports of Figure 2c.

- tree-enrich/

	- heat_map-legend.pdf - Raw export of part of Figure 4e.
	
	- hym-prop-changes-enriched.pdf - Raw export of part of Figure 4f.
	
	- hym-proportion-enriched-matrix.pdf - Raw export of part of Extended Data Figure 6.
	
	- hym-prop-time-enriched.pdf - Raw export of part of Figure 4f.
	
	- tree_enrich.pdf - Raw export of Figure 4a and 4e.

- tree-enrich-ai/

	- tree-enrich.ai - Illustrator file of post-processed composition of subfigures 4a, b, and e.
	
	- tree-enrich.pdf - Alternative export of post-processed composition of subfigures 4a, b, and e.

- lambda-mu-rate.pdf - Raw export of Extended Data Figure 5.

- morph_vs_age.pdf - Raw export of Supplementary Figure 7.

- morph_vs_mol.pdf - Raw export of Supplementary Figure 6.

- tree.pdf - Raw export of the phylogenetic tree used as a subplot in Extended Data Figures 2 and 5.


### make-BAMM/

Contains R script for processing data and raw data from BAMM analysis of [Blaimer et al. (2023)](https://doi.org/10.1038/s41467-023-36868-4).

- data/ - raw data from BAMM analyses of [Blaimer et al. (2023)](https://doi.org/10.1038/s41467-023-36868-4).

- make-BAMM-data.R - R script for processing raw BAMM data and plotting Supplementary Figure 8.


### ml-morpho/

Contains R script for preparing input data and the output data from the topology-constrained maximum-likelihood analysis of the adult characters modified from [Sharkey et al. (2012)](https://doi.org/10.1111/j.1096-0031.2011.00366.x)

- iqtree/ - output data from the maximum-likelihood analysis.

- make-data.R - R script for preparing input data for maximum-likelihood analysis.


### R/

Contains several R functions used in the analyses.

- ama-branch-specific.R - Script for testing the utils-branch-ama.R functions.

- helpers.R - General helper functions for model selection, amalgamation, stochastic mapping, plotting, etc.

- region_class.R - Defines R objects containing information about anatomical partitions.

- utils-branch-ama.R - Functions for calculating branch-specific rates.

- utils-edge-age.R - Utility function to extract edge ages.

- utils-maddfitz-root.R - Utility function to get root prior from corHMM.

- utils-ST.R - Utility functions to perform various analyses (rate correlations, pc calculation, neff calculation, etc.).


### stmaps/

Contains output data from the base ontophylo analyses - stochastic character mapping.

- larval/CH*.RDS - subfolder with R objects with the output of stochastic character mapping of the larval characters.

- CH*.RDS - R objects with the output of stochastic character mapping of the adult characters.


### stmaps_amalg/

Contains output data from the ontophylo analyses - character amalgamation.

- adult_anat_ent.RDS - R object with the output from character amalgamations of each of the 15 elementary anatomical entities of the adult phenome.

- adult_body_reg.RDS - R object with the output from character amalgamations of each of the 5 major body regions of the adult phenome.

- adult_phenome.RDS - R object with the output from character amalgamation of the adult phenome.

- larval.RDS - R object with the output from the amalgamation of the larval characters.



## Analysis scripts

- corr-mol-vs-morph.Rmd - R script to perform correlation analysis between morphological and molecular branch lengths (Supplementary Note 4) and plotting raw figures (Supplementary Figures 6 and 7).

- corr-rates-and-metrics.Rmd - R script to perform correlation analyses between diversification rates (PyRate, BAMM) and phenomic metrics (rates, Neff, Pc) and plotting raw figures (Extended Data Figure 3a).

- make-adult-marginal-rate.Rmd - R script for calculating branch-specific rates of the adult phenome marginalized across all lineages and plotting raw figures (Figure 4c).

- make-branch-specific-rates.Rmd - R script for calculating branch-specific rates for for each of the 15 elementary anatomical regions of the adult phenome.

- make-Neff-minus.Rmd - R script for calculating the Neff_minus metric for the adult phenome and plotting raw figures (Figure 4d).

- make-Neff-plus.Rmd - R script for calculating the Neff_plus metric for the adult phenome and plotting raw figures (Supplementary Figure 3).

- make-Neff-plus-sensitivity.Rmd - R script to perform sensitivity analyses of Neff metric to alternative thresholds (Supplementary Note 2) and plotting raw figures (Supplementary Figures 4 and 5).

- make-Pc.Rmd - R script for calculating the phenomic coordination (Pc) metric for the adult phenome and plotting raw figures (Figure 4d and Extended Data Figure 2).

- make-pyrate-data.Rmd - R script for preparing PyRate data for downstream analyses and plotting raw figures (Extended Data Figure 5).

- miscellaneous-analyses.Rmd - R script for calculating rate magnitudes and Hamming distances.

- ontophylo_base_analyses.Rmd - R script to perform base ontophylo analyses including model fitting, stochastic mapping, quality controls, and plotting raw figures (Supplementary Figure 1a-b).

- ontophylo_core_analyses.Rmd - R script to perform core ontophylo analyses including character amalgamation, rate estimation, morphospace reconstruction, and plotting raw figures (Figures 2a-c, 3a, Supplementary Figure 2).

- organize_data.Rmd - Starting R script to organize data for downstream analyses in the ontophylo workflow.

- pairwise-corr.Rmd - R script to perform correlation analyses among all elementary anatomical entities of the adult phenome and larval phenome, permutation tests, and plotting raw figures (Figure 3b).

- plot-enrich-over-16.Rmd - R script for plotting raw figures of enriched accelerated branches (Neff_plus) and tip enrichment (Figures 4a, e).

- plot-pyrate-Neff.Rmd - R script for plotting Neff_plus and PyRate curves-through-time (Figure 4b and Extended Data Figure 3b).

- plot-tip-changes.Rmd - R script for plotting tip enrichment data (Figure 4a, e, f and Extended Data Figure 6).

- plot-tree.Rmd - R script for plotting the phylogenetic tree used as a subplot in Extended Data Figures 2 and 5.

- toy-example.Rmd - R script for making Neff data used in the toy example of Extended Data Figure 4.

------------------------------

# Reproducibility instructions

## Required R packages
Install all required packages and dependencies (if needed) from within your R environment:

```r
install.packages("remotes")
remotes::install_github("phenoscape/rphenoscape")
remotes::install_github("uyedaj/rphenoscate")
install.packages(c("parallel", "ape", "phytools", "corHMM", "geiger", "ontophylo", "tidyverse", "ggpubr", "stringdist", "viridis", "deeptime", "ontologyIndex", "abind", "relaimpo", "reshape2", "mgcv", "bnpsd"))
```

## Running the full pipeline
Run the scripts below in the indicated order.

### Base phenomic analyses - ontophylo workflow (in order):
	- 01. organize_data.Rmd
	- 02. ontophylo_base_analyses.Rmd
	- 03. ontophylo_core_analyses.Rmd
	- 04. make-branch-specific-rates.Rmd
	
### Phenomic metrics (in order):
	- 05. make-Neff-plus.Rmd
	- 06. make-Neff-minus.Rmd
	- 07. make-Neff-plus-sensitivity.Rmd
	- 08. make-Pc.Rmd

### Evolutionary rates and diversification correlation analyses (in order):
	- 09. make-adult-marginal-rate.Rmd
	- 10. make-pyrate-data.Rmd
	- 11. pairwise-corr.Rmd
	- 12. corr-rates-and-metrics.Rmd

### Miscellanous analyses (no order):
	- corr-mol-vs-morph.Rmd
	- miscellaneous-analyses.Rmd

### Plots (no order):
	- toy-example.Rmd
	- plot-tree.Rmd
	- plot-tip-changes.Rmd
	- plot-pyrate-Neff.Rmd
	- plot-enrich-over-16.Rmd
