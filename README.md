# Oral–Gut Microbial Transmission Analysis

This repository contains the analysis code and processed data for the study of
strain-resolved oral–gut microbial transmission in gastrointestinal cancer
patients and its association with immunotherapy response.

The code reproduces the covariate-adjusted analyses and the figures reported in
the manuscript.

This repository is an **analysis pipeline**, not a standalone software package.
Checklist items intended for compiled software (binary, installer, demo
dataset) are therefore marked “not applicable”; the source code, full
processed data, environment specification and run instructions are provided.

## Repository structure

```
.
├── data/
│   └── clinical_features/          # harmonised clinical metadata used by the adjusted analyses
├── data processing/                # strain-level distance merging
├── fig1/                           # Figure 1 analyses (adjusted PCoA, alpha diversity, MaAsLin2)
├── fig2/                           # Figure 2 analyses (shared SGBs and cumulative abundance)
├── fig3/                           # Figure 3 analyses (transmitted strains and shared SGBs)
├── fig4/                           # Figure 4 analyses (PCoA, pooled LMM, PFS Cox models)
├── fig5/                           # Figure 5 analyses (functional and network analyses)
├── Extended Data Fig/              # Extended Data Figure analyses
├── environment.json                # Code Ocean environment (R/Python package versions)
├── requirements.txt                # Python packages used by the notebooks
├── run                             # master script (Figure 1–4 adjusted analyses)
├── LICENSE                         # MIT License
└── README.md
```

## System requirements

- **Operating system**: Linux (Ubuntu 22.04, Code Ocean image) or Windows 11.
- **R**: tested with R 4.3.2 (Windows 11) and R 4.4.2 (Code Ocean, Ubuntu 22.04).
- **Python**: 3.12.8 (only needed for the prediction notebooks).
- **R packages**: the complete list with versions is given in
  `environment.json`. Main packages include `Maaslin2` (1.20.0), `vegan`
  (2.7-2), `ggplot2` (4.0.0), `survival` (3.8-3), `lme4`, `coxphf`, `dplyr`,
  `tidyr`, `readxl`, `ggpubr`, `patchwork`, `microeco`, `randomForest` and
  `pROC`.
- **Python packages**: `numpy`, `pandas`, `scipy`, `scikit-learn`, `xgboost`
  (see `requirements.txt`).
- **Non-standard hardware**: none. All analyses run on a standard CPU; no GPU
  is required.

## Installation

### Option 1: Code Ocean environment (recommended)

The file `environment.json` defines the complete Code Ocean environment,
including the R/Python versions and package versions. Importing this
environment installs all dependencies automatically.

### Option 2: Local R installation

Install R (>= 4.3) and the required packages, for example:

```r
install.packages(c(
  "dplyr", "tidyr", "readxl", "ggplot2", "ggpubr", "patchwork",
  "vegan", "survival", "survminer", "lme4", "coxphf", "microeco",
  "randomForest", "pROC", "gghalves", "rstatix", "scales"
))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Maaslin2")
```

Typical installation time on a normal desktop computer is **15–30 minutes**,
depending on the internet connection and whether packages need compilation.

## Data

- `data/clinical_features/features_filled_with_metadata.csv` and
  `patient_level_metadata.csv` contain the harmonised clinical covariates
  (age, sex, BMI category, smoking, MMR/MSI status, tumour stage, treatment
  class, histology and PD-L1) used in the adjusted analyses.
- The processed microbiome data (distance matrices, SGB abundance tables and
  RData objects) are stored in the corresponding `fig1`–`fig5` folders.
- `data processing/distance_strainphlan_panphlan.RData` is the merged
  strain-level distance object used to derive nGD.

## Usage

Each analysis script expects to be run from its figure directory. For example:

```bash
# Figure 1
cd fig1
Rscript "adonis2_adjusted/run_adonis2_adjusted.R"
Rscript "adonis2_adjusted/fig1c-d PCoA adjusted.R"
Rscript "fig1bfg_adjusted/Fig1bfg_adjusted.R"
Rscript "maaslin2_adjusted/run_maaslin2_adjusted.R"
Rscript "maaslin2_adjusted/plot_maaslin2_adjusted.R"
Rscript "maaslin2_onevsrest_Diagnosis/run_maaslin2_onevsrest.R"
Rscript "maaslin2_onevsrest_Diagnosis/plot_maaslin2_onevsrest_align_bl.R"
```

```bash
# Figure 2
cd fig2
Rscript "fig2de_revised/run_fig2de_revised.R"

# Figure 3
cd fig3
Rscript "fig3cdef_revised/run_fig3cdef_revised.R"

# Figure 4
cd fig4
Rscript "fig4_adjusted/run_fig4a_adjusted.R"
Rscript "fig4_adjusted/run_fig4c_pooled.R"
Rscript "run_fig4d_revised.R"
Rscript "fig4_adjusted/run_fig4e_pooled.R"
Rscript "replot_fig4c_pooled.R"
```

The scripts write their outputs (PDF, PNG and CSV files) into the
corresponding `*_adjusted` or `*_revised` folder.

For Figure 5 and the Extended Data Figures, run the scripts from their
respective folders (for example, `cd "fig5/butyrate"` and
`Rscript butyrate_plot.R`).

To use your own data, replace the input files in the relevant figure folder
and update the input file names at the top of the corresponding script. All
paths are relative to the figure directory.

## Demo (not applicable)

This is an analysis pipeline rather than a standalone software tool, so no
separate demo dataset or installer is provided. The repository includes the
full real datasets and scripts; see **Usage** for the commands to run the
analyses on the study data.

## Reproduction

After installing the dependencies, the adjusted analyses for Figures 1–4 can
be run from the repository root with:

```bash
bash run
```

The master script changes into each figure directory and executes the
corresponding R scripts. The full run may take several hours because it
includes PERMANOVA, MaAsLin2, mixed models, bootstrap resampling and Cox
models.

A complete description of the code functionality and the statistical models
is provided in the **Methods** section of the manuscript.

## Data availability

All processed data required to reproduce the analyses are included in this
repository. The repository is publicly available at:

https://github.com/xu-research-lab/oral-gut-project

## License

This project is released under the MIT License. See `LICENSE` for details.

## Contact

For questions about the code or data, please open an issue on GitHub or
contact the corresponding author.
