# Figure 6 — Mouse validation analyses

This folder contains the analysis code and processed data for the mouse
validation experiments: tumour growth, flow cytometry, plasma SCFA
metabolomics and SGB-level LEfSe.

## Folder structure

```
fig6/
├── growth_flow/              # tumour growth curves and flow-cytometry panels
│   ├── growth_flow_plots.R
│   ├── mouse_growth_flow_data.xlsx
│   └── growth_curve.pdf/png, flow_3panel.pdf/png, ...
├── scfa_metabolomics/        # plasma SCFA metabolomics
│   ├── run_scfa_metabolomics.R
│   ├── group_config.R
│   ├── meta_intensity_fill_plasma.xls
│   └── scfa_meta_focus_4_metabolites_bar.pdf/png, ...
└── lefse/                    # SGB-level LEfSe: PD-1 vs PD-1 + VP
    ├── run_lefse.R
    ├── plot_lefse_mirrored.R
    ├── group_config.R        # contains the sample-to-group map
    ├── annotate_sgb.R
    ├── merged_abundance_table.tsv
    ├── lefse_sgb_only_res_diff.tsv   # LEfSe result used by the plotting script
    └── lefse_SGB_LDA_bar_mirrored_pd1_vs_vp.pdf
```

## How to run

Run each script from its analysis folder:

```bash
# Tumour growth and flow cytometry
cd fig6/growth_flow
Rscript growth_flow_plots.R

# Plasma SCFA metabolomics
cd fig6/scfa_metabolomics
Rscript run_scfa_metabolomics.R

# LEfSe and mirrored LDA plot
cd fig6/lefse
Rscript run_lefse.R
Rscript plot_lefse_mirrored.R
```

All input paths are relative to the corresponding analysis folder.

## Input files

- `growth_flow/mouse_growth_flow_data.xlsx`: mouse tumour-growth and
  flow-cytometry workbook (sheets 1–4).
- `scfa_metabolomics/meta_intensity_fill_plasma.xls`: plasma SCFA
  intensities (tab-separated).
- `scfa_metabolomics/group_config.R`: group order and colours for the SCFA
  plots.
- `lefse/merged_abundance_table.tsv`: merged MetaPhlAn abundance table used
  for LEfSe.
- `lefse/group_config.R`: sample-to-group map used by LEfSe.
- `lefse/annotate_sgb.R`: SGB label annotation helper.

## Main outputs

- `growth_flow/growth_curve.pdf` and `growth_flow/flow_3panel.pdf`
- `scfa_metabolomics/scfa_meta_focus_4_metabolites_bar.pdf`
- `lefse/lefse_SGB_LDA_bar_mirrored_pd1_vs_vp.pdf`

Intermediate tables and other diagnostic plots are regenerated when the
scripts are run and are not included here.

## Notes

- All filenames are in English and all paths are relative.
- Required R packages include `microeco`, `data.table`, `readxl`, `dplyr`,
  `tidyr`, `patchwork`, `lme4`, `emmeans`, `ggplot2`, `ggpubr`, `pheatmap`,
  `vegan` and `rstatix`.
