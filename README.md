# Oral-Gut Microbial Transmission Analysis Pipeline

A comprehensive pipeline for analyzing microbial strain transmission between oral and gut microbiomes using metagenomic data.

## Quick Start

### 1. Environment Setup
```bash
conda create -n strain-transmission -c bioconda metaphlan strainphlan panphlan phylophlan bowtie2 samtools metawrap drep r
conda activate strain-transmission
```

### 2. Sample Processing
**Input**: `sample.list` (list of sample names)

**Run processing:**
```bash
# SNV-based strain analysis
sbatch run_strainphlan_analysis.sh

# Assembly and binning
sbatch run_genome_binning.sh
```

### 3. Genome Processing
```bash
# Dereplicate MAGs
sbatch run_dereplication.sh

# Map bins to SGBs
find analysis_results/dereplication/ -name "*.fa" > bin_list.txt
sbatch process_all_bins.sh
```

### 4. Gene Content Analysis
```bash
# Expand database and profile
sbatch integrate_bins_to_panphlan.sh
sbatch run_panphlan_profiling.sh
```

### 5. Transmission Detection
**Input**: `distance.RData`, `kept_dis.RData`
```bash
sbatch run_ngd_analysis.sh
```