#!/bin/bash
#SBATCH --job-name=strainphlan_analysis
#SBATCH --output=logs/strainphlan_%A_%a.out
#SBATCH --error=logs/strainphlan_%A_%a.err
#SBATCH --array=1-$(wc -l < sample.list)
#SBATCH --cpus-per-task=28
#SBATCH --mem=256G
#SBATCH --time=24:00:00

# Load conda and activate environment
source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate mpa4

# Get current sample from sample list
SAMPLE_LIST="sample.list"
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
echo "Processing sample: $SAMPLE (Task ID: $SLURM_ARRAY_TASK_ID)"

# Set directory paths
RAW_DATA_DIR="/home/huangxc/oral_gut/host_removed"
OUTPUT_DIR="/beegfs/huangxiaochang/Oral_gut_cancer/mpa4/"
CONSENSUS_MARKER_DIR="${OUTPUT_DIR}/consensus_marker"
STRAINPHLAN_DIR="${OUTPUT_DIR}/strainphlan"
SGB_DIR="${OUTPUT_DIR}/SGB"

# 1. Extract marker genes for strain-level analysis
echo "Step 2: Extracting consensus markers for $SAMPLE"
sample2markers.py \
    -i "${OUTPUT_DIR}/${SAMPLE}.sam.bz2" \
    -o "${CONSENSUS_MARKER_DIR}" \
    -n 7 \
    --dominant_frq_threshold 0.7 \
    --min_reads_aligning 5 \
    -b 50

# 2. Extract markers for specific sample
echo "Step 3: Extracting sample-specific markers for $SAMPLE"
extract_markers.py \
    -c "$SAMPLE" \
    -o "$SGB_DIR"

# 3. Run StrainPhlAn for all samples
echo "Step 4: Running StrainPhlAn for all samples of $SAMPLE"
strainphlan \
    -s "${CONSENSUS_MARKER_DIR}/*.pkl" \
    -c "$SAMPLE" \
    -m "${SGB_DIR}/${SAMPLE}.fna" \
    --marker_in_n_samples 10 \
    --sample_with_n_markers 10 \
    -o "${STRAINPHLAN_DIR}/nonPMA" \
    -n 7
