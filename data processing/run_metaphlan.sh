#!/bin/bash
#SBATCH --job-name=strainphlan_analysis
#SBATCH --output=logs/strainphlan_%A_%a.out
#SBATCH --error=logs/strainphlan_%A_%a.err
#SBATCH --array=1-$(wc -l < sample.list)
#SBATCH --cpus-per-task=7
#SBATCH --mem=64G
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
OUTPUT_DIR="/beegfs/huangxiaochang/Oral_gut_cancer/mpa4/profile"

# 1. Run MetaPhlAn to get microbiome profile
echo "Step 1: Running MetaPhlAn for $SAMPLE"
metaphlan \
    "${RAW_DATA_DIR}/${SAMPLE}_R1.fastq,${RAW_DATA_DIR}/${SAMPLE}_R2.fastq" \
    --bowtie2out "${OUTPUT_DIR}/${SAMPLE}.bowtie2.txt" \
    --samout "${OUTPUT_DIR}/${SAMPLE}.sam.bz2" \
    --nproc 28 \
    --input_type fastq \
    > "${OUTPUT_DIR}/${SAMPLE}_profiled_metagenome.txt"
