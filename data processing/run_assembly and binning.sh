#!/bin/bash
#SBATCH --job-name=genome_binning
#SBATCH --output=logs/binning_%A_%a.out
#SBATCH --error=logs/binning_%A_%a.err
#SBATCH --array=1-$(wc -l < sample.list)
#SBATCH --cpus-per-task=28
#SBATCH --mem=64G
#SBATCH --time=48:00:00

# Get current sample from sample list
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample.list)
echo "Processing sample: $SAMPLE (Task ID: $SLURM_ARRAY_TASK_ID)"

# Set directory paths
RAW_DATA_DIR="/home/raw-data/2021-11_oral-gut-cancer/Sequencing_Rawdata/cleandata"
HOST_REF_DIR="/beegfs/huangxiaochang/database/H.sapiens/GRCh38_noalt_as"
OUTPUT_DIR="analysis_results"
HOST_REMOVED_DIR="${OUTPUT_DIR}/host_removed"
ASSEMBLY_DIR="${OUTPUT_DIR}/assembly"
BINNING_DIR="${OUTPUT_DIR}/binning"
QUANT_DIR="${OUTPUT_DIR}/quantification"

# Create output directories
mkdir -p "${HOST_REMOVED_DIR}" "${ASSEMBLY_DIR}" "${BINNING_DIR}" "${QUANT_DIR}" "logs"

###############################################################################
# Step 1: Remove host DNA from sequencing data
###############################################################################
echo "Step 1: Removing host DNA from $SAMPLE"

bowtie2 --threads 28 \
    -x "${HOST_REF_DIR}/GRCh38_noalt_as" \
    -1 "${RAW_DATA_DIR}/${SAMPLE}_1.clean.fq.gz" \
    -2 "${RAW_DATA_DIR}/${SAMPLE}_2.clean.fq.gz" \
    -S "${HOST_REMOVED_DIR}/${sample}_mapped_and_unmapped.sam"

# Convert SAM to BAM
samtools view -bS "${HOST_REMOVED_DIR}/${sample}_mapped_and_unmapped.sam" \
    > "${HOST_REMOVED_DIR}/${sample}_mapped_and_unmapped.bam"

# Extract unmapped reads (both reads unmapped)
samtools view -b -f 12 -F 256 "${HOST_REMOVED_DIR}/${sample}_mapped_and_unmapped.bam" \
    > "${HOST_REMOVED_DIR}/${sample}_bothReadsUnmapped.bam"

# Sort by read name for paired-end output
samtools sort -n -@ 4 "${HOST_REMOVED_DIR}/${sample}_bothReadsUnmapped.bam" \
    > "${HOST_REMOVED_DIR}/${sample}_bothReadsUnmapped_sorted.bam"

# Convert to FASTQ format
samtools fastq -@ 4 \
    "${HOST_REMOVED_DIR}/${sample}_bothReadsUnmapped_sorted.bam" \
    -1 "${HOST_REMOVED_DIR}/${sample}_host_removed_R1.fastq.gz" \
    -2 "${HOST_REMOVED_DIR}/${sample}_host_removed_R2.fastq.gz"

# Clean up intermediate files
rm -f "${HOST_REMOVED_DIR}/${sample}_mapped_and_unmapped.sam" \
      "${HOST_REMOVED_DIR}/${sample}_mapped_and_unmapped.bam" \
      "${HOST_REMOVED_DIR}/${sample}_bothReadsUnmapped.bam" \
      "${HOST_REMOVED_DIR}/${sample}_bothReadsUnmapped_sorted.bam"

###############################################################################
# Step 2: Metagenomic assembly
###############################################################################
echo "Step 2: Assembling metagenome for $SAMPLE"

metawrap assembly \
    -1 "${HOST_REMOVED_DIR}/${sample}_host_removed_R1.fastq.gz" \
    -2 "${HOST_REMOVED_DIR}/${sample}_host_removed_R2.fastq.gz" \
    -m 130 \
    -t 28 \
    --metaspades \
    -o "${ASSEMBLY_DIR}/${sample}_assembly"

###############################################################################
# Step 3: Metagenomic binning
###############################################################################
echo "Step 3: Binning metagenomic contigs for $SAMPLE"

# Note: This step requires pre-assembled contigs and quality-checked reads
# Adjust paths according to your specific data structure
metaWRAP binning \
    --metabat2 \
    --maxbin2 \
    --concoct \
    --universal \
    -t 28 \
    -l 1000 \
    -m 200 \
    -a "${ASSEMBLY_DIR}/${sample}_assembly/contigs.fasta" \
    -o "${BINNING_DIR}/${sample}_initial_bins" \
    "${HOST_REMOVED_DIR}/${sample}_host_removed_R1.fastq.gz" \
    "${HOST_REMOVED_DIR}/${sample}_host_removed_R2.fastq.gz"

###############################################################################
# Step 4: Bin refinement
###############################################################################
echo "Step 4: Refining bins for $SAMPLE"

metawrap bin_refinement \
    -o "${BINNING_DIR}/${sample}_refined_bins" \
    -c 70 \
    -x 10 \
    --skip-refinement \
    -t 28 \
    -A "${BINNING_DIR}/${sample}_initial_bins/metabat2_bins" \
    -B "${BINNING_DIR}/${sample}_initial_bins/concoct_bins" \
    -C "${BINNING_DIR}/${sample}_initial_bins/maxbin2_bins"

###############################################################################
# Step 5: Bin quantification
###############################################################################
echo "Step 5: Quantifying bin abundance for $SAMPLE"

metaWRAP quant_bins \
    -t 28 \
    -o "${QUANT_DIR}/${sample}_quant" \
    -b "${BINNING_DIR}/${sample}_refined_bins/metawrap_70_10_bins" \
    "${HOST_REMOVED_DIR}/${sample}_host_removed_R1.fastq.gz" \
    "${HOST_REMOVED_DIR}/${sample}_host_removed_R2.fastq.gz"

echo "Genome binning pipeline completed for sample: $SAMPLE"