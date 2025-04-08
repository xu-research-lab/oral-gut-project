#!/bin/bash
#SBATCH --job-name=dereplication
#SBATCH --output=logs/dereplication_%A.out
#SBATCH --error=logs/dereplication_%A.err
#SBATCH --cpus-per-task=28
#SBATCH --mem=256G


# Activate conda environment for dRep
source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate mpa

# Set directory paths
BINNING_DIR="analysis_results/binning"
DEREP_DIR="analysis_results/dereplication"

# Create output directory
mkdir -p "${DEREP_DIR}" "logs"

# Generate list of all bins from all samples
find "${BINNING_DIR}" -name "*.fa" -o -name "*.fasta" -o -name "*.fna" > "${DEREP_DIR}/bin_list.txt"

###############################################################################
# Step 6: Dereplicate genome bins
###############################################################################
echo "Step 6: Dereplicating genome bins"

dRep dereplicate \
    -g "${DEREP_DIR}/bin_list.txt" \
    --S_algorithm fastANI \
    --multiround_primary_clustering \
    --clusterAlg greedy \
    -ms 10000 \
    -pa 0.9 \
    -sa 0.98 \
    -nc 0.50 \
    -cm larger \
    -p 28 \
    -o "${DEREP_DIR}/dereplicated_bins"

echo "Dereplication completed. Results in: ${DEREP_DIR}/dereplicated_bins"