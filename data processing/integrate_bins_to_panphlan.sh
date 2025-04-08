#!/bin/bash
#SBATCH --job-name=panphlan_integration
#SBATCH --output=logs/panphlan_integration_%A.out
#SBATCH --error=logs/panphlan_integration_%A.err
#SBATCH --cpus-per-task=28
#SBATCH --mem=256G

# Activate conda environment for PanPhlAn
source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate mpa4

# Set directory paths
SGB_ASSIGNMENTS="analysis_results/phylophlan/MAG_mapping_to_SGB_list.txt"
DEREP_DIR="/beegfs/huangxiaochang/Oral_gut_cancer/Derep/dereplicate/dereplicated_genomes_0.98"
PANPHLAN_DB_DIR="/beegfs/huangxiaochang/database/pangenome"
SAMPLE_DIR="/home/huangxc/oral_gut/host_removed"
UNIREF_DB="/beegfs/huangxiaochang/database/UniRef"  

# Create temporary directories
mkdir -p "logs" "temp_species_bins"

###############################################################################
# Step 1: Parse Species assignments and group bins by Species
###############################################################################
echo "Step 1: Grouping bins by species from $SGB_ASSIGNMENTS"

# Read the mapping file and group bins by species
awk -F'\t' 'NR>1 && $2 != "" && $3 != "" {
    species = $2;
    bin_name = $1;
    # Remove path if present and get just the bin filename
    gsub(/.*\//, "", bin_name);
    print species, bin_name;
}' "${SGB_ASSIGNMENTS}" | while read -r species bin_name; do
    # Create species directory in temp location
    species_dir="temp_species_bins/${species}"
    mkdir -p "${species_dir}"
    
    # Find the bin file in dereplicated genomes
    bin_file=$(find "${DEREP_DIR}" -name "${bin_name}*" -o -name "${bin_name}*.fa" | head -1)
    
    if [ -f "$bin_file" ]; then
        # Copy bin to species directory
        cp "$bin_file" "${species_dir}/${bin_name}.fna"
        echo "Copied $bin_name to $species directory"
    else
        echo "Warning: Bin file not found for $bin_name"
    fi
done

###############################################################################
# Step 2: Expand PanPhlAn database for each species with new bins
###############################################################################
echo "Step 2: Expanding PanPhlAn databases with new bins"

# Use panphlan_exporter.py to add bins directly to the original database directory
python panphlan_exporter.py \
    --input "$species_dir" \
    --output "$original_db_dir" \
    --db_path "$UNIREF_DB" \
    --clade_name "$species" \
    --nprocs 14 \
    --tmp "${original_db_dir}/tmp"

###############################################################################
# Step 3: Run PanPhlAn gene content profiling with expanded database
###############################################################################
echo "Step 3: Running PanPhlAn gene content profiling"

# Get list of samples to process
SAMPLE_LIST="sample.list"

# Process each sample with each expanded species database
while read -r sample; do
    echo "Processing sample: $sample"
    
    sample_output_dir="analysis_results/panphlan_profiles/${sample}"
    mkdir -p "$sample_output_dir"
    
    # Process each species database
    find "${PANPHLAN_DB_DIR}" -mindepth 1 -type d | while read -r species_db_dir; do
        species=$(basename "$species_db_dir")
        
        # Check if database files exist
        if [ -f "${species_db_dir}/${species}_pangenome.tsv" ] && [ -f "${species_db_dir}/${species}.bt2" ]; then
            echo "Mapping sample $sample to $species database"
            
            # Run PanPhlAn mapping
            panphlan_map.py \
                -i "${SAMPLE_DIR}/${sample}_host_removed_R1.fastq.gz,${SAMPLE_DIR}/${sample}_host_removed_R2.fastq.gz" \
                --nproc 14 \
                -m 8 \
                --indexes "${species_db_dir}/${species}" \
                -p "${species_db_dir}/${species}_pangenome.tsv" \
                -o "${sample_output_dir}/${sample}_${species}_profile.tsv"
        fi
    done
    
done < "$SAMPLE_LIST"

###############################################################################
# Step 4: Merge results from each sample
###############################################################################
echo "Step 4: Merging results from all samples"

# Create output directory for merged matrices
mkdir -p "analysis_results/merged_panphlan_matrices"

# Use panphlan_profiling.py to merge results
panphlan_profiling.py \
    -p "${species_db_dir}/${species}_pangenome.tsv" \
    -i "analysis_results/panphlan_profiles" \
     --o_matrix "analysis_results/merged_panphlan_matrices/${species}_matrix"
        

