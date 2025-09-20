from micom.workflows import build, grow
from micom.qiime_formats import load_qiime_medium
from micom.viz import plot_association
from micom.measures import production_rates
import pandas as pd
import pickle
import numpy as np

def main():
    species_abundance = pd.read_csv('species_abundance_micom.csv')
    species_abundance['species'] = 's__' + species_abundance['species']

    manifest = build(species_abundance, "agora103_species.qza", "models_agora_all_level", 
                     cutoff=1e-2, threads=20)

    medium = load_qiime_medium("western_diet_gut.qza")

    filtered_df = manifest[manifest['sample_id'].str.startswith('Feces')]

    growth_results = grow(filtered_df, "models_agora_all_level", medium, tradeoff=0.5, threads=20)

    with open("growth_new.pickle", "wb") as f:
        pickle.dump(growth_results, f)

    phenotype = pd.read_table('phenotype.txt')
    phenotype = phenotype.rename(columns={'id': 'sample_id'})
    phenotype.set_index('sample_id', inplace=True)
    phenotype.iloc[:, 0] = pd.to_numeric(phenotype.iloc[:, 0])

    p1 = plot_association(
        growth_results,
        phenotype= phenotype['transmitted abundance in Gut'],
        variable_type="continuous",
        filename="association.html",
        fdr_threshold=0.25,
    )

if __name__ == '__main__':
    import multiprocessing
    multiprocessing.freeze_support()
    main()
