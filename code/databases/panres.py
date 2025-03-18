import pandas as pd
from owlready2 import *

import sys
sys.path.append('..')
from functions import get_instance, clean_gene_name, get_or_create_subclass

database2name = {
    'amrfinderplus': 'AMRFinderPlus',
    'argannot': 'ARGANNOT',
    'card_amr' :'CARD', 
    'csabapal': 'CsabaPal',
    'functional_amr': 'ResFinderFG',
    'megares': 'MegaRes',
    'metalres': 'MetalRes',
    'resfinder': 'ResFinder'
}

def add_panres_genes(file, onto, logger):
    '''
    Loop through panres genes and add them
    '''

    # Load file
    panres_metadata = pd.read_csv(file, sep=',' if file.endswith('.csv') else '\t', skiprows=1)
    panres_metadata['userGeneName'] = panres_metadata['userGeneName'].str.replace("_v1.0.0", "", regex=False)
    panres_metadata['database'] = panres_metadata['database'].str.replace("_genes", "")

    genes = panres_metadata.loc[
        (panres_metadata['database'].isin(['CsabaPal'])), #& ~(panres_metadata['fa_header'].str.contains('Drugs')), 
        'userGeneName'
    ].unique().tolist()

    # genes = panres_metadata['userGeneName'].unique().tolist()
    # genes = genes + ['pan_1875', 'pan_451', 'pan_2299', 'pan_5', 'pan_1204', 'pan_1']


    for gene in set(genes):
        new_gene = onto.PanGene(gene)

        # Get annotation data
        m = panres_metadata.loc[panres_metadata['userGeneName'] == gene, ]

        # Gene length
        if m['gene_len'].nunique() == 1:
            new_gene.has_length.append(int(m['gene_len'].values[0]))
        else:
            logger.warning(f"PanRes: {gene} has multiple gene lengths associated.")
        
        # Loop through the annotation data for this gene
        for _, row in m.iterrows():
            database_shortname = row['database']
            database_name = database2name[database_shortname.lower()]
            database_instance = get_instance(onto, database_name)
            new_gene.is_from_database.append(database_instance)

            fasta_header = row['fa_header'].replace("~~~", "|").replace("'", "")
            gene_name = clean_gene_name(fasta_header, database_shortname.lower()).strip()
            original_gene_instance = onto.OriginalGene(gene_name)
                        
            original_gene_instance.is_from_database.append(database_instance)
            original_gene_instance.original_fasta_header.append(fasta_header.strip() + '|' + database_name)

            new_gene.same_as.append(original_gene_instance)
            # original_gene_instance.subclass_of.append(new_gene)
    
    # logger.info(f"Adding {len(genes)} PanRes genes to the ontology.")
    logger.success(f"Added PanRes genes (n={len(genes)}) to the ontology.")