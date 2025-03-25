import pandas as pd
from owlready2 import *
import numpy as np
from collections import defaultdict

import sys
sys.path.append('..')
from functions import find_genes_from_database
from targets import gene_target

agg_funcs = {
    'class': lambda x: '/'.join(x.unique())
}

compound2compound = {'Quaternary Ammonium': 'Quaternary Ammonium Compounds (QACs)'}

def add_amrfinderplus_annotations(file: str, onto: Ontology, logger, db_name: str = 'AMRFinderPlus'):
    """Add AMRFinderPlus annotations to the

    Parameters
    ----------
    file : str
        Path to the file containing AMRFinderPlus annotations
    onto : Ontology
        The ontology object to which annotations will be added
    logger : _type_
        Logger object for logging messages
    db_name : str, optional
        Name of the database in the ontolgoy, by default 'AMRFinderPlus'
    """

    # Load the annotation file
    amrfinderplus_anno = pd.read_csv(file, sep='\t')

    # Replace nan strings with actual NaN values and fill with empty strings
    string_columns = amrfinderplus_anno.select_dtypes(include='object').columns
    amrfinderplus_anno[string_columns] = amrfinderplus_anno[string_columns].replace('nan', np.nan).fillna('')

    # Find the genes from the specified database in the ontology
    matched_genes = find_genes_from_database(onto, database_name=db_name)

    # Lists for storing failed matches (for logging purposes)
    failed_matches = []
    failed_ab_matches = defaultdict(list)

    # Loop through matched pan_ genes and original gene names to add the annotations
    for gene, og in matched_genes.items():

        # Extract the fasta header for the current gene
        fasta_header = [fh for fh in og.original_fasta_header if db_name in fh][0]

        # Match the annotation data based on the refseq protein accession
        m = amrfinderplus_anno.loc[amrfinderplus_anno["refseq_protein_accession"] == fasta_header.split('|')[1]]
        m = m.groupby("refseq_protein_accession").agg(agg_funcs).reset_index()

        # Log if no matches were foiund
        if m.shape[0] == 0:
            failed_matches.append(f"{gene.name} ({og.name})")
            continue

        # Split the antibiotic classes and add those annotations to the gene instances
        ab_classes = m['class'].item().split('/')
        for ab_class in ab_classes:
            ab_class = compound2compound.get(ab_class.title(), ab_class.title())
            success_match = gene_target(gene, og, target=ab_class, onto=onto, db_name=db_name)

            # log if no successful match
            if not success_match:
                failed_ab_matches[ab_class] = f"{gene.name} ({og.name})"

    # Output logging messages for any failed matches
    if len(failed_matches) > 0:
        logger.error(f"{db_name}: Failed to find the genes in the annotation file ({file}): {', '.join(failed_matches)}")    
    
    if len(failed_ab_matches) > 0:
        failed_ab_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k, v in failed_ab_matches.items()])
        logger.warning(f"{db_name}: Failed to find the annotations for the following antibiotics:\n" + failed_ab_matches_string)

    # Log the successful addition of annotations
    logger.success(f"Added {db_name} annotations to the PanRes ontology.")

