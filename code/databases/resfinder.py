from typing import Type
import pandas as pd
from owlready2 import *
import numpy as np
from collections import defaultdict
import re

import sys
sys.path.append('..')
from functions import find_genes_from_database
from targets import gene_target

agg_funcs = {
    'Class': lambda x: ', '.join(x.unique()),
    'Phenotype': lambda x: ', '.join(x.unique()),
    'PMID': 'first',
    'Mechanism of resistance': lambda x: ', '.join(x.unique()),
}

class2class = {'Ionophores': 'Ionophore'}

def add_resfinder_annotations(file: str, onto: Ontology, logger, db_name: str = 'ResFinder'):
    """Add ResFinder annotations to the ontology.

    Parameters
    ----------
    file : str
        Path to the file containing ResFinder annotations
    onto : Ontology
        The ontology object to which annotations will be added
    logger : loguru.logger
        Logger object for logging messages.
    db_name : str, optional
        Name of the database, by default 'ResFinder'
    """

    # Load the annotation file
    resfinder_annotations = pd.read_csv(file, sep='\t')

    # Clean strings
    resfinder_annotations['Gene_accession no.'] = resfinder_annotations['Gene_accession no.'].str.replace("'", "")
    string_columns = resfinder_annotations.select_dtypes(include='object').columns
    resfinder_annotations[string_columns] = resfinder_annotations[string_columns].replace(['nan'], np.nan).fillna('')

    # Lists for storing failed matches
    failed_matches = [] 
    failed_class_matches = defaultdict(list)
    failed_phenotype_matches = defaultdict(list)

    # Find genes from the specified database in the ontology
    matched_genes = find_genes_from_database(onto, database_name=db_name)

    # Loop through each matched gene and add annotations
    for gene, og in matched_genes.items():
        # Match the annotation data based on the gene accession number
        m = resfinder_annotations.loc[resfinder_annotations['Gene_accession no.'] == og.name]
        m = m.groupby('Gene_accession no.').agg(agg_funcs).reset_index()

        # Log failed matches
        if m.shape[0] == 0:
            failed_matches.append(f"{gene.name} ({og.name})")
            continue
        
        # Get class annotations and clean them a bit
        ab_classes = m['Class'].item().replace(" Unknown", "").title()

        for ab_class in ab_classes.split(','):
            ab_class = ab_class.strip()
            success_match = gene_target(gene, og, target=class2class.get(ab_class,ab_class), onto=onto, db_name=db_name)
            # log if failed match
            if not success_match:
                failed_class_matches[ab_class].append(f"{gene.name} ({og.name})")

        # Get phenotype annotations and clean them
        phenotypes = m['Phenotype'].item().split(',') 
        for phenotype in set(phenotypes):
            phenotype = phenotype.replace('Unknown', '').strip().title()
            phenotype = re.sub(r"s$", "", phenotype)
            success_match = gene_target(gene, og, target=phenotype, onto=onto, db_name=db_name)

            # log if failed match
            if not success_match:
                failed_phenotype_matches[phenotype].append(f"{gene.name} ({og.name})")
        
        # # Add DNA accession
        dna_acc = og.name.split('_')[-1].replace(f"|{db_name}", "")
        gene.accession.append(dna_acc)
        og.accession.append(dna_acc)
    
    # Output logging messages for any failed matches
    if len(failed_matches) > 0:
        logger.error(f"{db_name}: Failed to find the genes in the annotation file ({file}): {', '.join(failed_matches)}")    
    
    if len(failed_class_matches) > 0:
        failed_class_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k,v in failed_class_matches.items()])
        logger.warning(f"{db_name}: Failed to find classes for the following annotations:\n" + failed_class_matches_string)

    if len(failed_phenotype_matches) > 0:
        failed_phenotype_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k,v in failed_phenotype_matches.items()])
        logger.warning(f"{db_name}: Failed to find phenotypes for the following annotations:\n" + failed_phenotype_matches_string)
    
    # Log the successful addition of annotations
    logger.success(f"Added {db_name} annotations to the PanRes ontology.")
