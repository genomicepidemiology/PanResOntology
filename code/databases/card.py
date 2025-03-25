import pandas as pd
from owlready2 import *
import re
import numpy as np
from collections import defaultdict

import sys
sys.path.append('..')
from functions import find_genes_from_database
from targets import gene_target

agg_funcs = {
    'Drug Class': lambda x: ';'.join(x.unique()),
    'DNA Accession': lambda x: ';'.join(x.unique()),
    'Protein Accession': lambda x: ';'.join(x.unique()),
    'CVTERM ID': 'first'
}


def add_card_annotations(file: str, onto: Ontology, logger, db_name: str = 'CARD'):
    """Add CARD annotations to the ontology.

    Parameters
    ----------
    file : str
        Path to the file containing CARD annotations
    onto : Ontology
        The loaded ontology to which annotations will be added
    logger : loguru.logger
        Logger object for logging messages
    db_name : str, optional
        Name of the database, by default 'CARD'
    """

    # load the CARD annotation file
    card_annotations = pd.read_csv(file, sep='\t')

    # Compile a regular expression to match ARO numbers in the fasta headers
    p = re.compile(r"(ARO\:\d+)")

    # Find genes from the specified database in the ontology
    matched_genes = find_genes_from_database(onto, database_name=db_name)

    # Lists for logging failed matches
    failed_matches = []
    failed_phenotype_matches = defaultdict(list)

    # Loop through matched pan_ genes and original gene names to add the annotations
    for gene, og in matched_genes.items():

        # Extract the fasta header and extract ARO number
        fasta_header = og.original_fasta_header[0]
        regex_match = p.findall(fasta_header)

        # check if there was a match
        if regex_match:
            # Match the annotation data based on the ARO number
            aro_number = regex_match[0]
            m = card_annotations.loc[card_annotations["ARO Accession"] == aro_number]
            
            if m.shape[0] == 0:
                # Log if no matches were found
                failed_matches.append(f"{gene.name} ({og.name}, {aro_number})")
                continue
            
            # Aggregate the metadata
            m = m.groupby("ARO Accession").agg(agg_funcs).reset_index()

            # Extract the phenotypes and clean them
            phenotypes = m['Drug Class'].item().replace('-like', '').split(';')
            
            # Loop through phenotypes and add them to the ontology
            for phenotype in phenotypes:
                phenotype = phenotype.strip().replace(' antibiotic', '').title()
                success_match = gene_target(gene, og, target=phenotype, onto=onto, db_name=db_name)
                if not success_match:
                    failed_phenotype_matches[phenotype].append(f"{gene.name} ({og.name})")

            # Add DNA accession number to the pan_ gene and the original gene name
            dna_accessions = m['DNA Accession'].item().split(';')
            for dna_acc in dna_accessions:
                gene.accession.append(dna_acc)
                og.accession.append(dna_acc)

            # Add link to CARD information
            card_url=f"https://card.mcmaster.ca/ontology/{m['CVTERM ID'].item()}"
            gene.card_link.append(card_url)
    
    # Output logging messages for any failed matches
    if len(failed_matches) > 0:
        logger.error(f"{db_name}: Failed to find the genes in the annotation file ({file}): {', '.join(failed_matches)}")    

    if len(failed_phenotype_matches):
        failed_phenotype_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k,v in failed_phenotype_matches.items()])
        logger.warning(f"{db_name}: Could not match phenotypes annotations for the following:\n" + failed_phenotype_matches_string)
    
    # Log the successful addition of annotations
    logger.success(f"Added {db_name} annotations to the PanRes ontology.")
