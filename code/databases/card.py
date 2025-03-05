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


def add_card_annotations(file, onto, logger):

    card_annotations = pd.read_csv(file, sep='\t')

    p = re.compile(r"(ARO\:\d+)")

    matched_genes = find_genes_from_database(onto, database_name='CARD')

    failed_matches = []
    failed_phenotype_matches = defaultdict(list)

    for gene, og in matched_genes.items():
        
        fasta_header = og.original_fasta_header[0]
        regex_match = p.findall(fasta_header)
        if regex_match:
            aro_number = regex_match[0]
            m = card_annotations.loc[card_annotations["ARO Accession"] == aro_number]
            
            if m.shape[0] == 0:
                # logger.error(f"CARD: no match found for {gene.name} ({og.name})")
                failed_matches.append(f"{gene.name} ({og.name}, {aro_number})")
                continue

            m = m.groupby("ARO Accession").agg(agg_funcs).reset_index()

            # like = predicted phenotype
            phenotypes = m['Drug Class'].item().replace('-like', '').split(';')

            for phenotype in phenotypes:
                phenotype = phenotype.strip().replace(' antibiotic', '').title()
                success_match = gene_target(gene, og, target=phenotype, onto=onto)
                if not success_match:
                    failed_phenotype_matches[phenotype].append(f"{gene.name} ({og.name})")

            # Add DNA accession number
            dna_accessions = m['DNA Accession'].item().split(';')
            for dna_acc in dna_accessions:
                gene.dna_accession.append(dna_acc)

            # Add protein accession number

            # Add link to CARD information
            card_url=f"https://card.mcmaster.ca/ontology/{m['CVTERM ID'].item()}"
            gene.card_link.append(card_url)
    
    if len(failed_matches) > 0:
        logger.error(f"CARD: Failed to find the genes in the annotation file ({file}): {', '.join(failed_matches)}")    

    if len(failed_phenotype_matches):
        failed_phenotype_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k,v in failed_phenotype_matches.items()])
        logger.warning("CARD: Could not match phenotypes annotations for the following:\n" + failed_phenotype_matches_string)

    logger.success("Added CARD annotations to the PanRes ontology.")
