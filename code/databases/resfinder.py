from typing import Type
import pandas as pd
from owlready2 import *
import numpy as np
from collections import defaultdict
import re

import sys
sys.path.append('..')
from functions import check_targets, find_genes_from_database, find_target_annotation, get_or_create_instance

agg_funcs = {
    'Class': lambda x: ', '.join(x.unique()),
    'Phenotype': lambda x: ', '.join(x.unique()),
    'PMID': 'first',
    'Mechanism of resistance': lambda x: ', '.join(x.unique()),
}

def add_resfinder_annotations(file, onto, targetfile, logger):

    targets = check_targets(targetfile)

    resfinder_annotations = pd.read_csv(file, sep='\t')
    resfinder_annotations['Gene_accession no.'] = resfinder_annotations['Gene_accession no.'].str.replace("'", "")

    string_columns = resfinder_annotations.select_dtypes(include='object').columns
    resfinder_annotations[string_columns] = resfinder_annotations[string_columns].replace(['nan'], np.nan).fillna('')

    failed_matches = [] 
    failed_class_matches = defaultdict(list)
    failed_phenotype_matches = defaultdict(list)


    matched_genes = find_genes_from_database(onto, database_name='ResFinder')
    for gene, og in matched_genes.items():
        m = resfinder_annotations.loc[resfinder_annotations['Gene_accession no.'] == og.name]
        m = m.groupby('Gene_accession no.').agg(agg_funcs).reset_index()

        if m.shape[0] == 0:
            # logger.error(f"ResFinder: no match found for {gene.name} ({og.name})")
            failed_matches.append(f"{gene.name} ({og.name})")
            continue


        ab_classes = m['Class'].item().replace(" Unknown", "").title()

        for ab_class in ab_classes.split(','):
            ab_class = ab_class.strip()
            ab_class_match = find_target_annotation(target=ab_class, annotated_targets=targets)
        
            if ab_class_match is not None:
                ab_class_instance = get_or_create_instance(onto, onto.AntibioticResistanceClass, ab_class)
                gene.has_resistance_class.append(ab_class_instance)
            else:
                failed_class_matches[ab_class].append(f"{gene.name} ({og.name})")
                # logger.warning(f"ResFinder: No class match found for '{ab_class}' for {gene.name} ({og.name}) ")
            
        phenotypes = m['Phenotype'].item().split(',') #TODO: What does plus mean here?
        for phenotype in set(phenotypes):
            phenotype = phenotype.replace('Unknown', '').strip().title()
            phenotype = re.sub(r"s$", "", phenotype)
            ab_phenotype_match = find_target_annotation(target=phenotype, annotated_targets=targets)
            if ab_phenotype_match is not None:
                if ab_phenotype_match[0] == 'antibiotic':
                    phenotype_class = ab_phenotype_match[1]['class'].unique()[0]
                    ab_class_instance = get_or_create_instance(onto, onto.AntibioticResistanceClass, phenotype_class)

                    gene.has_resistance_class.append(ab_class_instance)

                    if phenotype != phenotype_class:
                        ab_phenotype_instance = get_or_create_instance(onto, onto.AntibioticResistancePhenotype, phenotype)            
                        ab_phenotype_instance.phenotype_is_class.append(ab_class_instance)
                        gene.has_predicted_phenotype.append(ab_phenotype_instance)
                elif ab_phenotype_match[0] == 'antibiotic_class':
                    ab_class_instance = get_or_create_instance(onto, onto.AntibioticResistanceClass, phenotype)
                    gene.has_resistance_class.append(ab_class_instance)
            elif '+' in phenotype:
                ab_phenotype_instance = get_or_create_instance(onto, onto.AntibioticResistancePhenotype, phenotype)
                ab_phenotype_instance.is_drug_combination.append(True)
                gene.has_predicted_phenotype.append(ab_phenotype_instance)
            else:
                failed_phenotype_matches[phenotype].append(f"{gene.name} ({og.name})")
                # logger.warning(f"ResFinder: No phenotype match found for '{phenotype}' {gene.name} ({og.name})")
        
        # Add DNA accession
        dna_acc = og.name.split('_')[-1].replace("|ResFinder", "")
        gene.dna_accession.append(dna_acc)

        # Add protein accession ?
        protein_acc = og.name.split('_')[-1].replace("|ResFinder", "")

    if len(failed_matches) > 0:
        logger.error(f"ResFinder: Failed to find the genes in the annotation file ({file}): {', '.join(failed_matches)}")    
    
    if len(failed_class_matches) > 0:
        failed_class_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k,v in failed_class_matches.items()])
    # logger.warning(f"ResFinder: Failed to find classes for the following annotations:\n{'\n'.join(failed_class_matches_string)}")
        logger.warning("ResFinder: Failed to find classes for the following annotations:\n" + failed_class_matches_string)

    if len(failed_phenotype_matches) > 0:
        failed_phenotype_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k,v in failed_phenotype_matches.items()])
        logger.warning("ResFinder: Failed to find phenotypes for the following annotations:\n" + failed_phenotype_matches_string)

    logger.success("Added ResFinder annotations to the PanRes ontology.")
