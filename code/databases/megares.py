from multiprocessing import Value
from re import split, sub
import pandas as pd
from owlready2 import *
from collections import defaultdict

import sys
sys.path.append('..')
from functions import check_targets, find_genes_from_database, find_target_annotation, get_or_create_instance

megares2class = {
    'Betalactam': 'Beta-Lactam',
    'Mls': 'Macrolide/Lincosamide/Streptogramin B',
    'Aluminum': 'Aluminium'
}

to_skip = [
    'Multi-Drug', 
    'Drug And Biocide', 
    'Multi-Metal', 
    'Drug And Biocide And Metal',
    'Multi-Biocide',
    'Biocide And Metal',
    'Drug And Metal'
]

def add_megares_annotations(onto, targetfile, logger):
    targets = check_targets(targetfile)


    matched_genes = find_genes_from_database(onto, database_name='MegaRes')

    failed_class_matches = defaultdict(list)
    failed_type_matches = defaultdict(list)

    p = re.compile(r"_resistance|\sresistance|_|s$")

    for gene, og in matched_genes.items():
        fasta_header = og.original_fasta_header[0].replace("|MegaRes", "")
        split_fasta_header = fasta_header.split('|')

        resistance_type = split_fasta_header[1]
        resistance_classes = split_fasta_header[2]

        # resistance_classes = sub(r's$', '', resistance_classes.replace('_resistance', '').replace('_',' ')).title()
        resistance_classes = p.sub(lambda m: ' ' if m.group(0) == '_' else '', resistance_classes.lower()).title()
        resistance_classes = megares2class.get(resistance_classes, resistance_classes)

        for resistance_class in resistance_classes.split('/'):
            if resistance_class in to_skip: 
                continue
            class_match = find_target_annotation(target=resistance_class, annotated_targets=targets)
            if class_match:
                if resistance_type == 'Drugs':
                    if class_match[0] == 'antibiotic_class':
                        ab_class_instance = get_or_create_instance(onto, onto.AntibioticResistanceClass, resistance_class)
                        gene.has_resistance_class.append(ab_class_instance)
                    elif class_match[0] == 'antibiotic':
                        ab_phenotype_instance = get_or_create_instance(onto, onto.AntibioticResistancePhenotype, resistance_class)
                        gene.has_predicted_phenotype.append(ab_phenotype_instance)

                        # class
                        ab_class_instance = get_or_create_instance(onto, onto.AntibioticResistanceClass, class_match[1]['class'].unique()[0])
                        gene.has_resistance_class.append(ab_class_instance)

                elif resistance_type == 'Metals' and class_match[0] == 'metal':
                    metal_instance = get_or_create_instance(onto, onto.Metal, resistance_class)
                    gene.has_predicted_metal_resistance.append(metal_instance)
                elif resistance_type == 'Biocides' and class_match[0] == 'biocide':
                    biocide_instance = get_or_create_instance(onto, onto.Biocide, resistance_class)
                    gene.has_predicted_biocide_resistance.append(biocide_instance)

                else:
                    # logger.error(f"MegaRes: Could not match the type of {resistance_class} to any known annotation for {gene.name} ({og.name}).")
                    failed_type_matches[resistance_class].append(f"{gene.name} ({og.name})")
            else:
                # logger.error(f"MegaRes: Could not figure out what '{resistance_class}' matches ({gene.name} ({og.name})).")
                failed_class_matches[resistance_class].append(f"{gene.name} ({og.name})")
        
    if len(failed_class_matches) > 0: 
        failed_class_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k,v in failed_class_matches.items()])
        logger.warning("MegaRes: Did not match the following to any type of resistance annotation:\n" + failed_class_matches_string)
    
    if len(failed_type_matches) > 0: 
        failed_type_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k,v in failed_type_matches.items()])
        logger.warning("MegaRes: Did not figure out what type of annotations for:\n" + failed_type_matches_string)
    
    logger.success("Added MegaRes annotations to the PanRes ontology.")
