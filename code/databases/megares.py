import re
import pandas as pd
from owlready2 import *
from collections import defaultdict

import sys
sys.path.append('..')
from functions import find_genes_from_database
from targets import gene_target

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

def add_megares_annotations(onto, logger):

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
            success_match = gene_target(gene, og, target=resistance_class, onto=onto)
            if not success_match:
                    failed_type_matches[resistance_class].append(f"{gene.name} ({og.name})")
            # else:
                # failed_class_matches[resistance_class].append(f"{gene.name} ({og.name})")
        
    if len(failed_class_matches) > 0: 
        failed_class_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k,v in failed_class_matches.items()])
        logger.warning("MegaRes: Did not match the following to any type of resistance annotation:\n" + failed_class_matches_string)
    
    if len(failed_type_matches) > 0: 
        failed_type_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k,v in failed_type_matches.items()])
        logger.warning("MegaRes: Did not figure out what type of annotations for:\n" + failed_type_matches_string)
    
    logger.success("Added MegaRes annotations to the PanRes ontology.")
