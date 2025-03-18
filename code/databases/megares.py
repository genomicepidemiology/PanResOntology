import re
import pandas as pd
from owlready2 import *
from collections import defaultdict

import sys
sys.path.append('..')
from functions import find_genes_from_database, get_instance
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

extra_approved_databases = [
     'BacMet'
]

def add_megares_annotations(onto: Ontology, mappingfile: str, logger, db_name: str = 'MegaRes'):

    matched_genes = find_genes_from_database(onto, database_name=db_name)

    megares_mappings = pd.read_csv(mappingfile, on_bad_lines='warn', sep=';')

    failed_class_matches = defaultdict(list)
    failed_type_matches = defaultdict(list)

    p = re.compile(r"_resistance|\sresistance|_|s$")

    for gene, og in matched_genes.items():
        fasta_header = og.original_fasta_header[0].replace(f"|{db_name}", "")
        split_fasta_header = fasta_header.split('|')

        resistance_type = split_fasta_header[1]
        resistance_classes = split_fasta_header[2]

        # resistance_classes = sub(r's$', '', resistance_classes.replace('_resistance', '').replace('_',' ')).title()
        resistance_classes = p.sub(lambda m: ' ' if m.group(0) == '_' else '', resistance_classes.lower()).title()
        resistance_classes = megares2class.get(resistance_classes, resistance_classes)

        for resistance_class in resistance_classes.split('/'):
            if resistance_class in to_skip: 
                continue
            success_match = gene_target(gene, og, target=resistance_class, onto=onto, db_name=db_name)
            if not success_match:
                    failed_type_matches[resistance_class].append(f"{gene.name} ({og.name})")
            # else:
                # failed_class_matches[resistance_class].append(f"{gene.name} ({og.name})")

        # get where MegaRes has pulled genes from
        megares_original = megares_mappings.loc[megares_mappings['MEGARes_header'].fillna('').str.startswith(og.name + '|')]
        if megares_original.shape[0]  == 1 and megares_original['Database'].item() in extra_approved_databases:
            megares_original_db_instance = get_instance(onto = onto, name = megares_original['Database'].item())
            og.is_from_database.append(megares_original_db_instance)
            gene.is_from_database.append(megares_original_db_instance)
        
    if len(failed_class_matches) > 0: 
        failed_class_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k,v in failed_class_matches.items()])
        logger.warning(f"{db_name}: Did not match the following to any type of resistance annotation:\n" + failed_class_matches_string)
    
    if len(failed_type_matches) > 0: 
        failed_type_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k,v in failed_type_matches.items()])
        logger.warning(f"{db_name}: Did not figure out what type of annotations for:\n" + failed_type_matches_string)
    
    logger.success(f"Added {db_name} annotations to the PanRes ontology.")
