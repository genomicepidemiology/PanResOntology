import pandas as pd
from owlready2 import *
import re
import numpy as np
from collections import defaultdict

import sys
sys.path.append('..')
from functions import find_genes_from_database
from targets import gene_target

acr2class = {
    'AGly': 'Aminoglycoside',
    'Bla': 'Beta-Lactam',
    'Fos': 'Fosfomycin',
    'Flq': 'Fluoroquinolone',
    'Gly': 'Glycopeptide',
    'MLS': 'Macrolide/Lincosamide/Streptogramin B',
    'Phe': 'Phenicol',
    'Rif': 'Rifampin',
    'Sul': 'Sulfonamide',
    'Tet': 'Tetracycline',
    'Tmt': 'Trimethoprim',
    'Col': 'Colistin'
}

def add_argannot_annotations(onto, logger):
    matched_genes = find_genes_from_database(onto, database_name='ARGANNOT')

    failed_regex_matches = list()
    failed_ab_matches = defaultdict(list)

    p = re.compile(r"\|\((\w{3,4})\)")

    for gene, og in matched_genes.items():
        fasta_header = [fh for fh in og.original_fasta_header if 'ARGANNOT' in fh][0]
        m = p.findall(fasta_header)
        ab_acronym = m[0] if m else None
        if ab_acronym is not None and ab_acronym in acr2class.keys():
            ab = acr2class[ab_acronym].title()
            success_match = gene_target(gene, og, target=ab, onto=onto)
            if not success_match:
                failed_ab_matches[ab].append(f"{gene.name} ({og.name})")
        else:
            failed_regex_matches.append(f"{gene.name} ({og.name})")
        
    
    if len(failed_regex_matches) > 0:
        logger.warning(f"ARGANNOT: Failed to extract class acronyms: {', '.join(failed_regex_matches)}")    

    if len(failed_ab_matches) > 0:
        failed_ab_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k, v in failed_ab_matches.items()])
        logger.warning("ARGANNOT: Failed to find target annotations for:\n" + failed_ab_matches_string)
    
    logger.success("Added ARG-ANNOT annotations to the PanRes ontology.")
