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
    'antibiotic': lambda x: "/".join(x)
}

acr2ab = {
    'APS': 'Apramycin',
    'CEF': 'Ceftobiprole',
    'CPR': 'Ciprofloxacin',
    'DEL': 'Delafloxacin',
    'DOX': 'Doxycycline',
    'ERA': 'Eravacycline',
    'FEP': 'Cefepime',
    'FID': 'Cefiderocol',
    'GEN': 'Gentamicin',
    'GPC': 'Gentamicin',
    'MER': 'Meropenem',
    'MOX': 'Moxifloxacin',
    'MOXa': 'Moxifloxacin',
    'OMA': 'Omadacycline',
    'POL': 'POL-7306',
    'PXB': 'Polymyxin-B',
    'SCH': 'SCH79797',
    'SPR': 'SPR-206',
    'SUL': 'Sulopenem',
    'TCS': 'Triclosan',
    'ZOL': 'Zoliflodacin',
}

def add_csabapal_annotations(file: str, onto: Ontology, logger, db_name: str = 'CsabaPal'):
    anno = pd.read_csv(file, sep=',')

    matched_genes = find_genes_from_database(onto, database_name=db_name)

    failed_matches = defaultdict(list)
    failed_ab_matches = defaultdict(list)

    for gene, og in matched_genes.items():
        m = anno.loc[anno['orf_unique'] == og.name,]
        m = m.groupby('orf_unique').agg(agg_funcs).reset_index()
        
        for ab_short in m['antibiotic'].item().split('/'):
            ab = acr2ab.get(ab_short, None)
            if ab is None:
                failed_matches[ab].append(f"{gene.name} ({og.name})")
                continue
            success_match = gene_target(gene, og, target=ab, onto = onto, db_name = db_name)
            if not success_match:
                failed_ab_matches[ab].append(f"{gene.name} ({og.name})")
            
    if len(failed_matches) > 0:
        logger.error(f"{db_name}: Failed to find acronyms translations for: {', '.join(failed_matches)}")    
    
    if len(failed_ab_matches) > 0:
        failed_ab_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k, v in failed_ab_matches.items()])
        logger.warning(f"{db_name}: Failed to find the annotations for the following antibiotics:\n" + failed_ab_matches_string)

    logger.success(f"Added {db_name} annotations to the PanRes ontology.")