from multiprocessing import Value
import re
import pandas as pd
from owlready2 import *
from collections import defaultdict

import sys
sys.path.append('..')
from functions import find_genes_from_database
from targets import gene_target

metal2metal = {
    'Zink': 'Zinc',
    'Aluminum': 'Aluminium'
}

def add_metalres_annotations(onto, logger):
    matched_genes = find_genes_from_database(onto, database_name='MetalRes')

    p = re.compile(r"\s(\S+)\sresistance")

    failed_class_matches = defaultdict(list)
    failed_type_matches = defaultdict(list)

    for gene, og in matched_genes.items():
        fasta_header = [fh for fh in og.original_fasta_header if fh.endswith('MetalRes')][0]
        m = p.findall(fasta_header)
        if len(m) == 1:
            metals = m[0].split('/')
            for metal in metals:
                metal = metal2metal.get(metal, metal)
                success_match = gene_target(gene=gene, og=og, target=metal, onto=onto)
                if not success_match:
                    failed_class_matches[metal].append(f"{gene.name} ({og.name})")
        else:
            failed_type_matches[metal].append(f"{gene.name} ({og.name})")
        
        # Add accession
        accession = og.name.split('_')[-1]
        gene.protein_accession.append(accession)
    

    if len(failed_class_matches) > 0: 
        failed_class_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k,v in failed_class_matches.items()])
        logger.warning("MetalRes: Did not match the following to any type of resistance annotation:\n" + failed_class_matches_string)
    
    if len(failed_type_matches) > 0: 
        failed_type_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k,v in failed_type_matches.items()])
        logger.warning("MetalRes: Did not extract metals for:\n" + failed_type_matches_string)
    
    logger.success("Added MetalRes annotations to the PanRes ontology.")
