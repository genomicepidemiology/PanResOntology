import pandas as pd
from owlready2 import *
import numpy as np
from collections import defaultdict

import sys
sys.path.append('..')
from functions import find_genes_from_database
from targets import gene_target

agg_funcs = {
    'class': lambda x: '/'.join(x.unique())
}

def add_amrfinderplus_annotations(file: str, onto: Ontology, logger, db_name: str = 'AMRFinderPlus'):
    amrfinderplus_anno = pd.read_csv(file, sep='\t')
    string_columns = amrfinderplus_anno.select_dtypes(include='object').columns
    amrfinderplus_anno[string_columns] = amrfinderplus_anno[string_columns].replace('nan', np.nan).fillna('')

    # amr_classes = amrfinderplus_anno['class'].dropna().unique().tolist()
    matched_genes = find_genes_from_database(onto, database_name=db_name)

    failed_matches = []
    failed_ab_matches = defaultdict(list)

    for gene, og in matched_genes.items():
        fasta_header = [fh for fh in og.original_fasta_header if db_name in fh][0]
        m = amrfinderplus_anno.loc[amrfinderplus_anno["refseq_protein_accession"] == fasta_header.split('|')[1]]
        m = m.groupby("refseq_protein_accession").agg(agg_funcs).reset_index()

        if m.shape[0] == 0:
            # logger.error(f"AMRFinderPlus: no match found for {gene.name} ({og.name})")
            failed_matches.append(f"{gene.name} ({og.name})")
            continue

        ab_classes = m['class'].item().split('/')
        for ab_class in ab_classes:
            ab_class = ab_class.title()
            success_match = gene_target(gene, og, target=ab_class, onto=onto, db_name=db_name)
            if not success_match:
                failed_ab_matches[ab_class] = f"{gene.name} ({og.name})"

    if len(failed_matches) > 0:
        logger.error(f"{db_name}: Failed to find the genes in the annotation file ({file}): {', '.join(failed_matches)}")    
    
    if len(failed_ab_matches) > 0:
        failed_ab_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k, v in failed_ab_matches.items()])
        logger.warning(f"{db_name}: Failed to find the annotations for the following antibiotics:\n" + failed_ab_matches_string)

    logger.success(f"Added {db_name} annotations to the PanRes ontology.")

