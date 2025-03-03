import pandas as pd
from owlready2 import *
import numpy as np
from collections import defaultdict

import sys
sys.path.append('..')
from functions import check_targets, find_genes_from_database, find_target_annotation, get_or_create_instance


agg_funcs = {
    'class': lambda x: '/'.join(x.unique())
}

def add_amrfinderplus_annotations(file, onto, targetfile, logger):
    targets = check_targets(targetfile)

    amrfinderplus_anno = pd.read_csv(file, sep='\t')
    string_columns = amrfinderplus_anno.select_dtypes(include='object').columns
    amrfinderplus_anno[string_columns] = amrfinderplus_anno[string_columns].replace('nan', np.nan).fillna('')

    # amr_classes = amrfinderplus_anno['class'].dropna().unique().tolist()
    matched_genes = find_genes_from_database(onto, database_name='AMRFinderPlus')

    failed_matches = []
    failed_ab_matches = defaultdict(list)

    for gene, og in matched_genes.items():
        fasta_header = [fh for fh in og.original_fasta_header if 'AMRFinderPlus' in fh][0]
        m = amrfinderplus_anno.loc[amrfinderplus_anno["refseq_protein_accession"] == fasta_header.split('|')[1]]
        m = m.groupby("refseq_protein_accession").agg(agg_funcs).reset_index()

        if m.shape[0] == 0:
            # logger.error(f"AMRFinderPlus: no match found for {gene.name} ({og.name})")
            failed_matches.append(f"{gene.name} ({og.name})")
            continue

        ab_classes = m['class'].item().split('/')
        for ab_class in ab_classes:
            ab_class = ab_class.title()
            matched_class = find_target_annotation(target=ab_class, annotated_targets=targets)
            if matched_class:
                if matched_class[0] == 'antibiotic_class':
                    ab_class_instance = get_or_create_instance(onto, onto.AntibioticResistanceClass, ab_class)
                    gene.has_resistance_class.append(ab_class_instance)
                elif matched_class[0] == 'antibotic':
                    ab_phenotype_instance = get_or_create_instance(onto, onto.AntibioticResistancePhenotype, ab_class)            
                    gene.has_predicted_phenotype.append(ab_phenotype_instance)

                    phenotype_class = matched_class[1]['class'].drop_duplicates().item() 
                    ab_class_instance = get_or_create_instance(onto, onto.AntibioticResistanceClass, phenotype_class)
                    ab_phenotype_instance.phenotype_is_class.append(ab_class_instance)
                    # matched_class_info = 
                    # logger.warning(f"AMRFinderPlus: This seems like an antibiotic, not a antibiotic class: {matched_class[1]}")
                    # ab_phenotype_instance = get_or_create_instance(onto, onto.AntibioticResistancePhenotype
            else:
                # logger.warning(f"AMRFinderPlus: No class match found for {ab_class} {gene.name} ({og.name})")
                failed_ab_matches[ab_class] = f"{gene.name} ({og.name})"

    if len(failed_matches) > 0:
        logger.error(f"AMRFinderPlus: Failed to find the genes in the annotation file ({file}): {', '.join(failed_matches)}")    
    
    if len(failed_ab_matches) > 0:
        failed_ab_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k, v in failed_ab_matches.items()])
        logger.warning("AMRFinderPlus: Failed to find the annotations for the following antibiotics:\n" + failed_ab_matches_string)


    logger.success("Added AMRFinderPlus annotations to the PanRes ontology.")

