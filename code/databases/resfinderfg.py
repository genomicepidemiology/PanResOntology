import pandas as pd
from owlready2 import *
from collections import defaultdict

import sys
sys.path.append('..')
from functions import check_targets, find_genes_from_database, find_target_annotation, get_or_create_instance

# Missing acronyms:
# SMZ: pan_19 (KY705325.1) 
#   - KY705325.1 = clone SCT-SMZ-28 
#       - SMZ : sulfamethazine
# AMP: pan_6 (MG586042.1),pan_2195 (KU606673.1),pan_2 (KU606237.1),pan_7 (MG586043.1),pan_17 (MG586023.1),pan_1693 (KU606242.1),pan_20 (KU606773.1),pan_15 (KU606489.1)
#   - MG586042.1 = clone AMP10_S_DS_C2
#       - AMP : ampicillin
#  - KU606673.1:  1-E2_AP_Contig_9
# AZM: pan_23 (MG585948.1)
#   - MG585948.1 = clone AZI1_F_WW_C1 
#       - AZM : Azithromycin
# CHL: pan_1875 (KU544508.1),pan_451 (KU544029.1)
#   - KU544508.1 = clone CH_OUT_B_Contig_4
#        - CHL: chloramphenicol ?
# AMC: pan_2299 (KU605848.1),pan_5 (KU605846.1),pan_1204 (KU607243.1)
#   - KU605848.1 = clone  5-B1_AXCL_Contig_21 
#       - AXCL: Amoxicillin + Clavulanate (Clavulanic acid)
# KAN: pan_2197 (KY705354.1)
#   - KY705354.1 = clone UTC-KAN-K1-CONTIG-B

fg2class = {
    'SMZ': 'Sulfamethazine',
    'AMP': 'Ampicillin',
    'AZM': 'Azithromycin',
    'CHL': 'Chloramphenicol',
    'AMC': 'Amoxicillin+Clavulanic acid',
    'KAN': 'Kanamycin'
}

def add_resfinderfg_annotations(file, onto, targetfile, logger):
    
    targets = check_targets(targetfile)

    with open(file, 'r') as f:
        lines = [l.strip() for l in f.readlines()]
    for l in lines:
        ll = l.split(':')
        if len(ll) == 2:
            fg2class[ll[0].strip()]=ll[1].strip()

    failed_acronym_matches = defaultdict(list)

    matched_genes = find_genes_from_database(onto, database_name='ResFinderFG')
    for gene, og in matched_genes.items():
        fasta_header = og.original_fasta_header[0].replace("|ResFinderFG", "")
        ab_class_acronym = fasta_header.split('|')[-1]
        ab_class_name = fg2class.get(ab_class_acronym, ab_class_acronym).title()

        dna_acc = fasta_header.split('|')[1]
        gene.dna_accession.append(dna_acc)
        
        ab_class_match = find_target_annotation(target=ab_class_name, annotated_targets=targets)
        
        if '+' in ab_class_name or ab_class_match is not None:
            ab_class_instance = get_or_create_instance(onto, onto.AntibioticResistanceClass, ab_class_name)
            gene.has_resistance_class.append(ab_class_instance)
            
        else:
            failed_acronym_matches[ab_class_acronym].append(f"{gene.name} ({og.name})")
    
    if len(failed_acronym_matches) > 0:
        failed_acronym_matches_string = "\n".join([f"{k}: {','.join(v)}" for k,v in failed_acronym_matches.items()])
        logger.warning("ResFinderFG: Failed to find acronym translations for:\n" + failed_acronym_matches_string)
    
    logger.success("Added ResFinderFG annotations to the PanRes ontology.")
