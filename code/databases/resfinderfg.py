import pandas as pd
from owlready2 import *
from collections import defaultdict

import sys
sys.path.append('..')
from functions import find_genes_from_database
from targets import gene_target

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
# CEF: pan_11851 (KM113772.1),pan_7473 (KM113767.1),pan_8469 (KM113770.1),pan_13710 (KM113771.1),pan_9978 (KM113768.1),pan_13450 (KM113773.1)
#   - KM113767, KM113768, KM113770, KM113771, KM113772, KM113773 - all from the same study PUBMED 25288759.
#       - KM113769 not in ResFinderFG?
# SDX: pan_10206 (KS307999.1),pan_9587 (KS307948.1), pan_12812 (KS307981.1)
#   - SUL: sulfadimethoxine 
# STX: pan_11836 (MK159021.1),pan_13458 (MK159019.1)
# SPT: pan_11835 (JN625762.1),pan_13997 (MN340016.1)
#   -  spectinomycin
# MEM: pan_10205 (KY705336.1),pan_10208 (KY705337.1)
#   - meropenem
# TZB: pan_12801 (KX125600.1)
#   - Tazobactam
# STR: pan_9976 (JX875573.1)
# ERY: pan_10870 (KY705334.1)

fg2class = {
    'SMZ': 'Sulfamethazine',
    'AMP': 'Ampicillin',
    'AZM': 'Azithromycin',
    'CHL': 'Chloramphenicol',
    'AMC': 'Amoxicillin+Clavulanic acid',
    'KAN': 'Kanamycin',
    'TZP': 'Piperacillin+Tazobactam',
    'CEP': 'Cephalothin',
    'SDX': 'Sulfadimethoxine',
    'SPT': 'Spectinomycin',
    'MEM': 'Meropenem',
    'TZB': 'Tazobactam',
    'ERY': 'Erythromycin' 

}

def add_resfinderfg_annotations(file: str, onto: Ontology, logger, db_name: str = 'ResFinderFG'):
    """Add ResFinderFG annotations to the ontology.

    Parameters
    ----------
    file : str
        Path to the file containing ResFinderFG annotations
    onto : Ontology
        The ontology object to which annotations will be added
    logger : loguru.logger
        Logger object for logging messages.
    db_name : str, optional
        Name of the database, by default 'ResFinderFG'
    """
    
    # Load the annotations by parsing the text file
    with open(file, 'r') as f:
        lines = [l.strip() for l in f.readlines()]
    for l in lines:
        ll = l.split(':')
        if len(ll) == 2:
            fg2class[ll[0].strip()]=ll[1].strip()

    # Lists for storing failed acronym matches
    failed_acronym_matches = defaultdict(list)
    
    # Find genes from the specified database in the ontology
    matched_genes = find_genes_from_database(onto, database_name=db_name)
    
    # Loop through each matched gene and add annotations
    for gene, og in matched_genes.items():
        # Extract and clean the fasta header for the current gene
        fasta_header = og.original_fasta_header[0].replace(f"|{db_name}", "")
        ab_class_acronym = fasta_header.split('|')[-1]
        ab_class_name = fg2class.get(ab_class_acronym, ab_class_acronym).title()

        # Add resistance links
        success_match = gene_target(gene, og, target=ab_class_name, onto=onto, db_name=db_name)
        # log if no match
        if not success_match:
            failed_acronym_matches[ab_class_acronym].append(f"{gene.name} ({og.name})")

        # Add accession numbers
        dna_acc = fasta_header.split('|')[1]
        gene.accession.append(dna_acc)
        og.accession.append(dna_acc)
    
    # Output logging messages for any failed matches
    if len(failed_acronym_matches) > 0:
        failed_acronym_matches_string = "\n".join([f"{k}: {','.join(v)}" for k,v in failed_acronym_matches.items()])
        logger.warning(f"{db_name}: Failed to find acronym translations for:\n" + failed_acronym_matches_string)
    
    # Log the successful addition of annotations
    logger.success(f"Added {db_name} annotations to the PanRes ontology.")
