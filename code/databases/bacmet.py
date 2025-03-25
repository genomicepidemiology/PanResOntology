import re
import pandas as pd
from owlready2 import *
from collections import defaultdict

import sys
sys.path.append('..')
from functions import get_instance, get_or_create_instance, accession_to_pubmed, find_genes_from_database
from targets import gene_target

agg_funcs = {
    'Compound': lambda x: ",".join(x),
    'Accession': lambda x: ",".join(x)
}

def add_bacmet_annotations(onto: Ontology, mappingfile: str, logger, db_name: str = 'BacMet'):
    """Adds BacMet annotations to the ontology.

    Parameters
    ----------
    onto : Ontology
        The ontology object to which annotations will be added
    mappingfile : str
        Path to the file containing BacMet annotations
    logger : _type_
        Logger object for logging messages  
    db_name : str, optional
        Name of the database in the ontology, by default 'BacMet'
    """

    # Load the annotation file
    exp_annotations = pd.read_csv(mappingfile, sep='\t')

    # Clean up strings
    exp_annotations['gene_name'] = exp_annotations['Gene_name'].str.lower().str.split('/')
    exp_annotations = exp_annotations.explode('gene_name')

    # Lists for storing failed matches
    failed_matches = []
    failed_type_matches = defaultdict(list)
    failed_class_matches = defaultdict(list)

    # Compile regular expressions to match compounds and classes
    p = re.compile(r"^((\w+)(\s\w+)?)\s")
    p_class = re.compile(r"^((\w+)(\s\w+)?)\s(\[class\W+((\w+)(\s\w+)?)\])")

    # Get the database instance from the ontology
    db_instance = get_instance(onto=onto, name=db_name)

    # Find the genes from the specified database in the ontology
    matched_genes = find_genes_from_database(onto = onto, database_name = db_name)

    # Loop through matched pan_ genes and original gene names to add the annotations   
    # Matching to the BacMet names is based on the MegaRes gene names
    for gene, meg_gene in matched_genes.items():
        meg_header = meg_gene.original_fasta_header[0]
        meg_gene_name = meg_header.split('|')[-2].lower()

        # Match the annotation data based on the gene name
        m = exp_annotations.loc[
            (exp_annotations['gene_name'] == meg_gene_name)
        ].groupby(['Gene_name']).agg(agg_funcs).reset_index()

        # Only allowing for one match
        if m.shape[0] == 1:
            # Find compounds and original gene names
            compounds = m['Compound'].item()
            original_name = m['Gene_name'].item()

            # Create or get the original gene name instance
            og = get_or_create_instance(onto = onto, cls=onto.OriginalGene, name=original_name)

            # Annotate where the original gene is from
            og.is_from_database.append(db_instance)
            gene.equivalent_to.append(og)

            # Set metal resistance annotations
            for compound in set(compounds.split(',')):
                # test if its a class or an actual compound
                if 'class' in compound:
                    m_compound = p_class.findall(compound.strip())
                else:
                    m_compound = p.findall(compound.strip())
                
                # Check if match is found and get it
                if len(m_compound) > 0:
                    c = m_compound[0][0]

                    # add the annotation and if it fails, log it
                    if not gene_target(gene=gene, og=og, target=c, onto=onto, db_name=db_name):
                        failed_type_matches[c].append(f"{gene.name} ({og.name})")
                    if len(m_compound[0]) > 3:
                        cc = m_compound[0][4]
                        if not gene_target(gene=gene, og=og, target=cc, onto=onto, db_name=db_name):
                            failed_class_matches[cc].append(f"{gene.name} ({og.name})")
            
            # Add gene accession numbers
            for accession in m['Accession'].item().split(','):
                gene.accession.append(accession)
                og.accession.append(accession)
        else:
            failed_matches.append(gene.name)
    
    # Output logging messages for any failed matches
    if len(failed_matches) > 0:
        logger.error(f"{db_name}: Failed to find the genes in the annotation file ({mappingfile}): {', '.join(failed_matches)}")

    if len(failed_type_matches) > 0: 
        failed_type_matches_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k,v in failed_type_matches.items()])
        logger.warning(f"{db_name}: Did not figure out what type of annotations for:\n" + failed_type_matches_matches_string)

    if len(failed_class_matches) > 0: 
        ffiled_class_matches_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k,v in failed_class_matches.items()])
        logger.warning(f"{db_name}: Did not figure out class annotations for:\n" + ffiled_class_matches_matches_string)
    
    # Log the successful addition of annotations
    logger.success(f"Added {db_name} annotations to the PanRes ontology.")
