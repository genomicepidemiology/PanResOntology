from owlready2 import *
import re
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

def add_argannot_annotations(onto: Ontology, logger, db_name: str = 'ARGANNOT'):
    """Add ARGANNOT annotations to the ontology.

    Parameters
    ----------
    onto : Ontology
        The ontology object to load annotations into
    logger : _type_
        Logger object for logging messages
    db_name : str, optional
        Name of the database, by default 'ARGANNOT'
    """

    # Find genes from the specified database in the ontology
    matched_genes = find_genes_from_database(onto, database_name=db_name)

    # Lists for storing failed matches
    failed_regex_matches = list()
    failed_ab_matches = defaultdict(list)

    # Compile a regular expression to match class acronyms in fasta headers
    p = re.compile(r"\|\((\w{3,4})\)")

    # Loop through matched pan_ genes and original gene names to add the annotations
    for gene, og in matched_genes.items():
        # Extract the fasta header for the current gene
        fasta_header = [fh for fh in og.original_fasta_header if db_name in fh][0]
        
        # Match the class acronym using the regular expression
        m = p.findall(fasta_header)
        ab_acronym = m[0] if m else None

        # Check if the acronym is valid and is in the acr2class dictionary
        if ab_acronym is not None and ab_acronym in acr2class.keys():
            antibiotics = acr2class[ab_acronym].title()

            # Ad the annotations for the gene in question
            for ab in antibiotics.split('/'):
                success_match = gene_target(gene, og, target=ab, onto=onto, db_name=db_name)

                # log if no successful match
                if not success_match:
                    failed_ab_matches[ab].append(f"{gene.name} ({og.name})")
        
        # log if no successful acronym match
        else:
            failed_regex_matches.append(f"{gene.name} ({og.name})")
        
    # Output logging messages if any failed matches
    if len(failed_regex_matches) > 0:
        logger.warning(f"{db_name}: Failed to extract class acronyms: {', '.join(failed_regex_matches)}")    

    if len(failed_ab_matches) > 0:
        failed_ab_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k, v in failed_ab_matches.items()])
        logger.warning("{db_name}: Failed to find target annotations for:\n" + failed_ab_matches_string)
    
    # Log the successful addition of annotations
    logger.success(f"Added {db_name} annotations to the PanRes ontology.")
