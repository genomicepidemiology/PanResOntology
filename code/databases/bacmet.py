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

    # Load file
    exp_annotations = pd.read_csv(mappingfile, sep='\t')
    exp_annotations['gene_name'] = exp_annotations['Gene_name'].str.lower().str.split('/')
    exp_annotations = exp_annotations.explode('gene_name')


    failed_matches = []
    failed_type_matches = defaultdict(list)
    failed_class_matches = defaultdict(list)

    p = re.compile(r"^((\w+)(\s\w+)?)\s")
    p_class = re.compile(r"^((\w+)(\s\w+)?)\s(\[class\W+((\w+)(\s\w+)?)\])")

    db_instance = get_instance(onto=onto, name=db_name)
    genes = [gene for gene in onto.PanGene.instances() if db_instance in gene.is_from_database]
    matched_genes = find_genes_from_database(onto = onto, database_name = db_name)
    for gene, meg_gene in matched_genes.items():
        meg_header = meg_gene.original_fasta_header[0]
        meg_gene_name = meg_header.split('|')[-2].lower()
        m = exp_annotations.loc[
            (exp_annotations['gene_name'] == meg_gene_name)
        ].groupby(['Gene_name']).agg(agg_funcs).reset_index()

        if m.shape[0] == 1:
            compounds = m['Compound'].item()
            original_name = m['Gene_name'].item()
            og = get_or_create_instance(onto = onto, cls=onto.OriginalGene, name=original_name)
            og.is_from_database.append(db_instance)
            gene.equivalent_to.append(og)

            # set metal resistance
            for compound in set(compounds.split(',')):
                if 'class' in compound:
                    m_compound = p_class.findall(compound.strip())
                else:
                    m_compound = p.findall(compound.strip())
                if len(m_compound) > 0:
                    c = m_compound[0][0]
                    if not gene_target(gene=gene, og=og, target=c, onto=onto, db_name=db_name):
                        failed_type_matches[c].append(f"{gene.name} ({og.name})")
                    if len(m_compound[0]) > 3:
                        cc = m_compound[0][4]
                        if not gene_target(gene=gene, og=og, target=cc, onto=onto, db_name=db_name):
                            failed_class_matches[cc].append(f"{gene.name} ({og.name})")
            
            for accession in m['Accession'].item().split(','):
                gene.accession.append(accession)
                og.accession.append(accession)

                # pubmed = accession_to_pubmed(accession=accession)
                # if pubmed is not None:
                #     for pb in pubmed:
                #         gene.pubmed.append(pb)
                #         og.pubmed.append(pb)
        else:
            failed_matches.append(gene.name)

    if len(failed_matches) > 0:
        logger.error(f"{db_name}: Failed to find the genes in the annotation file ({mappingfile}): {', '.join(failed_matches)}")    
    

    if len(failed_type_matches) > 0: 
        failed_type_matches_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k,v in failed_type_matches.items()])
        logger.warning(f"{db_name}: Did not figure out what type of annotations for:\n" + failed_type_matches_matches_string)

    if len(failed_class_matches) > 0: 
        ffiled_class_matches_matches_string = "\n".join([f"{k}: {', '.join(v)}" for k,v in failed_class_matches.items()])
        logger.warning(f"{db_name}: Did not figure out class annotations for:\n" + ffiled_class_matches_matches_string)

    logger.success(f"Added {db_name} annotations to the PanRes ontology.")
