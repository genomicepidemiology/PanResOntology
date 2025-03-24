from tkinter import N
import pandas as pd
from owlready2 import *
import subprocess

import sys
sys.path.append('..')
from functions import get_instance, clean_gene_name, get_or_create_instance, get_or_create_subclass
import re

database2name = {
    'amrfinderplus': 'AMRFinderPlus',
    'argannot': 'ARGANNOT',
    'card_amr' :'CARD', 
    'csabapal': 'CsabaPal',
    'functional_amr': 'ResFinderFG',
    'megares': 'MegaRes',
    'metalres': 'MetalRes',
    'resfinder': 'ResFinder'
}

def add_panres_genes(file: str, onto: Ontology, logger):
    '''
    Loop through panres genes and add them
    '''

    # Load file
    panres_metadata = pd.read_csv(file, sep=',' if file.endswith('.csv') else '\t', skiprows=1)
    panres_metadata['userGeneName'] = panres_metadata['userGeneName'].str.replace("_v1.0.0", "", regex=False)
    panres_metadata['chosenSeq'] = panres_metadata['chosenSeq'].str.replace("_v1.0.0", "", regex=False)
    panres_metadata['database'] = panres_metadata['database'].str.replace("_genes", "")

    genes = panres_metadata['userGeneName'].unique().tolist()


    for gene in set(genes):
        new_gene = onto.PanGene(gene)

        # Get annotation data
        m = panres_metadata.loc[panres_metadata['userGeneName'] == gene, ]

        # Gene length
        if m['gene_len'].nunique() == 1:
            new_gene.has_length.append(int(m['gene_len'].values[0]))
        else:
            logger.warning(f"PanRes: {gene} has multiple gene lengths associated.")
        
        # Loop through the annotation data for this gene
        for _, row in m.iterrows():
            database_shortname = row['database']
            database_name = database2name[database_shortname.lower()]
            database_instance = get_instance(onto, database_name)
            new_gene.is_from_database.append(database_instance)

            fasta_header = row['fa_header'].replace("~~~", "|").replace("'", "")
            gene_name = clean_gene_name(fasta_header, database_shortname.lower()).strip()
            original_gene_instance = onto.OriginalGene(gene_name)
                        
            original_gene_instance.is_from_database.append(database_instance)
            original_gene_instance.original_fasta_header.append(fasta_header.strip() + '|' + database_name)

            new_gene.same_as.append(original_gene_instance)
            # original_gene_instance.subclass_of.append(new_gene)

            # gene cluster
            new_cluster = get_or_create_instance(onto = onto, cls = onto.PanGeneCluster, name=row['chosenSeq'].replace('pan', 'panc'))
            new_gene.member_of.append(new_cluster)
    
    # logger.info(f"Adding {len(genes)} PanRes genes to the ontology.")
    logger.success(f"Added PanRes genes (n={len(genes)}) to the ontology.")

def add_panres_proteins(file: str, clstrs: str, onto: Ontology, logger):

    p = subprocess.run(f"grep '>' {file}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    protein_names = [l.split('_v1')[0].replace(">", "") for l in p.stdout.decode().split()]

    gene_cluster_instances = onto.PanGeneCluster.instances()
    for protein_name in protein_names:
        protein_instance = get_or_create_instance(onto = onto, cls = onto.PanProtein, name=protein_name.title())
        gene_cluster = get_instance(onto = onto, name = protein_name)
        if gene_cluster is not None:
            gene_cluster.translates_to.append(protein_instance)    
    
        
    
    logger.success("Added PanRes proteins to the ontology.")

    ## Read the CD-HIT output file
    p = re.compile(r">(pan\_\d+)")
    clusters = defaultdict(list)
    cluster_representative = None
    members = []
    with open(clstrs, 'r') as f:
        for l in f.readlines():
            l = l.strip()
            m = p.findall(l)
            if l.startswith('>Cluster'):
                if cluster_representative is not None:
                    clusters[cluster_representative] = members
                    members, cluster_representative = [], None
            elif m and l.endswith('*'):
                cluster_representative = m[0]
                members.append(m[0])
            else:
                members.append(m[0])

    for cl, cl_members in clusters.items():
        cluster_protein_instance = get_or_create_instance(onto = onto, cls = onto.PanProteinCluster, name=cl.replace('pan', 'Pan'))
        for cl_member in cl_members:
            cl_instance = get_or_create_instance(onto = onto, cls = onto.PanProtein, name = cl_member.title())
            if cl_instance is not None:
                cl_instance.member_of.append(cluster_protein_instance)
        
    
    logger.success("Added PanRes protein clusters to the ontology.")
