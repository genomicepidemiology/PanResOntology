from turtle import left
import pandas as pd
from owlready2 import *
import subprocess
from loguru import logger

import sys
sys.path.append('..')
from functions import get_instance, clean_gene_name, get_or_create_instance
import re
import os

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

def discarded_genes(discarded_list: str):
    """Read a list of discarded genes from a file and return them as a list.

    Parameters
    ----------
    discarded_list : str
        Path to the file containing the list of discarded genes.

    Returns
    -------
    list
        A list of discarded gene names.
    """
    
    with open(discarded_list, 'r') as f:
        discarded = [l.strip().split('_v')[0].replace('>', '') for l in f.readlines()]
    
    return discarded

def add_panres_genes(file: str, onto: Ontology, discarded: str, logger = logger):
    """Add PanRes genes to the ontology

    Parameters
    ----------
    file : str
        File containing metadata for the genes included in the first version of PanRes
    onto : Ontology
        The Ontology to load the genes into
    discarded : str
        File with gene names that have been discarded from the PanRes database, e.g. because their sequence didnt contain correct start and stop codons
    logger : loguru.logger
        Logger object for logging messages.
    """

    # Load the metadata file (.tsv or .csv)

    # Clean up gene names and database names in the metadata
    panres_metadata = pd.read_csv(file, sep=',' if file.endswith('.csv') else '\t', skiprows=1)
    panres_metadata['userGeneName'] = panres_metadata['userGeneName'].str.replace("_v1.0.0", "", regex=False)
    panres_metadata['chosenSeq'] = panres_metadata['chosenSeq'].str.replace("_v1.0.0", "", regex=False)
    panres_metadata['database'] = panres_metadata['database'].str.replace("_genes", "")


    # Get list of discarded genes
    discarded_genes_list = discarded_genes(discarded)
    discarded_df = pd.DataFrame({'pan_gene': discarded_genes_list, 'status': 'discarded', 'reason': ''})
    discarded_df = discarded_df.merge(
        panres_metadata.pivot_table(
            index='userGeneName',
            columns = 'database',
            values='fa_header',
            aggfunc=','.join
        ),
        left_on = 'pan_gene',
        right_index=True,
    )
    discarded_df.columns = discarded_df.columns.str.lower()#
    discarded_df = discarded_df.rename(columns = database2name)
    discarded_df.to_csv(os.path.join(os.path.dirname(discarded), 'panres_discarded_genes.csv'), index=False)

    # Remove discarded genes from the metadata
    panres_metadata = panres_metadata[~panres_metadata['userGeneName'].isin(discarded_genes_list)]
    
    # Get a list of the unique gene names
    genes = panres_metadata['userGeneName'].unique().tolist()
    # genes = panres_metadata.loc[panres_metadata['database'] == 'resfinder', 'userGeneName'].unique().tolist()
    
    # Loop through each unique gene, add it to the ontology and its annotations
    for gene in set(genes):

        # Get annotation data for the correct gene
        m = panres_metadata.loc[panres_metadata['userGeneName'] == gene, ]
        
        # Make new PanGene instance
        new_gene = onto.PanGene(gene)

        # add the gene length if its consistent across all entries
        if m['gene_len'].nunique() == 1:
            new_gene.has_length.append(int(m['gene_len'].values[0]))
        else:
            logger.warning(f"PanRes: {gene} has multiple gene lengths associated.")
        
        # Loop through the annotation data for this gene
        for _, row in m.iterrows():
            # Map the full database name from the short name
            database_shortname = row['database']
            database_name = database2name[database_shortname.lower()]
            
            # Get the instance and link the original gene name to the database
            database_instance = get_instance(onto, database_name)
            new_gene.is_from_database.append(database_instance)

            # Clean and format the fasta header
            fasta_header = row['fa_header'].replace("~~~", "|").replace("'", "")
            gene_name = clean_gene_name(fasta_header, database_shortname.lower()).strip()

            # Add the cleaned gene name as an individual of the class OriginalGene
            original_gene_instance = onto.OriginalGene(gene_name)
                        
            # Track where the original gene name is from
            original_gene_instance.is_from_database.append(database_instance)
            original_gene_instance.original_fasta_header.append(fasta_header.strip() + '|' + database_name)

            # Annotate that the pan gene name is the same as the original gene name
            new_gene.same_as.append(original_gene_instance)

            # Get which cluster the gene belongs to 
            new_cluster = get_or_create_instance(onto = onto, cls = onto.PanGeneCluster, name=row['chosenSeq'].replace('pan', 'panc'))
            new_gene.member_of.append(new_cluster)
    
    # logger.info(f"Adding {len(genes)} PanRes genes to the ontology.")
    logger.success(f"Added PanRes genes (n={len(genes)}) to the ontology.")

def add_panres_proteins(file: str, clstrs: str, onto: Ontology, logger):
    """Adds PanRes proteins and their clusters to the ontology.

    Parameters
    ----------
    file : str
        Path to the faa file containing PanRes protein sequences.
    clstrs : str
        Path to the clstr file containing PanRes protein clusters, as determined by CD-HIT.
    onto : Ontology
        The ontology to load the protein information into
    logger : loguru.logger
        Logger object for logging messages
    """

    # Grep all headers of the protein sequences and extract them
    p = subprocess.run(f"grep '>' {file}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    protein_names = [l.split('_v1')[0].replace(">", "") for l in p.stdout.decode().split()]

    # Loop through the protein names and add to the ontology
    for protein_name in protein_names:
        protein_instance = get_or_create_instance(onto = onto, cls = onto.PanProtein, name=protein_name.title())

        # Get the corresponding gene cluster instance and link the protein name to it
        gene_cluster = get_instance(onto = onto, name = protein_name)
        if gene_cluster is not None:
            gene_cluster.translates_to.append(protein_instance)    
    
    # Log the successfull addition  of proteins
    logger.success("Added PanRes proteins to the ontology.")

    # Read the CD-HIT output file to extract cluster information
    p = re.compile(r">(pan\_\d+)")
    clusters = defaultdict(list)
    cluster_representative = None
    members = []
    with open(clstrs, 'r') as f:
        for l in f.readlines():
            l = l.strip()
            m = p.findall(l)

            # Check if the line indicates a new cluster
            if l.startswith('>Cluster'):
                if cluster_representative is not None:
                    clusters[cluster_representative] = members
                    members, cluster_representative = [], None
            # Else identify the cluster representative and it to the member list
            elif m and l.endswith('*'):
                cluster_representative = m[0]
                members.append(m[0])
            # Add members to the current cluster
            else:
                members.append(m[0])

    # Loop through each cluster and add it to the ontology
    for cl, cl_members in clusters.items():
        cluster_protein_instance = get_or_create_instance(onto = onto, cls = onto.PanProteinCluster, name=cl.replace('pan', 'Pan'))

        # Link each protein member to its cluster representative
        for cl_member in cl_members:
            cl_instance = get_or_create_instance(onto = onto, cls = onto.PanProtein, name = cl_member.title())
            if cl_instance is not None:
                cl_instance.member_of.append(cluster_protein_instance)
        
    # Log the successful addition of protein clusters
    logger.success("Added PanRes protein clusters to the ontology.")
