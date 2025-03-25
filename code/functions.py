import subprocess
from owlready2 import *
import pandas as pd

def get_instance(onto: Ontology, name: str) -> Thing:
    """Retrieves an instance from the ontology based on its name

    Parameters
    ----------
    onto : Ontology
        The ontology object
    name : str
        The name of the instance to retrieve

    Returns
    -------
    Thing
        The instance if found, otherwise None
    """

    # Clean the name by replacing spaces and hyphens with underscores
    cleaned_name = name.replace(" ", "_").replace("-", "_")

    # Construct the full IRI for the instance
    full_iri = onto.base_iri + cleaned_name

    # Search for the instance
    instance = onto.search_one(iri = full_iri)
    return instance


def get_or_create_instance(onto: Ontology, cls: Thing, name: str) -> Thing:
    """Retrieves or creates an instance in the ontology

    Parameters
    ----------
    onto : Ontology
        The ontology object
    cls : Thing
        The class of hte instance to create
    name : str
        The name of the instance

    Returns
    -------
    Thing
        The instance
    """ 
    
    # Construct the full IRI for the instance
    full_iri = onto.base_iri + name
    
    # Search for the instance
    instance = onto.search_one(iri=full_iri)#"*{}".format(name))
    if instance is None:
        # Create the instance if it doesn't exist
        instance = cls(name)

    return instance

def get_or_create_subclass(onto: Ontology, parent_cls: Thing, subclass_name: str) -> Thing:
    """Retrieves or creates a subclass in the ontology

    Parameters
    ----------
    onto : Ontology
        The ontology object
    parent_cls : Thing
        The parent class of the subclass.
    subclass_name : str
        The name of the subclass.

    Returns
    -------
    Thing
        The subclass
    """
    
    # Clean the subclass name by replacing spaces and hyphens with underscores
    cleaned_name = subclass_name.replace(" ", "_").replace("-", "_")

    # Search for the subclass in the ontology
    subclass_instance = onto.search_one(iri = onto.base_iri + cleaned_name)
    
    # Create if it doesnt exist
    if subclass_instance is None:
        subclass_instance = types.new_class(cleaned_name, (parent_cls, ))

        # Set the label to the original name with spaces
        subclass_instance.label = [cleaned_name]
    
    return subclass_instance


def find_original_name(gene_instance: Thing, database_name: str) -> Thing:
    """
    Finds the original name of a gene instance from a specific database.

    Parameters
    ----------
    gene_instance : Thing
        The gene instance.
    database_name : str
        The name of the database.

    Returns
    -------
    Thing
        The original gene instance if found, otherwise None.
    """
    # Check if the gene has an original name
    for og in gene_instance.same_as:
        for ogname in og.is_from_database:
            if ogname.name == database_name:
                return og

def find_genes_from_database(onto: Ontology, database_name: str) -> dict:    
    """
    Finds genes from a specific database in the ontology

    Parameters
    ----------
    onto : Ontology
        The ontology object
    database_name : str
        The name of the database

    Returns
    -------
    dict
        A dictionary mapping genes to their original gene instances.
    """

    # Search for the database instance
    database_instance = onto.search_one(iri="*{}".format(database_name))
    
    if not database_instance:
        print(f"Database '{database_name}' not found in the ontology.")
        return []
    
    # Find all genes associated with the database
    genes = [gene for gene in onto.PanGene.instances() if database_instance in gene.is_from_database and onto.PanGene in gene.is_a]

    # Map genes to their original gene instances
    gene2og = {gene: find_original_name(gene, database_name) for gene in genes}
    return gene2og

def clean_gene_name(gene_name: str, db: str) -> str:
    """Cleans the gene name based on the database.

    Parameters
    ----------
    gene_name : str
        The gene name to clean.
    db : str
        The name of the database the gene name is from.

    Returns
    -------
    str
        The cleaned gene name.
    """
    
    # Remove database prefix from the gene name
    gene_name = gene_name.replace(db + '|','')
    db = db.lower()
    
    # Clean the gene name based on the database
    if db == 'amrfinderplus': 
        return gene_name.split('|')[5]
    elif db == 'card_amr':
        return gene_name.split('|')[5].split(' [')[0]
    elif db == 'megares':
        return gene_name.split('|')[4]
    elif db == 'argannot':
        return ")".join(gene_name.split('|')[0].split(')')[1:])
    elif db == 'functional_amr': 
        return gene_name.split('|')[1]
    elif db == 'metalres': 
        return gene_name.split(' ')[0]
    else:
        return gene_name    

def accession_to_pubmed(accession: str) -> list:
    """
    Retrieves PubMed IDs associated with a protein or gene accession

    Parameters
    ----------
    accession : str
        The protein accession

    Returns
    -------
    list
        A list of PubMed IDs
    """

    # Construct the command to retrieve PubMed IDs
    cmd = f"esearch -db protein -query {accession} | elink -target pubmed | efetch -format uid"

    # Run the command
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # Check if the command was successful
    if p.returncode == 0:
        return p.stdout.decode().strip().split()
    
    return None