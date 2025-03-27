import subprocess
from owlready2 import *
import pandas as pd
from graphviz import Digraph
from IPython.display import Image, display

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

def get_subclasses(onto: Thing, class_name: str) -> pd.DataFrame: 
    """Find subclasses of an Ontology class

    Parameters
    ----------
    onto : Thing
        The loaded ontology to query
    class_name : str
        Subclasses of class_name to find

    Returns
    -------
    pd.DataFrame
        Column containing subclass matches of the class_name
    """
    class_match = onto.search_one(iri=f"*{class_name}")
    subclasses = [sc.name for sc in list(class_match.subclasses())]
    
    return pd.DataFrame([class_match.name] + subclasses, columns=['match'])

def class_to_genes(onto: Thing, class_instance: Thing, annotation_property: Thing) -> list:
    """Find the Pan genes that are linked by the annotation_property to the class instance

    Parameters
    ----------
    onto : Thing
        The loaded ontology to query
    class_instance : Thing
        The class instance to match for
    annotation_property : Thing
        The annotaiton property that should link the pangene to the class instance

    Returns
    -------
    list
        list of pan genes matched
    """
    
    pan_instances = pd.DataFrame(onto.PanGene.instances(), columns=['instance'])
    pan_instances['has_resistance_to_class'] = pan_instances['instance'].apply(
        lambda x: class_instance in annotation_property[x]
    )
    pan_instances = pan_instances.loc[
        pan_instances.has_resistance_to_class == True,
    ].drop(columns='has_resistance_to_class').values.tolist()  
    return pan_instances

def summarise_classes(onto: Ontology, class_name: str) -> pd.DataFrame:
    """Summarise children of an ontology class and the genes associated with the class.

    Parameters
    ----------
    onto : Ontology
        The loaded ontology to query.
    class_name : str
        The name of the ontology class in question

    Returns
    -------
    pd.DataFrame
        Return a pandas DataFrame that shows the class names, 
        class types of the match and the children, 
        and the number of genes.
    """
    subclasses = get_subclasses(onto, class_name)
    subclasses['instance'] = [onto.search_one(iri=f"*{c}") for c in subclasses.iloc[:, 0]]
    subclasses['type'] = subclasses['instance'].apply(lambda x: x.is_a[0].name)
    subclasses['n_genes'] = 0
    
    if any(subclasses['type'].str.endswith('Class')):
        annotation_property = onto.has_resistance_class
        
        subclasses.loc[
            subclasses['type'].str.endswith('Class'), 'n_genes'
        ] = subclasses.loc[
            subclasses['type'].str.endswith('Class'), 'instance'
        ].apply(
            lambda x: class_to_genes(onto = onto, class_instance = x, annotation_property = annotation_property)
        ).apply(len)
        
        
    if any(subclasses['type'].str.endswith('Phenotype')):
        annotation_property = onto.has_predicted_phenotype
    
        subclasses.loc[
            subclasses['type'].str.endswith('Phenotype'), 'n_genes'
        ] = subclasses.loc[
            subclasses['type'].str.endswith('Phenotype'), 'instance'
        ].apply(
            lambda x: class_to_genes(onto = onto, class_instance = x, annotation_property = annotation_property)
        ).apply(len)
        
    
    return subclasses.drop(columns=['instance']).rename(columns = {'match': 'Class'})

def get_annotations_of_individual(individual: Thing):
    """Get all annotations of an individual in the the ontology"""
    return individual.get_annotations()

def visualize_specific_classes(onto: Ontology, class_names: list, output_file: str ="specific_classes_visualization"):
    """
    Visualizes the relationship between specific classes in the given ontology using Graphviz.

    Parameters
    ----------
    onto : Ontology
        The ontology to visualize
    class_names : list
        A list of class names to visualize.
    output_file : str, optional
        Name of file to save the graph to, by default "specific_classes_visualization"

    
    Examples
    ----------
    >>> specific_classees = ["PanGene, "OriginalGene", "Database"]
    >>> visualize_specific_classes(onto = onto, specific_classes)
    """
    dot = Digraph(comment='Specific Classes Visualization')
    
    def add_class_and_relationships(cls):
        dot.node(cls.name, cls.name)
        for parent in cls.is_a:
            if isinstance(parent, owlready2.ThingClass):
                dot.edge(parent.name, cls.name)
                dot.node(parent.name, parent.name)
    
    # Add nodes and edges for specified classes
    for class_name in class_names:
        cls = onto[class_name]
        add_class_and_relationships(cls)
    
    # Save and render the graph
    dot.format = 'png'
    dot.render(output_file, format='png', cleanup=True)
    print(f"Visualization saved as {output_file}.png")
    
    display(Image(filename=output_file + '.png'))