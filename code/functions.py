import subprocess
from numpy import full
from owlready2 import *
import pandas as pd

def get_instance(onto, name):
    # instance = onto.search_one(iri="*{}".format(name))
    cleaned_name = name.replace(" ", "_").replace("-", "_")
    full_iri = onto.base_iri + cleaned_name
    instance = onto.search_one(iri = full_iri)
    return instance

# Function to get or create an instance
def get_or_create_instance(onto: Ontology, cls: Thing, name: str) -> Thing:
    
    # Construct the full IRI for the instance
    full_iri = onto.base_iri + name
    
    # Search for the instance
    instance = onto.search_one(iri=full_iri)#"*{}".format(name))
    if instance is None:
        # Create the instance if it doesn't exist
        instance = cls(name)

    return instance

def get_or_create_subclass(onto: Ontology, parent_cls: Thing, subclass_name: str) -> Thing:
    
    # Remove spaces
    cleaned_name = subclass_name.replace(" ", "_").replace("-", "_")

    subclass_instance = onto.search_one(iri = onto.base_iri + cleaned_name)
    
    # Create if it doesnt exist
    if subclass_instance is None:
        subclass_instance = types.new_class(cleaned_name, (parent_cls, ))

        # Set the label to the original name with spaces
        subclass_instance.label = [cleaned_name]
    
    return subclass_instance


def find_original_name(gene_instance, database_name):
    # Check if the gene has an original name
    for og in gene_instance.equivalent_to:
        for ogname in og.is_from_database:
            if ogname.name == database_name:
                return og

def find_genes_from_database(onto, database_name):    
    # Search for the database instance
    database_instance = onto.search_one(iri="*{}".format(database_name))
    
    if not database_instance:
        print(f"Database '{database_name}' not found in the ontology.")
        return []
    
    # Find all genes associated with the database
    genes = [gene for gene in onto.PanGene.instances() if database_instance in gene.is_from_database and onto.PanGene in gene.is_a]

    gene2og = {gene: find_original_name(gene, database_name) for gene in genes}
    return gene2og

def clean_gene_name(gene_name, db):
    gene_name = gene_name.replace(db + '|','')
    db = db.lower()
    if db == 'amrfinderplus': #in ['amrfinderplus', 'card_amr']:
        return gene_name.split('|')[5]
    elif db == 'card_amr':
        return gene_name.split('|')[5].split(' [')[0]
    elif db == 'megares':
        return gene_name.split('|')[0]
    elif db == 'argannot':
        return ")".join(gene_name.split('|')[0].split(')')[1:])
    elif db == 'functional_amr': 
        return gene_name.split('|')[1]
    elif db == 'metalres': 
        return gene_name.split(' ')[0]
    else:
        return gene_name    
    
def check_targets(excelfile):
    antibiotics = pd.read_excel(excelfile, sheet_name='antibiotic')
    antibiotics['drug'] = antibiotics['drug'].str.strip().str.title()
    to_drop = [c for c in antibiotics.columns if c.startswith('Unnamed')]
    antibiotics['group'] = antibiotics['group'].str.replace(r's$', '', regex=True)
    antibiotics['class'] = antibiotics['class'].str.strip().str.title()

    metals = pd.read_excel(excelfile, sheet_name='metals')['Metal'].str.lower().str.title()
    biocides = pd.read_excel(excelfile, sheet_name='biocides')['Biocide'].str.lower().str.title()

    return {
        'antibiotics': antibiotics.drop(columns=to_drop), 
        'metals': metals.values.tolist(),
        'biocides': biocides.values.tolist()
    }

def find_target_annotation(target,annotated_targets):

    antibiotics = annotated_targets['antibiotics']
    drug_match = antibiotics['drug'] == target
    if drug_match.any():
        return 'antibiotic', antibiotics.loc[drug_match, ]
    
    class_match = antibiotics['class'] == target
    if class_match.any():
        return 'antibiotic_class', antibiotics.loc[class_match]

    subclass_match = antibiotics['group'] == target
    if subclass_match.any():
        return 'antibiotic_subclass', antibiotics.loc[subclass_match]
    

    metals = annotated_targets['metals']
    if target in metals:
        return ('metal', target)

    biocides = annotated_targets['biocides']
    if target in biocides:
        return ('biocide', target)
    

def accession_to_pubmed(accession: str):

    cmd = f"esearch -db protein -query {accession} | elink -target pubmed | efetch -format uid"
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if p.returncode == 0:
        return p.stdout.decode().strip().split()
    
    return None