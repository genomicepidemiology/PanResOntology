import argparse
from owlready2 import get_ontology, Thing, ThingClass
import pandas as pd
import os

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-f', '--file', 
        type=str, 
        required=True, 
        help='Path to the ontology file',
        dest='file'
    )
    parser.add_argument(
        '-o', '--output', 
        type=str, 
        required=True, 
        help='Path to the output file',
        dest='output'
    )

    parser.add_argument(
        '-c', '--columns',
        type=str,
        nargs='+',
        default=['name', 'accession', 'has_predicted_phenotype', 'has_resistance_class', 'is_from_database', 'same_as'],
        help='Features and classes to include in the output file',
        dest='columns',
        choices=[
            'name', 
            'accession', 
            'has_length',
            'pubmed', 
            'card_link', 
            'has_predicted_phenotype', 
            'has_resistance_class', 
            'is_from_database', 
            'same_as', 
            'gene_alt_name', 
            'member_of', 
            'has_members',
            'translates_to',
            'has_mechanism_of_resistance',
            'original_fasta_header',
            'folds_to'
        ]

    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    
    # Load the ontology
    onto = get_ontology(args.file).load()  

    # Get the data specified by the user
    data = []
    for gene in onto.PanGene.instances():
        row = {}
        for col in args.columns:
            value = getattr(gene, col)
            if isinstance(value, list):
                vv = []
                for v in value:
                    if isinstance(v, ThingClass) or isinstance(v, Thing):
                        vv.append(str(v.name))
                    else:
                        vv.append(str(v))
                value = ';'.join(vv)
            elif isinstance(value, Thing):
                value = str(value.name)
            row[col] = value
        data.append(row)
    
    # Convert data to pandas dataframe
    df = pd.DataFrame(data)

    # Replace empty strings with NaN
    df = df.replace('', pd.NA)
    
    # Drop empty rows in attribute columns
    columns = [col for col in args.columns if col not in ['name']]
    df = df.dropna(subset=columns, how='all')

    # Save the dataframe to a CSV file
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    df.to_csv(args.output, index=False)