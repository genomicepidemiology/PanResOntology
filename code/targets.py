import pandas as pd
from owlready2 import Thing, Ontology, destroy_entity
from functions import get_or_create_subclass, get_instance

def load_targets(excelfile: str, onto: Ontology, logger=None):

    # Antibiotics
    antibiotics = pd.read_excel(excelfile, sheet_name='antibiotic')
    antibiotics['drug'] = antibiotics['drug'].str.strip().str.title()
    antibiotics['group'] = antibiotics['group'].str.replace(r's$', '', regex=True)
    antibiotics['class'] = antibiotics['class'].str.strip().str.title()

    for _, row in antibiotics.iterrows():
        
        # Add or get class
        ab_class_instance = get_or_create_subclass(
            onto = onto,
            parent_cls=onto.AntibioticResistanceClass,
            subclass_name=row['class']
        )

        # Add or get phenotype
        ab_phenotype_instance = get_or_create_subclass(
            onto = onto,
            parent_cls = onto.AntibioticResistancePhenotype,
            subclass_name=row['drug']
        )

        try:# ab_class_instance not in ab_phenotype_instance.is_a:
            ab_phenotype_instance.is_a.append(ab_class_instance) 
        except TypeError:
            pass
    
    # Metals
    metals = pd.read_excel(excelfile, sheet_name='metals')#['Metal'].str.lower().str.title()
    metals[['symbol', 'note']] = metals[['symbol', 'note']].fillna('')
    for _ , row in metals.iterrows():
        metal_instance = get_or_create_subclass(
            onto = onto,
            parent_cls = onto.Metal,
            subclass_name = row['Metal'].strip().lower().title()
        )

        if len(row['symbol']) > 0:
            metal_instance.metal_symbol.append(row['symbol'])
        if len(row['note']) > 0:
            metal_instance.metal_commenent.append(f"{row['Metal']} {row['note']}") 

    # Metal
    biocides = pd.read_excel(excelfile, sheet_name='biocides')
    biocides['Biocide'] = biocides['Biocide'].str.strip().str.lower().str.title()
    biocides['Class'] = biocides['Class'].str.strip().str.lower().str.title()
    for _, row in biocides.iterrows():
        biocide_instance = get_or_create_subclass(
            onto = onto,
            parent_cls = onto.Biocide,
            subclass_name = row['Biocide']
        )

        biocide_class_instance = get_or_create_subclass(
            onto = onto,
            parent_cls = onto.BiocideClass,
            subclass_name = row['Class']
        )

        biocide_instance.is_a.append(biocide_class_instance)

    # Unclassified
    unclassified = pd.read_excel(excelfile, sheet_name='unclassified')
    for _, row in unclassified.iterrows():
        unclassified_instance = get_or_create_subclass(
            onto = onto,
            parent_cls = onto.UnclassifiedResistance,
            subclass_name = row['Compound'].strip().lower().title(),
        )

        try:
            unclassified_class_instance = get_or_create_subclass(
                onto = onto,
                parent_cls = onto.UnclassifiedResistanceClass,
                subclass_name = row['Class'].strip().lower().title(),
            )

            unclassified_instance.is_a.append(unclassified_class_instance)
        except: 
            pass

    if logger is not None:
        logger.success("Loaded drugs, biocide and metal targets into the ontology.")
    

def remove_unused_subclasses_with_property(onto: Ontology, parent_cls: Thing, property_name: str, logger):
    """
    Remove all unused subclasses of a given parent class in the ontology that do not have a specific object property.

    Parameters:
    onto (Ontology): The ontology object.
    parent_cls (Thing): The parent class whose unused subclasses will be removed.
    property_name (str): The name of the object property to check.
    """
    unused_subclasses = []

    if parent_cls is None and logger :
        logger.error(f"Parent class {parent_cls} is not defined in the ontology.")
        return
    
    for subclass in parent_cls.subclasses():
        has_property = False
        for instance in onto.PanGene.instances():
            # Check if the subclass has the specified object property
            if subclass in getattr(instance, property_name, []):
                has_property = True
                break
        
        if not has_property:
            unused_subclasses.append(subclass)
    
    # Remove unused subclasses
    for subclass in unused_subclasses:
        destroy_entity(subclass)
    
    logger.info(f"Destroyed {len(unused_subclasses)} subclasses of {parent_cls.name}")


def gene_target(gene: Thing, og: Thing, target: str, onto: Ontology, db_name: str = None):

    target = target.replace(" ", "_").replace("-", "_").title()

    # get class from the ontology
    target_instance = None
    if '+' in target:
        target_instance = get_or_create_subclass(onto = onto, parent_cls=onto.AntibioticResistancePhenotype, subclass_name=target)
        target_instance.is_drug_combination.append(True)
        # split +
        for ts in target.split('+'):
            ts_instance = get_instance(onto = onto, name = ts.replace(" ", "_"))
            if ts_instance is not None:
                target_instance.is_a.append(ts_instance)
    else:
        target_instance = get_instance(onto = onto, name = target)
    
    if target_instance is None:
        return False
    
    # Check if the target is a phenotype
    if onto.AntibioticResistancePhenotype in target_instance.is_a or onto.Metal in target_instance.is_a or onto.Biocide in target_instance.is_a or onto.UnclassifiedResistance in target_instance.is_a:
        gene.has_predicted_phenotype.append(target_instance)
        og.has_predicted_phenotype.append(target_instance)
    elif onto.AntibioticResistanceClass in target_instance.is_a or onto.BiocideClass in target_instance.is_a or onto.UnclassifiedResistanceClass in target_instance.is_a:
        gene.has_resistance_class.append(target_instance)
        og.has_resistance_class.append(target_instance)
    else:
        return False

    try:        
        db_instance = get_instance(onto = onto, name = db_name)
        target_instance.found_in.append(db_instance)
    except TypeError:
        pass
    
    return True