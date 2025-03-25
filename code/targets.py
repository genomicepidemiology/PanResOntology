import pandas as pd
from owlready2 import Thing, Ontology, destroy_entity
from functions import get_or_create_subclass, get_instance

def load_targets(excelfile: str, onto: Ontology, logger=None) -> None:
    """Load resistance targets from an Excel file into the ontology

    Parameters
    ----------
    excelfile : str
        Path to the Excel file  containing target data.
    onto : Ontology
        The ontology object to which the targets will be added.
    logger : loguru.logger, optional
        Logger object for logging messages, by default None
    """

    # Load antibiotic data from the Excel file
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

        if row['drug'] == row['class']:
            continue
        
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

    # Load metal data from the Excel file
    metals = pd.read_excel(excelfile, sheet_name='metals')
    metals[['symbol', 'note']] = metals[['symbol', 'note']].fillna('')
    for _ , row in metals.iterrows():
        # Add or get metal instance
        metal_instance = get_or_create_subclass(
            onto = onto,
            parent_cls = onto.Metal,
            subclass_name = row['Metal'].strip().lower().title()
        )

        # Add symbol and note if available
        if len(row['symbol']) > 0:
            metal_instance.metal_symbol.append(row['symbol'])
        if len(row['note']) > 0:
            metal_instance.metal_comment.append(f"{row['Metal']} {row['note']}") 

    # Load biocide data from the Excel file
    biocides = pd.read_excel(excelfile, sheet_name='biocides')
    biocides['Biocide'] = biocides['Biocide'].str.strip().str.lower().str.title()
    biocides['Class'] = biocides['Class'].str.strip().str.lower().str.title()
    for _, row in biocides.iterrows():
        # Add or get biocide class instance
        biocide_class_instance = get_or_create_subclass(
            onto = onto,
            parent_cls = onto.BiocideClass,
            subclass_name = row['Class']
        )
        
        if row['Class'] == row['Biocide']:
            continue
        
        # Add or get biocide instance
        biocide_instance = get_or_create_subclass(
            onto = onto,
            parent_cls = onto.Biocide,
            subclass_name = row['Biocide']
        )

        biocide_instance.is_a.append(biocide_class_instance)

    # # Load unclassified data from the Excel file
    unclassified = pd.read_excel(excelfile, sheet_name='unclassified')
    for _, row in unclassified.iterrows():
        # Add or get unclassified instance
        unclassified_instance = get_or_create_subclass(
            onto = onto,
            parent_cls = onto.UnclassifiedResistance,
            subclass_name = row['Compound'].strip().lower().title(),
        )

        try:
            # Add or get unclassified class instance
            unclassified_class_instance = get_or_create_subclass(
                onto = onto,
                parent_cls = onto.UnclassifiedResistanceClass,
                subclass_name = row['Class'].strip().lower().title(),
            )

            unclassified_instance.is_a.append(unclassified_class_instance)
        except: 
            pass
    
    # Log the successful loading of targets
    if logger is not None:
        logger.success("Loaded drugs, biocide and metal targets into the ontology.")
    

def remove_unused_subclasses_with_property(onto: Ontology, parent_cls: Thing, property_name: str, logger) -> None:
    """
    Remove all unused subclasses of a given parent class in the ontology that do not have a specific object property.

    Parameters
    ----------
    onto : Ontology
        The ontology object to remove class instances from.
    parent_cls : Thing
        The parent class whose unused subclasses will be removed.
    property_name : str
        The name of the object property to check.
    logger : loguru.logger, optional
        Logger object for logging messages, by default None
    """
    # List for storing unused subclasses
    unused_subclasses = []

    # Check if parent class even exists
    if parent_cls is None and logger :
        logger.error(f"Parent class {parent_cls} is not defined in the ontology.")
        return
    
    # Get all subclasses of the given parent class
    all_classes = list(parent_cls.subclasses())

    # Loop through them
    for subclass in all_classes:
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
    
    logger.info(f"Destroyed {len(unused_subclasses)}/{len(all_classes)} subclasses of {parent_cls.name}")


def gene_target(gene: Thing, og: Thing, target: str, onto: Ontology, db_name: str = None):
    """Assign a resistance target to a gene and its original gene instance in the ontology.

    Parameters
    ----------
    gene : Thing
        The gene instance
    og : Thing
        The original gene instance
    target : str
        The target to be assigned to the gene
    onto : Ontology
        The ontology object
    db_name : str, optional
        The name of the database, by default None

    Returns
    -------
    bool    
        True if target was found and assigned, False if not
    """

    # Replace spaces and dashes with underscores
    target = target.replace(" ", "_").replace("-", "_").title()

    # Get class from the ontology
    target_instance = None
    if '+' in target: # drug combination target
        target_instance = get_or_create_subclass(onto = onto, parent_cls=onto.AntibioticResistancePhenotype, subclass_name=target)
        target_instance.is_drug_combination.append(True)
        # Split the drug combination targets +
        for ts in target.split('+'):
            ts_instance = get_instance(onto = onto, name = ts.replace(" ", "_"))
            if ts_instance is not None:
                target_instance.is_a.append(ts_instance)
    else:
        target_instance = get_instance(onto = onto, name = target)
    
    # No match found!
    if target_instance is None:
        return False
    
    # Check if the target is a phenotype or a class
    type_check = {
        'phenotype': any(
            [
                onto.AntibioticResistancePhenotype in target_instance.is_a, 
                onto.Metal in target_instance.is_a, 
                onto.Biocide in target_instance.is_a, 
                onto.UnclassifiedResistance in target_instance.is_a
                ]
        ),
        'class': any(
            [
                onto.AntibioticResistanceClass in target_instance.is_a,
                onto.BiocideClass in target_instance.is_a,
                onto.MetalClass in target_instance.is_a,
                onto.UnclassifiedResistanceClass in target_instance.is_a
            ]
        )
    }

    # Assign the relations based on whether the drug target is a phenotype or a class object
    if type_check['phenotype']: 
        og.has_predicted_phenotype.append(target_instance)

        for target_relation in target_instance.is_a:
            if any([onto.AntibioticResistanceClass in target_relation.is_a, onto.BiocideClass in target_relation.is_a]):
                gene.has_resistance_class.append(target_relation)
                og.has_resistance_class.append(target_relation)
    elif type_check['class']: 
        gene.has_resistance_class.append(target_instance)
        og.has_resistance_class.append(target_instance)
    else:
        return False

    # Assign in which database the target annotation was encountered
    try:        
        db_instance = get_instance(onto = onto, name = db_name)
        target_instance.found_in.append(db_instance)
    except TypeError:
        pass
    
    return True

def reclassify_genes(onto: Ontology):
    """Reclassifies genes in the ontology based on their predicted phenotypes and resistance classes.

    Parameters
    ----------
    onto : Ontology
        The ontology object
    """

    # Get all PanGenes added to the ontology
    for gene in onto.search(type=onto.PanGene):
        # extract the associated resistance phenotypes and classes
        phenotypes = list(set(gene.has_predicted_phenotype))
        classes = list(set(gene.has_resistance_class))
        
        # If the gene has resistance to antimicrobials, its a antimicrobial resistance gene.
        if any([rt in ph.is_a for ph in phenotypes + classes for rt in [onto.AntibioticResistanceClass, onto.AntibioticResistancePhenotype]]):
            gene.is_a.append(onto.AntimicrobialResistanceGene)

        # If the gene has resistance to biocides, its a biocide resistance gene.
        if any([rt in ph.is_a for ph in phenotypes + classes for rt in [onto.BiocideClass, onto.Biocide]]):
            gene.is_a.append(onto.BiocideResistanceGene)

        # If the gene has resistance to metals, its a metal resistance gene.
        if any([rt in ph.is_a for ph in phenotypes + classes for rt in [onto.MetalClass, onto.Metal]]):
            gene.is_a.append(onto.MetalResistanceGene)
