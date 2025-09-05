from argparse import Namespace
from git import Object
from owlready2 import AnnotationProperty, Thing, FunctionalProperty, ObjectProperty, Or, AllDisjoint


def createModel(onto):
    '''
    Basic class definitions
    '''
    
    # To add descriptions to ontologies
    class description(AnnotationProperty): pass

    # There are one overall category: Resistance 
    # since its the PanRes database!

    class Resistance(Thing):
        namespace = onto
        description = "Overall category for resistance in the PanRes database"
    
    # And we are focusing on resistance genes
    class Gene(Resistance):
        description = "The major component of a reference sequence database: the genes."

    # Two categories of resistance genes:
    # panres names - PanRes identifier
    # original gene names - those from the various databases
    class PanGene(Gene): 
        description = "This is a gene identifier that follows the pan_ naming scheme."
    class OriginalGene(Gene): 
        description = "This is the original name extracted from the fasta header for each individual gene."

    class DiscardedPanGene(PanGene):
        description = "An identifier for PanGenes that have been discarded from the current version of the PanRes database"

    class PanGeneCluster(Gene): pass

    # Need a disjoint relationship between PanGene and OriginalGene individuals
    AllDisjoint([PanGene, OriginalGene])

    # Also working with databases - to trace where each gene is originally from
    class Database(Resistance): pass
    class AMRFinderPlus(Database): pass
    class ARGANNOT(Database): pass
    class CARD(Database): pass
    class CsabaPal(Database): pass
    class MegaRes(Database): pass
    class MetalRes(Database): pass
    class ResFinder(Database): pass
    class ResFinderFG(Database): pass
    class BacMet(Database): pass

    # There are three types of resistances in PanRes: antibiotic, metal and biocide
    class ResistanceType(Resistance): pass
    class AntibioticResistance(ResistanceType): pass
    class BiocideResistance(ResistanceType): pass
    class MetalResistance(ResistanceType): pass

    # Antibiotic resistance - class, phenotype, mechanism
    class AntibioticResistanceClass(AntibioticResistance): pass
    class AntibioticResistancePhenotype(AntibioticResistance): pass

    # Metal resistance - metal
    class MetalClass(MetalResistance): pass
    class Metal(MetalResistance): pass

    # Biocide resistance - biocide
    class BiocideClass(BiocideResistance): pass
    class Biocide(BiocideResistance): pass

    # Unclassified ?
    class UnclassifiedResistanceClass(ResistanceType): pass
    class UnclassifiedResistance(UnclassifiedResistanceClass): pass

    class ResistanceMechanism(ResistanceType): pass

    # In PanRes 2.0, we are pivoting into proteins as well
    class Protein(Resistance):
        description = "Genes are translated into proteins, which was added as a new component of the PanRes database in version 2.0."
    class Structure(Protein):
        description = "This is a class for 3D structures, which can be used to represent protein structures or other relevant structures in the context of resistance."

    class PanProteinCluster(Protein): pass
    class PanProtein(Protein): pass

    class PanStructureCluster(Structure): pass
    class PanStructure(Structure): pass

    '''
    Functional properties
    '''

    class has_length(AnnotationProperty): 
        domain = [Or([PanGene, PanProtein])]
        range = [int]
        namespace = onto
    class accession(AnnotationProperty): 
        domain = [Or([PanGene, OriginalGene, Protein])]
        range = [str]
        namespace = onto
    
    class pubmed(AnnotationProperty): 
        domain = [Or([PanGene, OriginalGene, Protein])]
        range = [str]
        namespace = onto

    class card_link(AnnotationProperty): 
        domain = [Or([PanGene, OriginalGene])]
        range = [str]
        namespace = onto

    class is_from_database(ObjectProperty): #PanGene >> Database):
        domain = [Or([PanGene, OriginalGene])]
        range = [Database]
        namespace = onto
    
    class member_of(ObjectProperty):
        domain = [Or([PanGene, PanProtein, PanStructure])]
        range = [Or([PanGeneCluster, PanProteinCluster, PanStructureCluster])]
        namespace = onto

    class has_members(AnnotationProperty):
        domain = [Or([PanGeneCluster, PanProteinCluster, PanStructureCluster])]
        range = [str]
        namespace = onto

    class original_fasta_header(AnnotationProperty):
        domain = [OriginalGene] 
        range = [str]
        namespace = onto

    class same_as(ObjectProperty):
        domain = [PanGene]
        range = [OriginalGene]
        namespace = onto
    class has_pan_name(ObjectProperty): 
        domain = [OriginalGene]
        range = [PanGene]
        namespace = onto

    has_pan_name.inverse_property = same_as
    class gene_alt_name(OriginalGene >> str): 
        namespace = onto

    class has_predicted_phenotype(ObjectProperty):
        domain = [Or([PanGene, OriginalGene])]
        range = [Or([AntibioticResistancePhenotype, Biocide, Metal, UnclassifiedResistance])]
        namespace = onto
        class_property_type = ["some"]
    class has_resistance_class(ObjectProperty):
        domain = [Or([PanGene, OriginalGene, AntibioticResistancePhenotype])]
        range = [Or([AntibioticResistanceClass, BiocideClass, MetalClass, UnclassifiedResistanceClass])]
        namespace = onto
        class_property_type = ["some"]

    class has_mechanism_of_resistance(ObjectProperty):
        domain = [Or([PanGene, OriginalGene])]
        range = [ResistanceMechanism]
        namespace = onto
        class_property_type = ["some"]
    class translates_to(ObjectProperty): #(Gene >> Protein):
        domain = [PanGeneCluster]
        range = [Protein]
        namespace = onto

    class folds_to(ObjectProperty): #(Protein >> Structure):
        domain = [PanProteinCluster]
        range = [PanStructure]
        namespace = onto

    class is_drug_combination(AnnotationProperty): #(AntibioticResistancePhenotype >> bool): 
        domain = [Or([AntibioticResistancePhenotype, Metal, Biocide])]
        range = [bool]
        namespace = onto
    class metal_symbol(Metal >> str):
        namespace = onto
    class metal_comment(Metal >> str):
        namespace = onto
    class found_in(ObjectProperty):
        domain = [Or([AntibioticResistancePhenotype, Metal, Biocide, AntibioticResistanceClass, MetalClass, BiocideClass])]
        range = [Database]
        namespace = onto

    class AntimicrobialResistanceGene(PanGene): pass
    class BiocideResistanceGene(PanGene): pass
    class MetalResistanceGene(PanGene): pass