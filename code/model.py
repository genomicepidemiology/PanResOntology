from owlready2 import AnnotationProperty, Thing, FunctionalProperty

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
    class OriginalGene(PanGene): 
        description = "This is the original name extracted from the fasta header for each individual gene."

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

    # There are three types of resistances in PanRes: antibiotic, metal and biocide
    class ResistanceType(Resistance): pass
    class AntibioticResistance(ResistanceType): pass
    class BiocideResistance(ResistanceType): pass
    class MetalResistance(ResistanceType): pass

    # Antibiotic resistance - class, phenotype, mechanism
    class AntibioticResistanceClass(AntibioticResistance): pass
    class AntibioticResistancePhenotype(AntibioticResistance): pass
    class AntibioticResistanceMechanism(AntibioticResistance): pass

    # Metal resistance - metal
    class Metal(MetalResistance): pass

    # Biocide resistance - biocide
    class Biocide(BiocideResistance): pass

    # In PanRes 2.0, we are pivoting into proteins as well
    class Protein(Resistance):
        description = "Genes are translated into proteins, which was added as a new component of the PanRes database in version 2.0."

    '''
    Functional properties
    '''

    class has_length_bp(PanGene >> int): 
        namespace = onto

    class has_length_aa(Protein >> int):
        namespace = onto

    class dna_accession(PanGene >> str): 
        namespace = onto

    class protein_accession(PanGene >> str):
        namespace = onto

    class card_link(PanGene >> str): 
        namespace = onto

    class is_from_database(PanGene >> Database):
        namespace = onto

    class original_fasta_header(OriginalGene >> str): 
        namespace = onto

    class gene_alt_name(OriginalGene >> str): 
        namespace = onto

    class original_gene_is_from_database(OriginalGene >> Database):
        namespace = onto

    class has_predicted_phenotype(PanGene >> AntibioticResistancePhenotype):
        namespace = onto

    class original_has_predicted_phenotype(OriginalGene >> AntibioticResistancePhenotype):
        namespace = onto

    class has_resistance_class(PanGene >> AntibioticResistanceClass):
        namespace = onto
    class original_has_resistance_class(OriginalGene >> AntibioticResistanceClass):
        namespace = onto

    # class phenotype_is_class(AntibioticResistancePhenotype >> AntibioticResistanceClass):
    #     namespace = onto
    
    class gene_translated_to_protein(Gene >> Protein):
        namespace = onto

    class is_drug_combination(AntibioticResistancePhenotype >> bool): 
        namespace = onto

    class has_predicted_metal_resistance(PanGene >> Metal):
        namespace = onto

    class original_has_predicted_metal_resistance(OriginalGene >> Metal):
        namespace = onto

    class has_predicted_biocide_resistance(PanGene >> Biocide):
        namespace = onto
    class original_has_predicted_biocide_resistance(OriginalGene >> Biocide):
        namespace = onto
    class phenotype_found_in(AntibioticResistancePhenotype >> Database):
        namespace = onto

    class class_found_in(AntibioticResistanceClass >> Database):
        namespace = onto
    '''
    Inferred relationships
    '''

    class AntimicrobialResistanceGene(PanGene):
        namespace = onto
        equivalent_to = [
            PanGene & 
            has_resistance_class.some(AntibioticResistanceClass) &
            original_gene_is_from_database.some(Database)
        ]