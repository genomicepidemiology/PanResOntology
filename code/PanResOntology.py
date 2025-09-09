from owlready2 import get_ontology, sync_reasoner
import model
from databases import panres, resfinder, resfinderfg, card, megares, amrfinderplus, argannot, metalres, bacmet, csabapal
from targets import *

from loguru import logger

logger.add("panres_messages.log")

onto = get_ontology("http://genepi.dk/PanResOntology.owl")
logger.success(f"Created an empty ontology: {onto.base_iri}")

model.createModel(onto)
logger.success("Created the ontology model.")

# Define targets
load_targets(excelfile='data/targets.xlsx', onto=onto, logger=logger)

# Load into data from first version of PanRes
panres.add_panres_genes("data/PanRes_data_v1.0.0.tsv", onto = onto, logger = logger, discarded='data/discarded/panres_removed_headers.txt')

# add proteins
panres.add_panres_proteins(file = 'data/proteins/panres_final_protein.faa', clstrs='data/proteins/panres_final_protein_50_90.faa.clstr', onto = onto, logger = logger)

# Load data about ResFinder genes
resfinder.add_resfinder_annotations("data/phenotypes.txt", onto, logger=logger)

# Load data about CARD genes
card.add_card_annotations("data/aro_index.tsv", onto, logger=logger)

# Load data about MegaRes genes
megares.add_megares_annotations(onto, logger=logger, mappingfile='data/megares_to_external_header_mappings_v3.00.csv')

# Load data about ResFinderFG genes
resfinderfg.add_resfinderfg_annotations("data/resfinderfg_anno.txt", onto, logger=logger)

# Load data about AMRFinderPlus genes
amrfinderplus.add_amrfinderplus_annotations("data/ReferenceGeneCatalog.txt", onto, logger=logger)

# Load data about ARGANNOT genes
argannot.add_argannot_annotations(onto, logger=logger)

# Load data about MetalRes genes
metalres.add_metalres_annotations(onto, logger=logger)

# Load data about BacMet genes
bacmet.add_bacmet_annotations(onto, mappingfile='data/BacMet_EXP.704.mapping.txt', logger=logger)

# Load data about CsabaPal genes
csabapal.add_csabapal_annotations(onto = onto, file='data/QSX_607_CsabaPal_metagenomics_metadata_final_notfiltered.csv', logger=logger)

# Remove unusued classes of AntimicrobialResistanceClass and AntimicrobialResistancePhenotype 
remove_unused_subclasses_with_property(onto=onto, parent_cls=onto.AntibioticResistanceClass, property_name='has_resistance_class', logger=logger)
remove_unused_subclasses_with_property(onto=onto, parent_cls=onto.AntibioticResistancePhenotype, property_name='has_predicted_phenotype', logger=logger)
remove_unused_subclasses_with_property(onto=onto, parent_cls=onto.Metal, property_name='has_predicted_phenotype', logger=logger)
remove_unused_subclasses_with_property(onto=onto, parent_cls=onto.Biocide, property_name='has_predicted_phenotype', logger=logger)
remove_unused_subclasses_with_property(onto=onto, parent_cls=onto.BiocideClass, property_name='has_resistance_class', logger=logger)
remove_unused_subclasses_with_property(onto=onto, parent_cls=onto.UnclassifiedResistance, property_name='has_predicted_phenotype', logger=logger)
remove_unused_subclasses_with_property(onto=onto, parent_cls=onto.ResistanceMechanism, property_name='has_mechanism_of_resistance', logger=logger)

reclassify_genes(onto)

logger.info("Syncing ontology reasonings..")
sync_reasoner(debug=0, infer_property_values = True)
# sync_reasoner_pellet(infer_property_values = False, infer_data_property_values = False, debug=0)

# Save the ontology to a file
ont_file = 'ontology/panres_v2.owl'
onto.save(file=ont_file, format="rdfxml")
logger.info(f"Saved ontology to file: {ont_file}.")