from owlready2 import get_ontology, sync_reasoner
import model
from databases import panres, resfinder, resfinderfg, card, megares, amrfinderplus, argannot, metalres
from targets import *

from loguru import logger
import sys  

logger.add("panres_messages.log")

onto = get_ontology("http://myonto.com/PanResOntology.owl")
logger.success(f"Created an empty ontology: {onto.base_iri}")

model.createModel(onto)
logger.success("Created the ontology model.")

# Define targets
load_targets(excelfile='data/targets.xlsx', onto=onto, logger=logger)

# Load into data from first version of PanRes
panres.add_panres_genes("data/PanRes_data_v1.0.0.tsv", onto, logger)

# Load data about ResFinder genes
resfinder.add_resfinder_annotations("data/resfinder_db/phenotypes.txt", onto, logger=logger)

card.add_card_annotations("data/aro_index.tsv", onto, logger=logger)

megares.add_megares_annotations(onto, logger=logger)

resfinderfg.add_resfinderfg_annotations("data/resfinderfg_anno.txt", onto, logger=logger)

amrfinderplus.add_amrfinderplus_annotations("data/ReferenceGeneCatalog.txt", onto, logger=logger)

argannot.add_argannot_annotations(onto, logger=logger)

metalres.add_metalres_annotations(onto, logger=logger)


# Remove unusued classes of AntimicrobialResistanceClass and AntimicrobialResistancePhenotype 
# remove_unused_subclasses_with_property(onto=onto, parent_cls=onto.AntibioticResistanceClass, property_name='has_resistance_class')
# remove_unused_subclasses_with_property(onto=onto, parent_cls=onto.AntibioticResistancePhenotype, property_name='has_predicted_phenotype')

sync_reasoner(debug=0, infer_property_values = True)

# Save the ontology to a file
ont_file = 'panres_v2.owl'
onto.save(file=ont_file, format="rdfxml")
logger.info(f"Saved ontology to file: {ont_file}.")