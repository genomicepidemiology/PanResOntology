from owlready2 import get_ontology, sync_reasoner
import model
from databases import panres, resfinder, resfinderfg, card, megares, amrfinderplus, argannot

from loguru import logger
import sys  

logger.add("panres_messages.log")

onto = get_ontology("http://myonto.com/PanResOntology.owl")
logger.success(f"Created an empty ontology: {onto.base_iri}")

model.createModel(onto)
logger.success("Created the ontology model.")

panres.add_panres_genes("data/PanRes_data_v1.0.0.tsv", onto, logger)

resfinder.add_resfinder_annotations("data/resfinder_db/phenotypes.txt", onto, targetfile='data/targets.xlsx', logger=logger)


card.add_card_annotations("data/aro_index.tsv", onto, targetfile='data/targets.xlsx', logger=logger)

megares.add_megares_annotations(onto, targetfile='data/targets.xlsx', logger=logger)

resfinderfg.add_resfinderfg_annotations("data/resfinderfg_anno.txt", onto, targetfile='data/targets.xlsx', logger=logger)

amrfinderplus.add_amrfinderplus_annotations("data/ReferenceGeneCatalog.txt", onto, targetfile='data/targets.xlsx', logger=logger)

argannot.add_argannot_annotations(onto, targetfile='data/targets.xlsx', logger=logger)


sync_reasoner(debug=0, infer_property_values = True)

# Save the ontology to a file
ont_file = 'panres_v2.owl'
onto.save(file=ont_file, format="rdfxml")
logger.info(f"Saved ontology to file: {ont_file}.")