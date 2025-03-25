# PanRes Ontology
This GitHub repository contains the "ontologisation" of the annotations for the genes and proteins part of the PanRes collection. 

## Building The Ontology
In order to build the ontology with the genes included in the first version of PanRes, make sure to have the data files downloaded as well (see [data/](/data)).

To built the ontology from scratch, the following libaries are required:
```
owlready2==0.44
loguru==0.7.3
pandas==1.5.3
numpy==1.24.3
```

then run the script:
```
python code/PanResOntology.py
```

## PanRes API reference
The module in [model.py](/code/model.py) defines the ontology schema for the PanRes database using `owlready2`. It includes classes for various types of resistance genes, databases, and resistance types, as well as functional properties to describe relationships and attributes.

Extended documentation on the classes and properties are described in the [model.md](model.md) file.