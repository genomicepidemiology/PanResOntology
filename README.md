# PanRes Ontology

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

### Classes
#### Class: `Resistance`
**Description**: Overall category for resistance in the PanRes database.

#### Class: `Gene`
**Description**: The major component of a reference sequence database: the genes. <br>
**Parent Class**: `Resistance`

#### Class: `PanGene`
**Description**: This is a gene identifier that follows the pan_ naming scheme. <br>
**Parent Class**: `Gene`

#### Class: `OriginalGene`
**Description**: This is the original name extracted from the fasta header for each individual gene.<br>
**Parent Class**: `Gene`

#### Class: `PanGeneCluster`
**Description**: Represents a cluster of PanGenes.<br>
**Parent Class**: `Gene`

#### Class: `Database`
**Description**: Represents a database from which genes are originally sourced.<br>
**Parent Class**: `Resistance`<br>
**Children**: `AMRFinderPlus`, `ARGANNOT`, `CARD`, `CsabaPal`, `MegaRes`, `MetalRes`, `ResFinder`, `ResFinderFG`, `BacMet`.<br>

#### Class: `ResistanceType`
**Description**: Represents different types of resistance.<br>
**Parent Class**: `Resistance`

#### Class: `AntibioticResistance`
**Description**: Represents antibiotic resistance.<br>
**Parent Class**: `ResistanceType` <br>
**Children**:
- Class: `AntibioticResistanceClass`
  - **Description**: Represents classes of antibiotic resistance.<br>
- Class: `AntibioticResistancePhenotype`
  - **Description**: Represents phenotypes of antibiotic resistance.<br>
- Class `AntibioticResistanceMechanism`
  - **Description**: Represents mechanisms of antibiotic resistance.<br>

#### Class: `BiocideResistance`
**Description**: Represents biocide resistance.<br>
**Parent Class**: `ResistanceType`
**Children**:
- Class: `BiocideClass`
  - **Description**: Represents classes of biocide resistance.<br>
- Class: `Biocide`
  - **Description**: Represents biocides involved in resistance.<br>

#### Class: `MetalResistance`
**Description**: Represents metal resistance.<br>
**Parent Class**: `ResistanceType`
**Children**:
- Class: `MetalClass`
  -  **Description**: Represents classes of metal resistance.<br>
- Class: `Metal`
  - **Description**: Represents metals involved in resistance.<br>


#### Class: `UnclassifiedResistanceClass`
**Description**: Represents unclassified resistance classes.<br>
**Parent Class**: `ResistanceType`
**Children**:
- Class: `UnclassifiedResistance`
  - **Description**: Represents unclassified resistance.<br>

#### Class: `Protein`
**Description**: Genes are translated into proteins, which was added as a new component of the PanRes database in version 2.0.
**Parent Class**: `Resistance`
**Children**:
- Class: `PanProtein`
  - **Description**: Represents a PanProtein.
- Class: `OriginalProtein`
  - **Description**: Represents an original protein.
- Class: `PanProteinCluster`
  - **Description**: Represents a cluster of PanProteins.

