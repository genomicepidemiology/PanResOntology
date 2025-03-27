### Classes
#### Class: `Resistance`
- **Description**: Overall category for resistance in the PanRes database.

#### Class: `Gene`
- **Description**: The major component of a reference sequence database: the genes. 
- **Parent Class**: `Resistance`
- **Children**:
  - Class: `PanGene`
    - **Description**: This is a gene identifier that follows the pan_ naming scheme. 
  - Class: `OriginalGene`
    - **Description**: This is the original name extracted from the fasta header for each individual gene.
  - Class: `PanGeneCluster`
    - **Description**: Represents a cluster of PanGenes.

#### Class: `Database`
- **Description**: Represents a database from which genes are originally sourced.
- **Parent Class**: `Resistance`
- **Children**: 
  - Class: `AMRFinderPlus`
  - Class: `ARGANNOT`
  - Class: `CARD`
  - Class: `CsabaPal`
  - Class: `MegaRes`
  - Class: `MetalRes`
  - Class: `ResFinder`
  - Class: `ResFinderFG`
  - Class: `BacMet`

#### Class: `ResistanceType`
- **Description**: Represents different types of resistance.
- **Parent Class**: `Resistance`
- **Children**:
  - Class: `AntibioticResistance`
    - **Description**: Represents antibiotic resistance.
    - **Children**:
      - Class: `AntibioticResistanceClass`
        - **Description**: Represents classes of antibiotic resistance.
      - Class: `AntibioticResistancePhenotype`
        - **Description**: Represents phenotypes of antibiotic resistance.
      - Class `AntibioticResistanceMechanism`
        - **Description**: Represents mechanisms of antibiotic resistance.
  - Class: `BiocideResistance`
    - **Description**: Represents biocide resistance.
    - **Children**:
      - Class: `BiocideClass`
        - **Description**: Represents classes of biocide resistance.
      - Class: `Biocide`
        - **Description**: Represents biocides involved in resistance.
  - Class: `MetalResistance`
    - **Description**: Represents metal resistance.
    - **Children**:
    - Class: `MetalClass`
      -  **Description**: Represents classes of metal resistance.
    - Class: `Metal`
      - **Description**: Represents metals involved in resistance.
  - Class: `UnclassifiedResistanceClass`
    - **Description**: Represents unclassified resistance classes.
    - **Children**:
      - Class: `UnclassifiedResistance`
        - **Description**: Represents unclassified resistance.

#### Class: `Protein`
- **Description**: Genes are translated into proteins, which was added as a new component of the PanRes database in version 2.0.
- **Parent Class**: `Resistance`
- **Children**:
  - Class: `PanProtein`
    - **Description**: Represents a PanProtein.
  - Class: `OriginalProtein`
    - **Description**: Represents an original protein.
  - Class: `PanProteinCluster`
    - **Description**: Represents a cluster of PanProteins.

### Functional Properties

#### Property: `has_length` [int]
- **Description**: Annotation property to specify the length of PanGenes and PanProteins.
- **Domain**: `PanGene`, `PanProtein`

#### Property: `accession` [str]
- **Description**: Annotation property to specify the accession number of genes and proteins.
- **Domain**: `PanGene`, `OriginalGene`, `Protein`

#### Property: `accession` [str]
- **Description**:
- **Domain**: