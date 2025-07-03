# SAVI-Space

![Title Image](images/SAVISpace.png)

Welcome to the **SAVI-Space** repository, which contains the Python-based code and tools for the preparation of building blocks used in the SAVI-Space project. This repository focuses on preprocessing, validation, and preparation of molecular fragment space generation.

## CHMTRN/PATRAN parser

The CHMTRN/PATRAN parser is a tool that reads the CHMTRN/PATRAN freaction pattern and converts them into reaction SMARTS. The tool is written in Python.

The tool is able to parse the CHMTRN/PATRAN string and convert it into a reaction SMARTS string:
```bash
python src/transpiler/transpiler.py <PARTAN/CHMTRN string>
```
or alternatively a json file containing the CHMTRN/PATRAN string:
```bash
python src/transpiler/transpiler.py -f <json file> -o <output file>
```

## Space creation


The code is organized into three main components:

---

### 1. **Building Block Molecule Standardization and Filtering Tool**

This module handles the preprocessing of molecular building blocks, ensuring consistency and quality. The key functionalities include:

#### **Standardization**
- Automatically cleans and standardizes molecular structures.
- Normalizes charges, removes invalid atoms, and ensures uniformity by removing stereochemistry and atom mappings.

#### **File Format Conversion**
- Converts `.sdf` files to `.smi` format for enhanced compatibility with downstream tools.

#### **Validation**
- Filters valid molecules using the RDKit library.
- Offers optional validation using an external `naomi_check` binary.

#### **Disconnected SMILES Fragmentation**
- Decomposes disconnected SMILES strings into individual fragments for independent analysis.

#### **Filtering Options**
- Filters molecules based on:
  - **Heavy Atom Count:** `min_heavy_atoms`, `max_heavy_atoms`.
  - **Molecular Weight:** `min_molecular_weight`, `max_molecular_weight`.
  - **Bertz Complexity Index:** `max_bertz_ct`.

#### **Parallel Processing**
- Leverages multi-core processing for faster execution.

#### **Duplicate Removal**
- Groups and deduplicates SMILES entries, combining associated names.

#### **Logging**
- Provides detailed logs for every processing step, ensuring transparency and traceability.

---

### 2. **Reactant (Building Block) Validation Tool**

This component validates building blocks against specified criteria to ensure their suitability for reaction-based applications. Its main features include:

#### **SMARTS Matching**
- Matches reactants to predefined SMARTS patterns to identify valid reaction sites.
- Includes advanced options to:
  - Exclude exocyclic double bonds.
  - Ensure unique matches across reactants.

#### **Kill Logic**
- Implements customizable "kill logic" to filter out unwanted reactants using SMARTS-based rules.

#### **Extensible Input/Output**
- Accepts SMILES strings as input.
- Outputs detailed CSV files with:
  - Match statistics.
  - Filtered results based on "kill" criteria.

---

### 3. **Space Generation Preparation**

This module prepares preprocessed building blocks for space generation, catering to two distinct toolkits:

**Colibri Toolkit and Naomi Toolkit**

---

## Getting Started

### **Prerequisites**
- Python 3.x
- Required libraries:
  - RDKit
  - Pandas
  - NumPy
  - Additional dependencies detailed in `requirements.txt`.

### **Installation**
1. Clone the repository:
   ```bash
   git clone https://github.com/rareylab/SAVI-Space.git
   ```
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

### **Usage**
#### **Building Block Molecule Standardization and Filtering Tool**

This tool standardizes and filters building blocks in `.sdf` or `.smi` format. To use the tool:

Run the standardization script:
   ```bash
   python src/standardize_bbs.py --input building_blocks.sdf
   ```

The standardized output will be saved as `building_blocks_standardized.smi`.

The filtering options can be customized using command-line arguments:
e.g., filtering by heavy atom count of minimum 10 and a maximum molecular weight of 500:
   ```bash
   python src/standardize_bbs.py --input building_blocks.sdf --min_heavy_atoms 10 --max_molecular_weight 500
   ```
A coma-separated list of filtering values is stored as `building_blocks.csv`.

Moreover, the tool can split disconnected SMILES strings into fragments:
   ```bash
   python src/standardize_bbs.py --input building_blocks.sdf --divide_disconnected_smiles
   ```
   The output will be saved as `building_blocks_fragments.smi` and `building_blocks_fragments_standardized.smi`.


#### **Reactant (Building Block) Validation Tool**

This tool validates building blocks against predefined SMARTS patterns and kill logic. To use the tool:

Run the validation script:
   ```bash
   python src/space_builder/reactants_validation.py -n SAVI-Space-2024 -r building_blocks.smi --apply_kill
   ```
The output will be saved at `SAVI-Space-2024/`. You can specify a path using the `-o` argument.
   ```bash
   python src/space_builder/reactants_validation.py -n SAVI-Space-2024 -r building_blocks.smi --apply_kill -o path/to/output
   ```
   Will be saved at `path/to/output/SAVI-Space-2024/`.

Additionaly the match logics can be customized using the `--unify_matches`:
- Unify matches on the same atoms.

or the `--unique_matches` flag:
- Unify matches based on the given indices.

#### **Space Generation Preparation**

This module prepares building blocks and reactions for space generation using the Colibri or Naomi toolkit. To use the tool:

- For Colibri:
   ```bash
   python src/space_builder/Colibri_space_creation.py -n SAVI-Space-2024
   ```

- For Naomi:
   ```bash
   python src/space_builder/Naomi_space_creation.py -n SAVI-Space-2024
   ```
A bash script will be generated to run the toolkit.
In the script, the paths to the executables should be specified.

---

## Legal Note

This material is part of the SAVI-Space project, created by Malte Korn and Matthias Rarey at the Center for Bioinformatics, University of Hamburg, with valuable support from Marc Nicklaus (NIH), Phil Judson(Self-employed), Raphael Klein (BioSolveIT GmbH), and Christian Lemmen (BioSolveIT GmbH).

Philip Judson and Marc Nicklaus provided the LHASA transform rules and contributed their expertise regarding SAVI-2020.

The SAVI-Space project, including all its components and this file/directory, is licensed under the Creative Commons Attribution 4.0 International License (CC-BY 4.0). Under this license, you are free to share, copy, and redistribute the material in any medium or format, as well as adapt, transform, and build upon the material, under the following terms:

    Attribution: You must give appropriate credit, provide a link to the license, and indicate if changes were made.

The building blocks used within this project are sourced from Enamine Ltd..

This work is based on SAVI-2020, as described in the publication: https://doi.org/10.1038/s41597-020-00727-4.

For further information about SAVI-Space or licensing, please contact the Center for Bioinformatics, University of Hamburg.

---

## License

This project is licensed under the `CC-BY 4.0` License. See the `LICENSE` file for details.

---

## References

Publication:
Korn, M., Judson, P., Klein, R. et al. SAVI Space—combinatorial encoding of the billion-size synthetically accessible virtual inventory. Sci Data 12, 1064 (2025). https://doi.org/10.1038/s41597-025-05384-z

Dataset:
Korn, Malte. (2025). SAVI-Space (Version 1.0.0) [Data set]. http://doi.org/10.25592/uhhfdm.15990
