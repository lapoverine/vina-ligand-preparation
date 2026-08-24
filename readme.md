# Ligand Preparation for AutoDock Vina

Solution for the task of the BIOCAD computational chemistry internship test case: batch preparation of small-molecule ligands for AutoDock Vina.

The script converts ligands from SMI or SDF files to individual PDBQT files using RDKit and Meeko.

## Features

* reads multiple ligands from `.smi` and `.sdf` files;
* preserves existing 3D coordinates from SDF files;
* generates 3D coordinates with RDKit ETKDGv3 when the input contains only 2D coordinates or SMILES;
* adds explicit hydrogen atoms when necessary;
* prepares ligands for AutoDock Vina with Meeko;
* writes one `.pdbqt` file per ligand using the ligand name from the input file;
* allows non-aromatic rings containing 6 or more atoms to be treated as flexible during ligand preparation;
* reports ligands that could not be read or converted without stopping processing of the remaining molecules.

## Environment

The Conda environment used for the project is provided in `environment.yml`.

To recreate it:

```bash
conda env create -f environment.yml
conda activate ligand-preparation
```

## Usage

The script accepts one SMI or SDF file containing multiple ligands.

### SMI input

```bash
python prepare_ligands.py vina_data/example.smi -o output_smi
```

### SDF input

```bash
python prepare_ligands.py vina_data/example.sdf -o output_sdf
```

The output directory is created automatically if it does not exist.

For example, a ligand named:

```text
CHEMBL4959907
```

is written as:

```text
CHEMBL4959907.pdbqt
```

## Ligand preparation

For molecules that already contain 3D coordinates, the existing geometry is retained and hydrogen atoms are added if required.

For SMILES and 2D SDF structures, a 3D conformer is generated using RDKit's ETKDGv3 method.

The resulting molecule is then prepared with Meeko and written in PDBQT format.

To implement the optional part of the test task, Meeko is configured with:

```python
MoleculePreparation(
    rigid_macrocycles=False,
    min_ring_size=6,
)
```

This enables flexible treatment of suitable non-aromatic rings containing more than five atoms.

## Project structure

```text
.
├── prepare_ligands.py
├── environment.yml
├── vina_data/
│   ├── example.sdf
│   └── example.smi
├── output_sdf/
├── output_smi/
└── README.md
```

`vina_data` contains the example input files provided with the original test task. The output directories contain the generated PDBQT files.
