# SwiftTCR: Efficient computational docking protocol of TCRpMHC-I complexes using restricted rotation matrices

## Overview
**SwiftTCR** is a fast fourier transform based rigid-body docking tool designed to predict bindings between T-cell receptors (TCR) and peptide-MHC complexes.

Link to the paper: [SwiftTCR](https://www.biorxiv.org/content/10.1101/2024.05.27.596020v2.full)

## Features
- Predict binding interactions between TCRs and peptide-MHC.
- User-friendly command-line interface.
- Efficient clustering algorithms for data analysis.
- With 12 CPU cores it takes around 260 seconds.

## Getting Started

To get started with SwiftTCR, follow these steps:

1. **Clone or Download** this repository.
2. Navigate into the SwiftTCR folder.

### Installation

To quickly install all the necessary packages, you can use the provided `swifttcr_install.yml` file. Run the following commands:

```
conda env create -f swifttcr_install.yml
conda activate swifttcr 
```

### Running SwiftTCR
Use the following command to execute SwiftTCR:

```bash
python3 scripts/swift_tcr.py -r /your/input/peptide-mhc -l /your/input/tcr -o output_directory -op output_prefix -c number_of_cores -t clustering_threshold (default=9) -m amount_of_models_generated
```
<br />

**Example command:**
```bash
python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/3w0w/3w0w_pmhc_renumbered.pdb -l example/input/benchmark_dataset/3w0w/3w0w_tcr.pdb -o example/output/ -op first_test -c 6 -t 9 -m 100
```

## Dependencies:
* Python 3.9.12 and python 2.7
* [Pymol open source: 3.0.0](https://github.com/schrodinger/pymol-open-source)
* [anarci: 2021.02.04](https://github.com/oxpig/ANARCI) 
* [gradpose: 0.1](https://github.com/X-lab-3D/GradPose)
* [pdb-tools: 2.5.0](http://www.bonvinlab.org/pdb-tools/)
* [torch: 2.4.1](https://pytorch.org/)
* [pdb2sql: 0.5.3](https://github.com/DeepRank/pdb2sql)
* [Biopython: 1.84](https://biopython.org/)

## Output SwiftTCR

### Output Structure Naming Convention initial placement

#### Peptide-MHC Chains
- **A** = MHC (Not IMGT numbered)
- **B** = β2m (Not IMGT numbered)
- **C** = Peptide (Not IMGT numbered)

#### TCR Chains
- **D** = TCR Alpha Chain (IMGT numbered)
- **E** = TCR Beta Chain (IMGT numbered)

### Output Structure Naming Convetion SwiftTCR

#### TCR-Peptide-MHC
- **A** = The ABC Chains of original Peptide-MHC combined
- **D** = The Alpha and Beta chains of TCR combined

### Structure of output folder
The output is a folder, named using the specified output prefix, created within the designated output directory. This folder contains the following files and subfolders:
```
output
    └── 3w0w
        ├── 3w0w
        ├── 3w0w_pmhc_renumbered_pnon.ms
        ├── 3w0w_pmhc_renumbered_pnon.pdb
        ├── 3w0w_pmhc_renumbered_pnon_rename.pdb
        ├── clustering.txt
        ├── ft.000.00
        ├── ft.001.00
        ├── ft.002.00
        ├── ft.003.00
        ├── ft.004.00
        ├── ft.005.00
        ├── ft.006.00
        ├── ft.007.00
        ├── irmsd.csv
        ├── merged
        │   └── merged_0.pdb
        ├── renumbered_3w0w_tcr.pdb
        ├── renumbered_3w0w_tcr_pnon.ms
        ├── renumbered_3w0w_tcr_pnon.pdb
        └── rotated
            └── 3w0w.0.pdb
```

#### Merged folder
A Folder that contains the predicted structures of the TCR-peptide-MHC.

#### Rotated folder
A folder that contains the rotated structures that where made with the use of the PIPER energies.

#### irmsd.csv
A .csv file that contains the iRMSD between all the merged files.

#### ft. files
These files contain the energies calculated by PIPER and are sorted on lowest energies first.

#### pnon.pdb files
Input structures that have been prepared and been aligned to a reference structure. They also have the same chainIDs as the reference structures.

#### .ms files
These files have added attractions that is used to run PIPER.

#### clustering.txt
This file is the output that contains the top ranked models and how many neigbours that where found with near the model.
---

# Installing Python 2.7 and Python 3.12.9 in a Conda Environment

Because of `prepare.py`, this program requires both Python 2.7 and Python 3.12.9. Follow these steps to install and configure both versions.

## Step-by-Step Installation Guide

### 1. Install Conda (if not already installed)

If you don't have Conda installed, you can download and install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution).

### 2. Create a Conda Environment for Python 3.12.9

Run the following command to create a new Conda environment named `swifttcr` with Python 3.12.9:

conda create -n swifttcr python=3.12.9

### 3. Activate the Conda Environment

Activate the newly created environment:

```bash
conda activate swifttcr
```

### 4. Install Required Packages

Install any required packages for your project within this environment. For example, if you need `numpy` and `pandas`, you can do so with:

```bash
conda install numpy pandas
```

### 5. Install Python 2.7

To install Python 2.7, you can use the following commands:

#### 5.1. Download Python 2.7

Run this command to download the Python 2.7 source code:

```bash
wget https://www.python.org/ftp/python/2.7.18/Python-2.7.18.tgz
```

#### 5.2. Extract the Source Code

Extract the downloaded tarball:

```bash
tar -xvf Python-2.7.18.tgz
```

#### 5.3. Move to the Python Directory

Change into the newly created Python directory:

```bash
cd Python-2.7.18
```

#### 5.4. Configure the Installation

Run the configure script with the `--prefix` option to specify the installation directory:

```bash
./configure --prefix=$HOME/.local
```

#### 5.5. Build and Install Python 2.7

Compile and install Python 2.7 with:

```bash
make
make install
```

### 6. Add Python 2.7 to Your PATH

To make Python 2.7 accessible from the command line, add it to your PATH by running:

```bash
echo 'export PATH=$HOME/.local/bin:$PATH' >> ~/.bashrc
```

### 7. Reload Your Profile

Apply the changes to your current session:

```bash
source ~/.bashrc
```

### 8. Check Installations

Verify the installations by checking the versions:

```bash
python2.7 --version
python3 --version  # This will point to Python 3.12.9 in the active Conda environment.
```