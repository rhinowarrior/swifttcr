# SwiftTCR
## Overview
SwiftTCR is tool that predict bindings between TCR and peptide-MHC

## Features

## Running SwiftTCR
To get started on using SwiftTCR the easiest way is to clone or download this repository. After that you can move into the folder and should be able to run it using the command under here.

When in the SwiftTCR folder the following command should be used
```
python scripts/swift_tcr.py -r /your/input/peptide-mhc -l /your/input/tcr -o output_directory -op output-prefix
```
<br />

**Example command:**
```
python scripts/swift_tcr.py -r example/input/3w0w/3w0w_pmhc.pdb -l example/input/3w0w/3w0w_tcr.pdb -o example/output/ -op first_test
```

## Dependencies:
* [Pymol open source: 3.0.0](https://github.com/schrodinger/pymol-open-source)
* [anarci: 2021.02.04](https://github.com/oxpig/ANARCI)
* [gradpose: 0.1](https://github.com/X-lab-3D/GradPose)
* [pdb-tools: 2.5.0](http://www.bonvinlab.org/pdb-tools/)
* [ProDy: 2.4.1](https://github.com/prody/ProDy)
* [torch: 2.4.1](https://pytorch.org/)
* [pdb2sql: 0.5.3](https://github.com/DeepRank/pdb2sql)

## Chains of the output structures
The peptide-MHC chains will be named in the following way
* A = MHC
* B = b2m
* C = peptide

The TCR chains will be named in the following chains
* D = TCR alpha chain (IMGT numbering)
* E = TCR beta chain (IMGT numbering)