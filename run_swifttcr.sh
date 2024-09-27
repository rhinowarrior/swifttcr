#!/bin/bash
##$ -l h_rt=03:00:00
#$ -cwd
#$ -o pipeline.out
#$ -e pipeline.err
#$ -V
#$ -pe smp 12

python3 scripts/swift_tcr.py -r example/input/5c0c/5c0c_pmhc.pdb -l example/input/5c0c/5c0c_tcr.pdb -o example/output/  -op 5c0c -c 12 -t 9