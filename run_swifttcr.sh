#!/bin/bash
##$ -l h_rt=03:00:00
#$ -cwd
#$ -o pipeline.out
#$ -e pipeline.err
#$ -V
#$ -pe smp 12

python3 scripts/swift_tcr.py -r /home/nils/swifttcr/example/input/6amu/6amu_pmhc.pdb -l /home/nils/swifttcr/example/input/6amu/6amu_tcr.pdb -o example/output/  -op bash_test