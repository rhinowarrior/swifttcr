#!/bin/bash
##$ -l h_rt=03:00:00
#$ -cwd
#$ -o pipeline.out
#$ -e pipeline.err
#$ -V
#$ -pe smp 12

python scripts/swift_tcr.py -r example/input/benchmark_dataset/3w0w/3w0w_pmhc_renumbered.pdb -l example/input/benchmark_dataset/3w0w/3w0w_tcr.pdb -o example/output/ -op first_test -c 6 -t 9
