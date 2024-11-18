#!/bin/bash
##$ -l h_rt=03:00:00
#$ -cwd
#$ -o pipeline.out
#$ -e pipeline.err
#$ -V
#$ -pe smp 20

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/3w0w/3w0w_pmhc_renumbered.pdb -l example/input/benchmark_dataset/3w0w/3w0w_tcr.pdb -o example/output/ -op first_test -c 20 -t 3 -m 100
