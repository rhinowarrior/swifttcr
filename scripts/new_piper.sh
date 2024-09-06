#!/bin/bash
##$ -l h_rt=03:00:00
#$ -cwd
#$ -o pipeline.out
#$ -e pipeline.err
#$ -V
#$ -pe smp 30

## qsub -q all.q@narrativum.umcn.nl cluster_jobs/run_pipeline.sh experiments/generated_rotations_from_angle/input_new/6eqb_l_u_pnon.pdb experiments/generated_rotations_from_angle/input_new/6eqb_r_u_pnon.pdb /home/jaarts/pipeline/6eqb_60/ /home/jaarts/Restraints/restraintsDE.json /home/jaarts/experiments/filtered_cr_in_60.prm
## Example usage qsub -q all.q@narrativum.umcn.nl run_pipeline receptor ligand outputdir restraints rotations
export OMP_NUM_THREADS=2
##Prepare script
receptorpath=$1
ligandpath=$2
echo "${receptorpath}"
receptor=$(basename ${receptorpath})
echo "${receptor}"
receptorbase=$(basename "$receptor" .pdb)
echo "${receptorbase}"
receptorms="${receptorbase}.ms"
echo "${receptorms}"
ligand=$(basename ${ligandpath})
ligandbase=$(basename "$ligand" .pdb)
ligandms="${ligandbase}.ms"
outputdir=$3
restraints=$4
rotations=$5
out_prefix=$(basename ${outputdir})

# echo "placement start at $(date)"
# python3 scripts/initial_placement.py "${receptorpath}" "${ligandpath}" "${outputdir}" input/2bnr_l_u.pdb input/2bnr_r_u.pdb  && cd "${outputdir}"
# wait
# python3 /home/nils/stage/scripts/pdb2ms.py "${outputdir}"
# wait
# echo "piper start at $(date)"
# /home/nils/stage/piper/piper_attr -k1 --msur_k=1.0 --maskr=1.0 -T FFTW_EXHAUSTIVE -p /home/nils/stage/piper/atoms04.prm -f /home/nils/stage/piper/coeffs04.prm -r "${rotations}" "${outputdir}""${receptorms}" "${outputdir}""${ligandms}"
# wait
# python3 scripts/postfilter.py "${outputdir}"ft.000.00 "${rotations}" "${restraints}" "${outputdir}""${receptor}" "${outputdir}""${ligand}" "${outputdir}"filtered_ft.000.00 && mkdir "${outputdir}"rotated && cd "${outputdir}"rotated
# wait
# echo pwd
# echo "python3 scripts/apply_results.py --limit 1000" "${outputdir}" "ft.000.00 ""${rotations}" "${outputdir}""${ligand}"" --out-prefix" "${out_prefix}""_r && mkdir" "${outputdir}""merged"
# python3 scripts/apply_results.py --limit 1000 "${outputdir}"filtered_ft.000.00 "${rotations}" "${outputdir}""${ligand}" --out-prefix "${out_prefix}"_r && mkdir "${outputdir}"merged 

# wait
# python3 scripts/merge_pdbs.py "${outputdir}""${receptor}" "${outputdir}"rotated/ "${outputdir}"merged/
#python3 /home/jaarts/python_programs/merge_pdbs.py /home/jaarts/pipeline/1ao7/1ao7_l_u_pnon.pdb /home/jaarts/pipeline/1ao7/rotated/ /home/jaarts/pipeline/1ao7/merged/
wait
python scripts/pairwise_rmsd_backup.py "${outputdir}"merged "${outputdir}"irmsd.csv A D interface 2
wait

# python3 scripts/clustering.py "${outputdir}"irmsd.csv > "${outputdir}"clustering.txt


