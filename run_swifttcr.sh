#!/bin/bash
##$ -l h_rt=03:00:00
#$ -cwd
#$ -o pipeline.out
#$ -e pipeline.err
#$ -V
#$ -pe smp 12

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/1ao7/1ao7_pmhc_renumbered.pdb -l example/input/benchmark_dataset/1ao7/1ao7_tcr.pdb -o example/output -op 1ao7 -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/1mi5/1mi5_pmhc_renumbered.pdb -l example/input/benchmark_dataset/1mi5/1mi5_tcr.pdb -o example/output -op 1mi5 -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/1mwa/1mwa_pmhc_renumbered.pdb -l example/input/benchmark_dataset/1mwa/1mwa_tcr.pdb -o example/output -op 1mwa -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/1oga/1oga_pmhc_renumbered.pdb -l example/input/benchmark_dataset/1oga/1oga_tcr.pdb -o example/output -op 1oga -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/2bnr/2bnr_pmhc_renumbered.pdb -l example/input/benchmark_dataset/2bnr/2bnr_tcr.pdb -o example/output -op 2bnr -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/2ckb/2ckb_pmhc_renumbered.pdb -l example/input/benchmark_dataset/2ckb/2ckb_tcr.pdb -o example/output -op 2ckb -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/2nx5/2nx5_pmhc_renumbered.pdb -l example/input/benchmark_dataset/2nx5/2nx5_tcr.pdb -o example/output -op 2nx5 -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/2oi9/2oi9_pmhc_renumbered.pdb -l example/input/benchmark_dataset/2oi9/2oi9_tcr.pdb -o example/output -op 2oi9 -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/2pye/2pye_pmhc_renumbered.pdb -l example/input/benchmark_dataset/2pye/2pye_tcr.pdb -o example/output -op 2pye -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/3dxa/3dxa_pmhc_renumbered.pdb -l example/input/benchmark_dataset/3dxa/3dxa_tcr.pdb -o example/output -op 3dxa -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/3h9s/3h9s_pmhc_renumbered.pdb -l example/input/benchmark_dataset/3h9s/3h9s_tcr.pdb -o example/output -op 3h9s -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/3kpr/3kpr_pmhc_renumbered.pdb -l example/input/benchmark_dataset/3kpr/3kpr_tcr.pdb -o example/output -op 3kpr -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/3kps/3kps_pmhc_renumbered.pdb -l example/input/benchmark_dataset/3kps/3kps_tcr.pdb -o example/output -op 3kps -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/3pwp/3pwp_pmhc_renumbered.pdb -l example/input/benchmark_dataset/3pwp/3pwp_tcr.pdb -o example/output -op 3pwp -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/3qdg/3qdg_pmhc_renumbered.pdb -l example/input/benchmark_dataset/3qdg/3qdg_tcr.pdb -o example/output -op 3qdg -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/3qdj/3qdj_pmhc_renumbered.pdb -l example/input/benchmark_dataset/3qdj/3qdj_tcr.pdb -o example/output -op 3qdj -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/3sjv/3sjv_pmhc_renumbered.pdb -l example/input/benchmark_dataset/3sjv/3sjv_tcr.pdb -o example/output -op 3sjv -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/3utt/3utt_pmhc_renumbered.pdb -l example/input/benchmark_dataset/3utt/3utt_tcr.pdb -o example/output -op 3utt -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/3vxr/3vxr_pmhc_renumbered.pdb -l example/input/benchmark_dataset/3vxr/3vxr_tcr.pdb -o example/output -op 3vxr -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/3vxs/3vxs_pmhc_renumbered.pdb -l example/input/benchmark_dataset/3vxs/3vxs_tcr.pdb -o example/output -op 3vxs -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/3w0w/3w0w_pmhc_renumbered.pdb -l example/input/benchmark_dataset/3w0w/3w0w_tcr.pdb -o example/output -op 3w0w -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/4jfd/4jfd_pmhc_renumbered.pdb -l example/input/benchmark_dataset/4jfd/4jfd_tcr.pdb -o example/output -op 4jfd -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/4jff/4jff_pmhc_renumbered.pdb -l example/input/benchmark_dataset/4jff/4jff_tcr.pdb -o example/output -op 4jff -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/5c0a/5c0a_pmhc_renumbered.pdb -l example/input/benchmark_dataset/5c0a/5c0a_tcr.pdb -o example/output -op 5c0a -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/5c0b/5c0b_pmhc_renumbered.pdb -l example/input/benchmark_dataset/5c0b/5c0b_tcr.pdb -o example/output -op 5c0b -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/5c0c/5c0c_pmhc_renumbered.pdb -l example/input/benchmark_dataset/5c0c/5c0c_tcr.pdb -o example/output -op 5c0c -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/5c07/5c07_pmhc_renumbered.pdb -l example/input/benchmark_dataset/5c07/5c07_tcr.pdb -o example/output -op 5c07 -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/5c08/5c08_pmhc_renumbered.pdb -l example/input/benchmark_dataset/5c08/5c08_tcr.pdb -o example/output -op 5c08 -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/5c09/5c09_pmhc_renumbered.pdb -l example/input/benchmark_dataset/5c09/5c09_tcr.pdb -o example/output -op 5c09 -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/5hhm/5hhm_pmhc_renumbered.pdb -l example/input/benchmark_dataset/5hhm/5hhm_tcr.pdb -o example/output -op 5hhm -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/5hyj/5hyj_pmhc_renumbered.pdb -l example/input/benchmark_dataset/5hyj/5hyj_tcr.pdb -o example/output -op 5hyj -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/5ivx/5ivx_pmhc_renumbered.pdb -l example/input/benchmark_dataset/5ivx/5ivx_tcr.pdb -o example/output -op 5ivx -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/5nme/5nme_pmhc_renumbered.pdb -l example/input/benchmark_dataset/5nme/5nme_tcr.pdb -o example/output -op 5nme -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/5nmf/5nmf_pmhc_renumbered.pdb -l example/input/benchmark_dataset/5nmf/5nmf_tcr.pdb -o example/output -op 5nmf -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/5nmg/5nmg_pmhc_renumbered.pdb -l example/input/benchmark_dataset/5nmg/5nmg_tcr.pdb -o example/output -op 5nmg -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/6amu/6amu_pmhc_renumbered.pdb -l example/input/benchmark_dataset/6amu/6amu_tcr.pdb -o example/output -op 6amu -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/6avf/6avf_pmhc_renumbered.pdb -l example/input/benchmark_dataset/6avf/6avf_tcr.pdb -o example/output -op 6avf -c 12 -t 3

wait

python3 scripts/swift_tcr.py -r example/input/benchmark_dataset/6eqb/6eqb_pmhc_renumbered.pdb -l example/input/benchmark_dataset/6eqb/6eqb_tcr.pdb -o example/output -op 6eqb -c 12 -t 3